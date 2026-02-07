package instance;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

/**
 * Parse MVPRP instance txt (delivery-style) and map it to your pickup+production model parameters.
 *
 * =========================
 * Mapping / Assumptions
 * =========================
 * 1) Pickup-generation mapping:
 *    - For each pickup node i, set s_{it} := Demand_i for all t in T (1..l)
 *    - And set s_{i,l+1} := 0 to support dummy period in T' = T ∪ {l+1}
 *
 * 2) External finished-goods demand d_t:
 *    - Default: d_t := (sum_i s_{it}) / k
 *
 * 3) Parameters NOT provided by txt (set via Options):
 *    - k, P00, hp, Lp, bigM
 *    - distance rounding mode for c_ij
 *
 * Node indexing:
 *   0 = depot/factory
 *   1..n = pickup nodes
 *   n+1 = return depot (copy of node 0)
 */
public final class Instance {

    // -----------------------
    // Options / Defaults
    // -----------------------
    public static final class Options {
        public enum DistanceMode { EUCLIDEAN_FLOAT, EUCLIDEAN_ROUND, EUCLIDEAN_FLOOR }

        /** raw->finished conversion factor k */
        public double k = 1.0;

        /** finished-goods initial inventory P00 */
        public double P00 = 0.0;

        /** finished-goods holding cost hp */
        public double hp = 0.0;

        /** finished-goods capacity Lp (if <=0, default to L0) */
        public double Lp = -1.0;

        /** Big-M for MTZ (if <=0, default to Q) */
        public double bigM = -1.0;

        /** how to build c_ij from coordinates */
        public DistanceMode distanceMode = DistanceMode.EUCLIDEAN_FLOAT;

        /** If true, d_t := sum_i s_{it}/k; else d_t := 0 */
        public boolean autoSetDt = true;

        public static Options defaults() { return new Options(); }
    }

    // -----------------------
    // Basic instance fields
    // -----------------------
    public final int n;          // # pickup nodes
    public final int l;          // # periods
    public final int K;          // max vehicles per period
    public final double Q;       // vehicle capacity

    // Nodes: 0 = depot/factory, 1..n = pickup nodes, n+1 = return depot
    public final int nodeCount;  // = n + 2
    public final double[] x;     // size nodeCount
    public final double[] y;     // size nodeCount

    // Factory/supplier (node 0) params from txt
    public final double I00;     // factory raw initial inventory
    public final double L0;      // factory raw capacity
    public final double C;       // production capacity per period
    public final double h0;      // factory raw holding cost
    public final double u;       // unit production cost
    public final double f;       // setup cost

    // Pickup nodes i=1..n params from txt
    public final double[] Ii0;   // initial inventory at pickup node i
    public final double[] Li;    // capacity at pickup node i
    public final double[] hi;    // holding cost at pickup node i
    public final double[] demand;// Demand_i from txt

    // -----------------------
    // Your model extra params (some auto-set)
    // -----------------------
    public final double k;       // conversion factor
    public final double P00;     // finished initial inventory
    public final double hp;      // finished holding cost
    public final double Lp;      // finished capacity
    public final double bigM;    // Big-M

    // s_{it}: pickup generation. index: i=1..n, t=1..l+1
    public final double[][] s;        // [n+1][l+2], ignore index 0

    // external finished-goods demand d_t. index t=1..l
    public final double[] dt;         // [l+1], ignore index 0

    // Distance/cost matrix on V={0..n+1}
    public final double[][] c;        // [nodeCount][nodeCount]

    // prefixS[i][t] = sum_{r=1..t} s_{ir}
    private final double[][] prefixS;     // [n+1][l+2]
    // cumPrefix[i][t] = sum_{j=0..t} prefixS[i][j]
    private final double[][] cumPrefix;   // [n+1][l+2]

    // pi(i,t) for i=1..n, t=1..l+1
    public final int[][] pi;          // [n+1][l+2]
    // mu(i,t) for i=1..n, t=0..l
    public final int[][] mu;          // [n+1][l+1]

    private Instance(Parsed parsed, Options opt) {
        this.n = parsed.n;
        this.l = parsed.l;
        this.K = parsed.K;
        this.Q = parsed.Q;

        this.nodeCount = n + 2;
        this.x = new double[nodeCount];
        this.y = new double[nodeCount];

        // Supplier node 0
        this.x[0] = parsed.supplierX;
        this.y[0] = parsed.supplierY;

        // Copy depot to node n+1
        this.x[n + 1] = this.x[0];
        this.y[n + 1] = this.y[0];

        this.I00 = parsed.supplierInit;
        this.L0  = parsed.supplierMax;
        this.C   = parsed.prodCap;
        this.h0  = parsed.supplierHold;
        this.u   = parsed.varCost;
        this.f   = parsed.fixCost;

        // Pickup nodes
        this.Ii0 = new double[n + 1];
        this.Li  = new double[n + 1];
        this.hi  = new double[n + 1];
        this.demand = new double[n + 1];

        for (int i = 1; i <= n; i++) {
            Node nd = parsed.retailers.get(i);
            if (nd == null) {
                throw new IllegalArgumentException("Missing Retailer " + i + " in txt.");
            }
            this.x[i] = nd.x;
            this.y[i] = nd.y;
            this.Ii0[i] = nd.init;
            this.Li[i]  = nd.max;
            this.hi[i]  = nd.hold;
            this.demand[i] = nd.demand;
        }

        // Options
        this.k = opt.k;
        this.P00 = opt.P00;
        this.hp = opt.hp;
        this.Lp = (opt.Lp > 0) ? opt.Lp : this.L0;
        this.bigM = (opt.bigM > 0) ? opt.bigM : this.Q;

        // Build s_{it} from demand, and s_{i,l+1}=0
        this.s = new double[n + 1][l + 2];
        for (int i = 1; i <= n; i++) {
            for (int t = 1; t <= l; t++) {
                this.s[i][t] = this.demand[i];
            }
            this.s[i][l + 1] = 0.0;
        }

        // Build dt
        this.dt = new double[l + 1];
        if (opt.autoSetDt) {
            double sumPerPeriod = 0.0;
            for (int i = 1; i <= n; i++) sumPerPeriod += this.s[i][1];
            double per = sumPerPeriod / this.k;
            for (int t = 1; t <= l; t++) this.dt[t] = per;
        } else {
            for (int t = 1; t <= l; t++) this.dt[t] = 0.0;
        }

        // Build cost matrix c_ij
        this.c = new double[nodeCount][nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            for (int j = 0; j < nodeCount; j++) {
                if (i == j) { c[i][j] = 0.0; continue; }
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dist = Math.sqrt(dx * dx + dy * dy);
                c[i][j] = applyDistanceMode(dist, opt.distanceMode);
            }
        }

        // Prefix sums
        this.prefixS = new double[n + 1][l + 2];
        this.cumPrefix = new double[n + 1][l + 2];

        for (int i = 1; i <= n; i++) {
            prefixS[i][0] = 0.0;
            for (int t = 1; t <= l + 1; t++) {
                prefixS[i][t] = prefixS[i][t - 1] + s[i][t];
            }
            cumPrefix[i][0] = prefixS[i][0]; // 0
            for (int t = 1; t <= l + 1; t++) {
                cumPrefix[i][t] = cumPrefix[i][t - 1] + prefixS[i][t];
            }
        }

        // Compute pi and mu
        this.pi = new int[n + 1][l + 2];
        this.mu = new int[n + 1][l + 1];

        // pi(i,t): t in 1..l+1
        for (int i = 1; i <= n; i++) {
            for (int t = 1; t <= l + 1; t++) {
                int best = -1;
                for (int v = 0; v <= t - 1; v++) {
                    double g = g(i, v, t);
                    if (g <= Li[i] + 1e-9) { best = v; break; }
                }
                if (best < 0) {
                    throw new IllegalStateException("No feasible pi for node " + i + " at t=" + t);
                }
                pi[i][t] = best;
            }
        }

        // mu(i,t): t in 0..l
        for (int i = 1; i <= n; i++) {
            for (int t = 0; t <= l; t++) {
                int best = -1;
                for (int v = l + 1; v >= t + 1; v--) {
                    double g = g(i, t, v);
                    if (g <= Li[i] + 1e-9) { best = v; break; }
                }
                if (best < 0) {
                    throw new IllegalStateException("No feasible mu for node " + i + " at t=" + t);
                }
                mu[i][t] = best;
            }
        }
    }

    // -----------------------
    // Public factory methods
    // -----------------------
    public static Instance fromFile(String filePath) throws IOException {
        return fromFile(Paths.get(filePath), Options.defaults());
    }

    public static Instance fromFile(String filePath, Options opt) throws IOException {
        return fromFile(Paths.get(filePath), opt);
    }

    public static Instance fromFile(Path path) throws IOException {
        return fromFile(path, Options.defaults());
    }

    public static Instance fromFile(Path path, Options opt) throws IOException {
        BufferedReader br = Files.newBufferedReader(path, StandardCharsets.UTF_8);
        try {
            Parsed p = parse(br);
            return new Instance(p, opt);
        } finally {
            br.close();
        }
    }

    public static Instance fromString(String txtContent) throws IOException {
        return fromString(txtContent, Options.defaults());
    }

    public static Instance fromString(String txtContent, Options opt) throws IOException {
        BufferedReader br = new BufferedReader(new StringReader(txtContent));
        Parsed p = parse(br);
        return new Instance(p, opt);
    }

    // -----------------------
    // Core model helper functions
    // -----------------------

    /** g_{ivt}: pickup amount at node i when visited at t and previous visit at v (v<t, t in 1..l+1) */
    public double g(int i, int v, int t) {
        if (i < 1 || i > n) throw new IllegalArgumentException("i must be in 1..n");
        if (t < 1 || t > l + 1) throw new IllegalArgumentException("t must be in 1..l+1");
        if (v < 0 || v >= t) throw new IllegalArgumentException("v must satisfy 0<=v<t");

        if (v == 0) {
            // g_{i0t} = I_{i0} + sum_{j=1..t} s_{ij}
            return Ii0[i] + prefixS[i][t];
        } else {
            // g_{ivt} = sum_{j=v+1..t} s_{ij}
            return prefixS[i][t] - prefixS[i][v];
        }
    }

    /** e_{ivt}: holding cost accumulated from period v to t-1 (your definition). */
    public double e(int i, int v, int t) {
        if (i < 1 || i > n) throw new IllegalArgumentException("i must be in 1..n");
        if (t < 1 || t > l + 1) throw new IllegalArgumentException("t must be in 1..l+1");
        if (v < 0 || v >= t) throw new IllegalArgumentException("v must satisfy 0<=v<t");

        double h = this.hi[i];

        if (v == 0) {
            // e_{i0t} = h_i * sum_{j=0..t-1} (I_{i0} + sum_{r=1..j} s_ir)
            //         = h_i * [ t*I_{i0} + sum_{j=0..t-1} prefixS[i][j] ]
            double sumPrefix0To = cumPrefix[i][t - 1];
            return h * (t * Ii0[i] + sumPrefix0To);
        } else {
            // e_{ivt} = h_i * sum_{j=v+1..t-1} (sum_{r=v+1..j} s_ir)
            //         = h_i * [ sum_{j=v+1..t-1} prefixS[i][j] - (t-v-1)*prefixS[i][v] ]
            if (t - 1 < v + 1) return 0.0;
            double sumPrefix = cumPrefix[i][t - 1] - cumPrefix[i][v];
            double count = (t - v - 1);
            return h * (sumPrefix - count * prefixS[i][v]);
        }
    }

    /** feasible previous v range: v in [pi(i,t), t-1] */
    public int[] feasibleVs(int i, int t) {
        if (i < 1 || i > n) throw new IllegalArgumentException("i must be in 1..n");
        if (t < 1 || t > l + 1) throw new IllegalArgumentException("t must be in 1..l+1");
        int start = pi[i][t];
        int end = t - 1;
        int[] vs = new int[end - start + 1];
        for (int idx = 0; idx < vs.length; idx++) vs[idx] = start + idx;
        return vs;
    }

    public double cost(int i, int j) {
        return c[i][j];
    }

    private static double applyDistanceMode(double dist, Options.DistanceMode mode) {
        if (mode == Options.DistanceMode.EUCLIDEAN_FLOAT) return dist;
        if (mode == Options.DistanceMode.EUCLIDEAN_FLOOR) return Math.floor(dist);
        // EUCLIDEAN_ROUND (use Math.rint for nearest int in double)
        return Math.rint(dist);
    }

    // -----------------------
    // Parsing
    // -----------------------

    private static final class Node {
        final double x, y;
        final double init, max, demand, hold;
        Node(double x, double y, double init, double max, double demand, double hold) {
            this.x = x; this.y = y;
            this.init = init; this.max = max;
            this.demand = demand; this.hold = hold;
        }
    }

    private static final class Parsed {
        int n = -1;
        int l = -1;
        int K = -1;
        double Q = Double.NaN;

        // Supplier
        double supplierX, supplierY;
        double supplierInit, supplierMax, prodCap, supplierHold, varCost, fixCost;

        // Retailers
        Map<Integer, Node> retailers = new HashMap<Integer, Node>();
    }

    private static Parsed parse(BufferedReader br) throws IOException {
        Parsed p = new Parsed();
        String[] retailerHeader = null;

        String line;
        while ((line = br.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;
            if (line.startsWith("#") || line.startsWith("//")) continue;

            String[] tok = line.split("\\s+");
            if (tok.length < 2) continue;

            String first = normalizeKey(tok[0]);

            if ("PeriodNumber".equalsIgnoreCase(first)) {
                p.l = (int) parseDouble(tok[1]);
            } else if ("CustomerNumber".equalsIgnoreCase(first)) {
                p.n = (int) parseDouble(tok[1]);
            } else if ("VechileNumber".equalsIgnoreCase(first) || "VehicleNumber".equalsIgnoreCase(first)) {
                p.K = (int) parseDouble(tok[1]);
            } else if ("VehicleCapacity".equalsIgnoreCase(first)) {
                p.Q = parseDouble(tok[1]);
            } else if ("Supplier".equalsIgnoreCase(first)) {
                Map<String, String> kv;
                if (tok.length >= 2 && !isNumericToken(tok[1])) {
                    // Table style:
                    // Supplier COORDX COORDY ...
                    // 0        ...
                    String[] header = subArray(tok, 1);
                    String dataLine = readNextDataLine(br);
                    if (dataLine == null) {
                        throw new IllegalArgumentException("Missing supplier data row after Supplier header.");
                    }
                    String[] row = dataLine.split("\\s+");
                    kv = parseHeaderRow(header, row);
                } else {
                    // Key-value style:
                    // Supplier 0 COORDX ... COORDY ...
                    kv = parseKeyValues(tok, 2);
                }
                fillSupplier(p, kv);
            } else if ("Retailer".equalsIgnoreCase(first)) {
                if (tok.length >= 2 && !isNumericToken(tok[1])) {
                    // Table header style:
                    // Retailer COORDX COORDY InitLevel MaxLevel Demand HoldCost
                    retailerHeader = subArray(tok, 1);
                } else {
                    // Key-value style:
                    // Retailer 1 COORDX ... COORDY ...
                    int idx = (int) parseDouble(tok[1]);
                    Map<String, String> kv = parseKeyValues(tok, 2);
                    addRetailer(p, idx, kv);
                }
            } else if (retailerHeader != null && isNumericToken(first)) {
                // Table data style:
                // 1 x y init max demand hold
                int idx = (int) parseDouble(first);
                Map<String, String> kv = parseHeaderRow(retailerHeader, tok);
                addRetailer(p, idx, kv);
            }
        }

        // Validate
        if (p.n <= 0) throw new IllegalArgumentException("CustomerNumber not found or invalid.");
        if (p.l <= 0) throw new IllegalArgumentException("PeriodNumber not found or invalid.");
        if (p.K <= 0) throw new IllegalArgumentException("VechileNumber/VehicleNumber not found or invalid.");
        if (!(p.Q > 0)) throw new IllegalArgumentException("VehicleCapacity not found or invalid.");

        if (p.retailers.size() < p.n) {
            throw new IllegalArgumentException("Retailer lines are incomplete: expected " + p.n + ", got " + p.retailers.size());
        }

        return p;
    }

    private static void fillSupplier(Parsed p, Map<String, String> kv) {
        p.supplierX = parseDoubleRequired(kv, "COORDX");
        p.supplierY = parseDoubleRequired(kv, "COORDY");
        p.supplierInit = parseDoubleRequired(kv, "InitLevel");
        p.supplierMax  = parseDoubleRequired(kv, "MaxLevel");
        p.prodCap      = parseDoubleRequired(kv, "ProdCapacity");
        p.supplierHold = parseDoubleRequired(kv, "HoldCost");
        p.varCost      = parseDoubleRequired(kv, "VarCost");
        p.fixCost      = parseDoubleRequired(kv, "FixCost");
    }

    private static void addRetailer(Parsed p, int idx, Map<String, String> kv) {
        double x = parseDoubleRequired(kv, "COORDX");
        double y = parseDoubleRequired(kv, "COORDY");
        double init = parseDoubleRequired(kv, "InitLevel");
        double max  = parseDoubleRequired(kv, "MaxLevel");
        double dem  = parseDoubleRequired(kv, "Demand");
        double hold = parseDoubleRequired(kv, "HoldCost");
        p.retailers.put(idx, new Node(x, y, init, max, dem, hold));
    }

    private static Map<String, String> parseHeaderRow(String[] header, String[] row) {
        int offset;
        if (row.length == header.length + 1) {
            offset = 1; // the first token is typically Supplier/Retailer index
        } else if (row.length == header.length) {
            offset = 0;
        } else {
            throw new IllegalArgumentException("Row length does not match header. header=" + header.length + ", row=" + row.length);
        }

        Map<String, String> kv = new HashMap<String, String>();
        for (int i = 0; i < header.length; i++) {
            kv.put(header[i], row[i + offset]);
        }
        return kv;
    }

    private static String[] subArray(String[] arr, int startInclusive) {
        String[] out = new String[arr.length - startInclusive];
        for (int i = startInclusive; i < arr.length; i++) {
            out[i - startInclusive] = arr[i];
        }
        return out;
    }

    private static boolean isNumericToken(String s) {
        try {
            Double.parseDouble(s);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    private static String normalizeKey(String key) {
        if (key == null) {
            return "";
        }
        return key.endsWith(":") ? key.substring(0, key.length() - 1) : key;
    }

    private static String readNextDataLine(BufferedReader br) throws IOException {
        String line;
        while ((line = br.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;
            if (line.startsWith("#") || line.startsWith("//")) continue;
            return line;
        }
        return null;
    }

    private static Map<String, String> parseKeyValues(String[] tok, int startIdx) {
        Map<String, String> kv = new HashMap<String, String>();
        int i = startIdx;
        while (i + 1 < tok.length) {
            kv.put(tok[i], tok[i + 1]);
            i += 2;
        }
        return kv;
    }

    private static double parseDoubleRequired(Map<String, String> kv, String key) {
        String v = kv.get(key);
        if (v == null) throw new IllegalArgumentException("Missing key: " + key);
        return parseDouble(v);
    }

    private static double parseDouble(String s) {
        return Double.parseDouble(s);
    }

    // -----------------------
    // Optional demo
    // -----------------------
    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            System.out.println("Usage: java Instance <path_to_instance_txt>");
            return;
        }

        Options opt = Options.defaults();
        opt.distanceMode = Options.DistanceMode.EUCLIDEAN_FLOAT;
        opt.k = 1.0;
        opt.autoSetDt = true;

        Instance data = Instance.fromFile(args[0], opt);

        System.out.println("n=" + data.n + ", l=" + data.l + ", K=" + data.K + ", Q=" + data.Q);
        System.out.println("Depot: (" + data.x[0] + "," + data.y[0] + "), L0=" + data.L0 + ", C=" + data.C);
        System.out.println("d_t(t=1)=" + data.dt[1]);
        System.out.println("Example: node 1 pi(t=3)=" + data.pi[1][3] + ", mu(t=0)=" + data.mu[1][0]);
        System.out.println("Example: g(1,0,2)=" + data.g(1,0,2) + ", e(1,0,2)=" + data.e(1,0,2));
        System.out.println("Cost c(0,1)=" + data.cost(0,1));
    }
}
