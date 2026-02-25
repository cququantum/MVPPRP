package lbbdModel;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import instance.Instance;
import model.CplexConfig;

import java.util.Arrays;
import java.util.HashMap;

public final class PeriodCvrpSubproblemSolver implements AutoCloseable {
    private static final boolean LOG_TO_CONSOLE = false;
    private static final double EPS = 1e-9;
    private static final int EXACT_DP_MAX_N = 18;

    private final HashMap<Integer, ExactPeriodModel> exactModels = new HashMap<Integer, ExactPeriodModel>();
    private final HashMap<Integer, LpPeriodModel> lpModels = new HashMap<Integer, LpPeriodModel>();
    private Instance subsetDpInstance = null;
    private ExactSubsetDpOracle subsetDpOracle = null;

    public static final class Result {
        public final boolean feasible;
        public final boolean optimal;
        public final String status;
        public final double objective;

        public Result(boolean feasible, boolean optimal, String status, double objective) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.status = status;
            this.objective = objective;
        }
    }

    public static final class DualCutResult {
        public final boolean feasible;
        public final boolean optimal;
        public final String status;
        public final double lpObjective;
        public final double constant;
        public final double[] coeffZ;
        public final double[] coeffQ;

        public DualCutResult(
                boolean feasible,
                boolean optimal,
                String status,
                double lpObjective,
                double constant,
                double[] coeffZ,
                double[] coeffQ
        ) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.status = status;
            this.lpObjective = lpObjective;
            this.constant = constant;
            this.coeffZ = coeffZ;
            this.coeffQ = coeffQ;
        }
    }

    public Result solve(Instance ins, int t, double[] qBar, int[] zBar) {
        validateInput(ins, t, qBar, zBar);

        double totalPickup = totalPickup(ins, qBar);
        int visitCount = visitCount(ins, zBar);
        if (visitCount == 0 && totalPickup <= EPS) {
            return new Result(true, true, "TrivialZero", 0.0);
        }
        if (totalPickup > ins.Q * ins.K + EPS) {
            return new Result(false, false, "InfeasibleByTotalCapacity", Double.NaN);
        }

        Result dpResult = trySolveBySubsetDp(ins, qBar, zBar);
        if (dpResult != null) {
            return dpResult;
        }

        try {
            ExactPeriodModel model = exactModels.get(t);
            if (model == null) {
                model = new ExactPeriodModel(ins, t);
                exactModels.put(t, model);
            }
            model.setSolveThreads(recommendedExactThreads(visitCount));
            return model.solve(qBar, zBar, totalPickup);
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve CVRP subproblem at t=" + t, e);
        }
    }

    public DualCutResult solveLpDualCut(Instance ins, int t, double[] qBar, int[] zBar) {
        validateInput(ins, t, qBar, zBar);

        double totalPickup = totalPickup(ins, qBar);
        int visitCount = visitCount(ins, zBar);
        if (visitCount == 0 && totalPickup <= EPS) {
            return new DualCutResult(
                    true,
                    true,
                    "TrivialZeroLP",
                    0.0,
                    0.0,
                    new double[ins.n + 1],
                    new double[ins.n + 1]
            );
        }
        if (totalPickup > ins.Q * ins.K + EPS) {
            return new DualCutResult(false, false, "InfeasibleByTotalCapacity", Double.NaN, Double.NaN, null, null);
        }

        try {
            LpPeriodModel model = lpModels.get(t);
            if (model == null) {
                model = new LpPeriodModel(ins, t);
                lpModels.put(t, model);
            }
            return model.solve(ins, qBar, zBar, totalPickup);
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve LP CVRP dual-cut subproblem at t=" + t, e);
        }
    }

    @Override
    public void close() {
        for (ExactPeriodModel model : exactModels.values()) {
            model.close();
        }
        exactModels.clear();
        for (LpPeriodModel model : lpModels.values()) {
            model.close();
        }
        lpModels.clear();
    }

    private static void validateInput(Instance ins, int t, double[] qBar, int[] zBar) {
        if (t < 1 || t > ins.l) {
            throw new IllegalArgumentException("t must be in 1..l");
        }
        if (qBar == null || qBar.length < ins.n + 1) {
            throw new IllegalArgumentException("qBar length must be at least n+1");
        }
        if (zBar == null || zBar.length < ins.n + 1) {
            throw new IllegalArgumentException("zBar length must be at least n+1");
        }
    }

    private static double totalPickup(Instance ins, double[] qBar) {
        double total = 0.0;
        for (int i = 1; i <= ins.n; i++) {
            total += qBar[i];
        }
        return total;
    }

    private static int visitCount(Instance ins, int[] zBar) {
        int cnt = 0;
        for (int i = 1; i <= ins.n; i++) {
            if (zBar[i] != 0) {
                cnt++;
            }
        }
        return cnt;
    }

    private static boolean isPickupNode(int node, int n) {
        return node >= 1 && node <= n;
    }

    private Result trySolveBySubsetDp(Instance ins, double[] qBar, int[] zBar) {
        if (ins.n > EXACT_DP_MAX_N) {
            return null;
        }
        if (subsetDpOracle == null || subsetDpInstance != ins) {
            subsetDpOracle = new ExactSubsetDpOracle(ins);
            subsetDpInstance = ins;
        }
        return subsetDpOracle.solve(qBar, zBar);
    }

    private static int recommendedExactThreads(int visitCount) {
        if (visitCount <= 6) {
            return 1;
        }
        if (visitCount <= 10) {
            return 2;
        }
        if (visitCount <= 13) {
            return 3;
        }
        return 4;
    }

    private static void configure(IloCplex cplex) throws IloException {
        cplex.setParam(IloCplex.Param.TimeLimit, CplexConfig.TIME_LIMIT_SEC);
        cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, CplexConfig.MIP_GAP);
        cplex.setParam(IloCplex.Param.Threads, 1);
        if (!LOG_TO_CONSOLE) {
            cplex.setOut(null);
            cplex.setWarning(null);
        }
    }

    /**
     * Exact CVRP cost for fixed-period demand via:
     * 1) global Held-Karp route-cost precomputation over all customer subsets (distance only)
     * 2) per-call partition DP with capacity and at-most-K routes.
     *
     * This is exact and much faster than solving a MIP repeatedly for n <= 15.
     */
    private static final class ExactSubsetDpOracle {
        private static final double INF = 1e100;

        private final Instance ins;
        private final int n;
        private final int maskCount;
        private final int fullMask;
        private final double[] routeCostByMask;
        private final double[] subsetLoad;
        private final double[][] memoByK;
        private final boolean[][] seenByK;

        ExactSubsetDpOracle(Instance ins) {
            this.ins = ins;
            this.n = ins.n;
            this.maskCount = 1 << n;
            this.fullMask = maskCount - 1;
            this.routeCostByMask = precomputeRouteCosts(ins, n, maskCount);
            this.subsetLoad = new double[maskCount];
            this.memoByK = new double[ins.K + 1][maskCount];
            this.seenByK = new boolean[ins.K + 1][maskCount];
        }

        Result solve(double[] qBar, int[] zBar) {
            int activeMask = 0;
            for (int i = 1; i <= n; i++) {
                if (zBar[i] != 0) {
                    activeMask |= (1 << (i - 1));
                } else if (qBar[i] > EPS) {
                    return new Result(false, false, "InfeasibleByPositivePickupWithoutVisit", Double.NaN);
                }
            }
            if (activeMask == 0) {
                return new Result(true, true, "TrivialZeroDP", 0.0);
            }

            buildSubsetLoads(qBar);

            // If any visited customer alone exceeds capacity, infeasible.
            int bits = activeMask;
            while (bits != 0) {
                int bit = bits & -bits;
                int idx = Integer.numberOfTrailingZeros(bit);
                if (qBar[idx + 1] > ins.Q + EPS) {
                    return new Result(false, false, "InfeasibleBySingleCustomerCapacity", Double.NaN);
                }
                bits ^= bit;
            }

            clearMemo();
            double best = solvePartition(activeMask, ins.K);
            if (best >= INF / 2) {
                return new Result(false, false, "InfeasibleByRoutePartitionDP", Double.NaN);
            }
            return new Result(true, true, "OptimalDP", best);
        }

        private void buildSubsetLoads(double[] qBar) {
            subsetLoad[0] = 0.0;
            for (int mask = 1; mask <= fullMask; mask++) {
                int bit = mask & -mask;
                int idx = Integer.numberOfTrailingZeros(bit); // 0-based customer local index
                subsetLoad[mask] = subsetLoad[mask ^ bit] + qBar[idx + 1];
            }
        }

        private void clearMemo() {
            for (int k = 0; k <= ins.K; k++) {
                Arrays.fill(seenByK[k], false);
            }
        }

        private double solvePartition(int mask, int kLeft) {
            if (mask == 0) {
                return 0.0;
            }
            if (kLeft <= 0) {
                return INF;
            }
            if (subsetLoad[mask] > kLeft * ins.Q + EPS) {
                return INF;
            }
            if (kLeft == 1) {
                return (subsetLoad[mask] <= ins.Q + EPS) ? routeCostByMask[mask] : INF;
            }
            if (seenByK[kLeft][mask]) {
                return memoByK[kLeft][mask];
            }
            seenByK[kLeft][mask] = true;

            double best = INF;
            int anchor = mask & -mask; // symmetry breaking: first route must contain anchor
            int sub = mask;
            while (sub != 0) {
                if ((sub & anchor) != 0 && subsetLoad[sub] <= ins.Q + EPS) {
                    double rest = solvePartition(mask ^ sub, kLeft - 1);
                    if (rest < INF / 2) {
                        double val = routeCostByMask[sub] + rest;
                        if (val < best) {
                            best = val;
                        }
                    }
                }
                sub = (sub - 1) & mask;
            }

            memoByK[kLeft][mask] = best;
            return best;
        }

        private static double[] precomputeRouteCosts(Instance ins, int n, int maskCount) {
            double[] routeCost = new double[maskCount];
            double[][] pathToLast = new double[maskCount][n];
            for (int mask = 0; mask < maskCount; mask++) {
                Arrays.fill(pathToLast[mask], INF);
            }
            routeCost[0] = 0.0;

            for (int j = 0; j < n; j++) {
                int mask = 1 << j;
                int gj = j + 1;
                pathToLast[mask][j] = ins.c[0][gj];
                routeCost[mask] = ins.c[0][gj] + ins.c[gj][ins.n + 1];
            }

            for (int mask = 1; mask < maskCount; mask++) {
                if ((mask & (mask - 1)) == 0) {
                    continue; // singleton already handled
                }

                int bits = mask;
                while (bits != 0) {
                    int bitJ = bits & -bits;
                    int j = Integer.numberOfTrailingZeros(bitJ);
                    int prevMask = mask ^ bitJ;
                    int gj = j + 1;

                    double best = INF;
                    int prevBits = prevMask;
                    while (prevBits != 0) {
                        int bitI = prevBits & -prevBits;
                        int i = Integer.numberOfTrailingZeros(bitI);
                        double prev = pathToLast[prevMask][i];
                        if (prev < INF / 2) {
                            double cand = prev + ins.c[i + 1][gj];
                            if (cand < best) {
                                best = cand;
                            }
                        }
                        prevBits ^= bitI;
                    }
                    pathToLast[mask][j] = best;
                    bits ^= bitJ;
                }

                double bestRoute = INF;
                int endBits = mask;
                while (endBits != 0) {
                    int bitEnd = endBits & -endBits;
                    int last = Integer.numberOfTrailingZeros(bitEnd);
                    double path = pathToLast[mask][last];
                    if (path < INF / 2) {
                        double cand = path + ins.c[last + 1][ins.n + 1];
                        if (cand < bestRoute) {
                            bestRoute = cand;
                        }
                    }
                    endBits ^= bitEnd;
                }
                routeCost[mask] = bestRoute;
            }
            return routeCost;
        }
    }

    private static final class ExactPeriodModel {
        private final Instance ins;
        private final int t;
        private final int n;
        private final int nodeCount;
        private final IloCplex cplex;
        private final IloRange vehCap;
        private final IloRange[] visitLink;
        private final IloRange[][] mtzRows;

        ExactPeriodModel(Instance ins, int t) throws IloException {
            this.ins = ins;
            this.t = t;
            this.n = ins.n;
            this.nodeCount = ins.nodeCount;
            this.cplex = new IloCplex();
            configure(cplex);

            IloNumVar m = cplex.intVar(0, ins.K, "m_" + t);
            IloNumVar[] u = new IloNumVar[nodeCount];
            IloNumVar[][] x = new IloNumVar[nodeCount][nodeCount];
            visitLink = new IloRange[n + 1];
            mtzRows = new IloRange[nodeCount][nodeCount];

            for (int i = 0; i < nodeCount; i++) {
                u[i] = cplex.numVar(0.0, Double.MAX_VALUE, "u_" + i + "_" + t);
            }
            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    x[i][j] = cplex.boolVar("x_" + i + "_" + j + "_" + t);
                }
            }

            IloLinearNumExpr obj = cplex.linearNumExpr();
            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    obj.addTerm(ins.c[i][j], x[i][j]);
                }
            }
            cplex.addMinimize(obj);

            IloLinearNumExpr vehCapExpr = cplex.linearNumExpr();
            vehCapExpr.addTerm(-ins.Q, m);
            vehCap = cplex.addLe(vehCapExpr, 0.0, "VehicleCap_" + t);

            for (int i = 1; i <= n; i++) {
                IloLinearNumExpr out = cplex.linearNumExpr();
                IloLinearNumExpr in = cplex.linearNumExpr();
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    out.addTerm(1.0, x[i][j]);
                    in.addTerm(1.0, x[j][i]);
                }
                cplex.addEq(out, in, "FlowBalance_" + i + "_" + t);
            }

            IloLinearNumExpr depart = cplex.linearNumExpr();
            for (int j = 1; j <= n; j++) {
                depart.addTerm(1.0, x[0][j]);
            }
            cplex.addEq(depart, m, "DepartFactory_" + t);

            IloLinearNumExpr back = cplex.linearNumExpr();
            for (int i = 1; i <= n; i++) {
                back.addTerm(1.0, x[i][n + 1]);
            }
            cplex.addEq(back, m, "ReturnFactory_" + t);

            for (int j = 1; j <= n; j++) {
                IloLinearNumExpr incoming = cplex.linearNumExpr();
                for (int i = 0; i < nodeCount; i++) {
                    if (i == j) {
                        continue;
                    }
                    incoming.addTerm(1.0, x[i][j]);
                }
                visitLink[j] = cplex.addEq(incoming, 0.0, "VisitLink_" + j + "_" + t);
            }
            cplex.addLe(m, ins.K, "VehicleCount_" + t);

            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    IloLinearNumExpr mtz = cplex.linearNumExpr();
                    mtz.addTerm(1.0, u[j]);
                    mtz.addTerm(-1.0, u[i]);
                    mtz.addTerm(-ins.bigM, x[i][j]);
                    mtzRows[i][j] = cplex.addGe(mtz, -ins.bigM, "MTZ_" + i + "_" + j + "_" + t);
                }
            }

            for (int i = 0; i < nodeCount; i++) {
                cplex.addLe(u[i], ins.Q, "LoadCap_" + i + "_" + t);
            }
            cplex.addEq(u[0], 0.0, "DepotLoad_" + t);
        }

        Result solve(double[] qBar, int[] zBar, double totalPickup) throws IloException {
            vehCap.setUB(-totalPickup);

            for (int j = 1; j <= n; j++) {
                double rhs = zBar[j];
                visitLink[j].setBounds(rhs, rhs);
            }
            updateMtzRhs(qBar);

            boolean solved = cplex.solve();
            String status = cplex.getStatus().toString();
            boolean optimal = status.startsWith("Optimal");
            if (!solved) {
                return new Result(false, false, status, Double.NaN);
            }
            return new Result(true, optimal, status, cplex.getObjValue());
        }

        void setSolveThreads(int threads) throws IloException {
            cplex.setParam(IloCplex.Param.Threads, Math.max(1, threads));
        }

        private void updateMtzRhs(double[] qBar) throws IloException {
            for (int i = 0; i < nodeCount; i++) {
                double rhs = -ins.bigM;
                if (isPickupNode(i, n)) {
                    rhs += qBar[i];
                }
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    mtzRows[i][j].setLB(rhs);
                }
            }
        }

        void close() {
            cplex.end();
        }
    }

    private static final class LpPeriodModel {
        private final Instance ins;
        private final int t;
        private final int n;
        private final int nodeCount;
        private final IloCplex cplex;
        private final IloRange vehCap;
        private final IloRange[] visitLink;
        private final IloRange[][] mtzRows;

        LpPeriodModel(Instance ins, int t) throws IloException {
            this.ins = ins;
            this.t = t;
            this.n = ins.n;
            this.nodeCount = ins.nodeCount;
            this.cplex = new IloCplex();
            configure(cplex);
            cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Dual);

            IloNumVar m = cplex.numVar(0.0, ins.K, "m_lp_" + t);
            IloNumVar[] u = new IloNumVar[nodeCount];
            IloNumVar[][] x = new IloNumVar[nodeCount][nodeCount];
            visitLink = new IloRange[n + 1];
            mtzRows = new IloRange[nodeCount][nodeCount];

            for (int i = 0; i < nodeCount; i++) {
                u[i] = cplex.numVar(0.0, Double.MAX_VALUE, "u_lp_" + i + "_" + t);
            }
            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    x[i][j] = cplex.numVar(0.0, 1.0, "x_lp_" + i + "_" + j + "_" + t);
                }
            }

            IloLinearNumExpr obj = cplex.linearNumExpr();
            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    obj.addTerm(ins.c[i][j], x[i][j]);
                }
            }
            cplex.addMinimize(obj);

            IloLinearNumExpr vehCapExpr = cplex.linearNumExpr();
            vehCapExpr.addTerm(-ins.Q, m);
            vehCap = cplex.addLe(vehCapExpr, 0.0, "VehicleCapLP_" + t);

            for (int i = 1; i <= n; i++) {
                IloLinearNumExpr out = cplex.linearNumExpr();
                IloLinearNumExpr in = cplex.linearNumExpr();
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    out.addTerm(1.0, x[i][j]);
                    in.addTerm(1.0, x[j][i]);
                }
                cplex.addEq(out, in, "FlowBalanceLP_" + i + "_" + t);
            }

            IloLinearNumExpr depart = cplex.linearNumExpr();
            for (int j = 1; j <= n; j++) {
                depart.addTerm(1.0, x[0][j]);
            }
            cplex.addEq(depart, m, "DepartFactoryLP_" + t);

            IloLinearNumExpr back = cplex.linearNumExpr();
            for (int i = 1; i <= n; i++) {
                back.addTerm(1.0, x[i][n + 1]);
            }
            cplex.addEq(back, m, "ReturnFactoryLP_" + t);

            for (int j = 1; j <= n; j++) {
                IloLinearNumExpr incoming = cplex.linearNumExpr();
                for (int i = 0; i < nodeCount; i++) {
                    if (i == j) {
                        continue;
                    }
                    incoming.addTerm(1.0, x[i][j]);
                }
                visitLink[j] = cplex.addEq(incoming, 0.0, "VisitLinkLP_" + j + "_" + t);
            }

            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    IloLinearNumExpr mtz = cplex.linearNumExpr();
                    mtz.addTerm(1.0, u[j]);
                    mtz.addTerm(-1.0, u[i]);
                    mtz.addTerm(-ins.bigM, x[i][j]);
                    mtzRows[i][j] = cplex.addGe(mtz, -ins.bigM, "MTZLP_" + i + "_" + j + "_" + t);
                }
            }

            for (int i = 0; i < nodeCount; i++) {
                cplex.addLe(u[i], ins.Q, "LoadCapLP_" + i + "_" + t);
            }
            cplex.addEq(u[0], 0.0, "DepotLoadLP_" + t);
        }

        DualCutResult solve(Instance ins, double[] qBar, int[] zBar, double totalPickup) throws IloException {
            vehCap.setUB(-totalPickup);

            for (int j = 1; j <= n; j++) {
                double rhs = zBar[j];
                visitLink[j].setBounds(rhs, rhs);
            }
            updateMtzRhs(qBar);

            boolean solved = cplex.solve();
            String status = cplex.getStatus().toString();
            boolean optimal = status.startsWith("Optimal");
            if (!solved) {
                return new DualCutResult(false, false, status, Double.NaN, Double.NaN, null, null);
            }

            double[] coeffZ = new double[n + 1];
            double[] coeffQ = new double[n + 1];

            double dualVehCap = cplex.getDual(vehCap);
            for (int i = 1; i <= n; i++) {
                coeffQ[i] += -dualVehCap;
            }

            for (int j = 1; j <= n; j++) {
                coeffZ[j] = cplex.getDual(visitLink[j]);
            }

            for (int i = 1; i <= n; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    coeffQ[i] += cplex.getDual(mtzRows[i][j]);
                }
            }

            double lpObj = cplex.getObjValue();
            double constant = lpObj;
            for (int i = 1; i <= n; i++) {
                constant -= coeffZ[i] * zBar[i];
                constant -= coeffQ[i] * qBar[i];
            }

            return new DualCutResult(true, optimal, status, lpObj, constant, coeffZ, coeffQ);
        }

        private void updateMtzRhs(double[] qBar) throws IloException {
            for (int i = 0; i < nodeCount; i++) {
                double rhs = -ins.bigM;
                if (isPickupNode(i, n)) {
                    rhs += qBar[i];
                }
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    mtzRows[i][j].setLB(rhs);
                }
            }
        }

        void close() {
            cplex.end();
        }
    }
}
