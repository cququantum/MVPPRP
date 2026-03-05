package lbbdModel.rmp;

import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import instance.Instance;
import model.CplexConfig;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

/**
 * Solve the LP relaxation of the set-partitioning CVRP subproblem via column generation,
 * and return duals used in lbbd.tex Eq. (3-18).
 */
public final class PeriodRouteMasterLpSolver {
    private static final boolean LOG_TO_CONSOLE = false;
    private static final double RC_EPS = 1e-8;
    private static final double ART_EPS = 1e-7;
    private static final int EXACT_SUBSET_LP_MAX_ACTIVE_CUSTOMERS = 15;

    // Adaptive column generation parameters based on instance size
    private final int maxColumnsPerPricingSmall;
    private final int maxColumnsPerPricingMedium;
    private final int maxColumnsPerPricingLarge;

    private ExactSubsetRouteCostOracle subsetRouteOracle;
    private Instance subsetRouteOracleInstance;

    public static final class Result {
        public final boolean feasible;
        public final boolean optimal;
        public final boolean pricingProvedOptimal;
        public final boolean artificialClean;
        public final String status;
        public final double lpObjective;
        public final double[] dualUiByGlobal; // size n+1
        public final double[][] dualUiByPrev; // size [n+1][l+2], only period-t entries used
        public final double dualU0;
        public final int generatedColumns;
        public final int customerCount;

        public Result(
                boolean feasible,
                boolean optimal,
                boolean pricingProvedOptimal,
                boolean artificialClean,
                String status,
                double lpObjective,
                double[] dualUiByGlobal,
                double[][] dualUiByPrev,
                double dualU0,
                int generatedColumns,
                int customerCount
        ) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.pricingProvedOptimal = pricingProvedOptimal;
            this.artificialClean = artificialClean;
            this.status = status;
            this.lpObjective = lpObjective;
            this.dualUiByGlobal = dualUiByGlobal;
            this.dualUiByPrev = dualUiByPrev;
            this.dualU0 = dualU0;
            this.generatedColumns = generatedColumns;
            this.customerCount = customerCount;
        }
    }

    public PeriodRouteMasterLpSolver(Instance ins) {
        // Adaptive column generation parameters based on instance size
        if (ins.n <= 12) {
            // Small scale: preserve original parameters
            this.maxColumnsPerPricingSmall = 8;
            this.maxColumnsPerPricingMedium = 12;
            this.maxColumnsPerPricingLarge = 16;
        } else if (ins.n <= 20) {
            // Medium scale: moderate increase
            this.maxColumnsPerPricingSmall = 10;
            this.maxColumnsPerPricingMedium = 15;
            this.maxColumnsPerPricingLarge = 20;
        } else {
            // Large scale: aggressive increase
            this.maxColumnsPerPricingSmall = 12;
            this.maxColumnsPerPricingMedium = 18;
            this.maxColumnsPerPricingLarge = 24;
        }
    }

    public Result solve(Instance ins, int t, double[] qBar, int[] zBar, int[] prevVisitBySupplier) {
        if (t < 1 || t > ins.l) {
            throw new IllegalArgumentException("t must be in 1..l");
        }
        int[] activeCustomers = collectActiveCustomers(ins, qBar, zBar);
        int m = activeCustomers.length;
        if (m == 0) {
            return new Result(
                    true,
                    true,
                    true,
                    true,
                    "TrivialZeroRmp",
                    0.0,
                    new double[ins.n + 1],
                    new double[ins.n + 1][ins.l + 2],
                    0.0,
                    0,
                    0
            );
        }

        double[] qLocal = new double[m];
        int[] localIndexByGlobal = new int[ins.nodeCount];
        Arrays.fill(localIndexByGlobal, -1);
        for (int k = 0; k < m; k++) {
            qLocal[k] = qBar[activeCustomers[k]];
            localIndexByGlobal[activeCustomers[k]] = k;
            if (qLocal[k] > ins.Q + 1e-9) {
                return new Result(false, false, false, false, "InfeasibleBySingleCustomerCapacity", Double.NaN,
                        null, null, Double.NaN, 0, m);
            }
        }

        if (m <= EXACT_SUBSET_LP_MAX_ACTIVE_CUSTOMERS) {
            return solveBySubsetColumnsLp(ins, t, activeCustomers, qLocal, m, prevVisitBySupplier);
        }

        final double artificialPenalty = computeArtificialPenalty(ins, m);
        PricingEspprcSolver pricingSolver = new PricingEspprcSolver();
        try (IloCplex cplex = new IloCplex()) {
            configureLp(cplex);

            IloObjective obj = cplex.addMinimize();
            IloRange[] coverEq = new IloRange[m];
            for (int k = 0; k < m; k++) {
                coverEq[k] = cplex.addEq(cplex.linearNumExpr(), 1.0, "RMP_Cover_" + t + "_" + k);
            }
            IloRange vehLimit = cplex.addLe(cplex.linearNumExpr(), ins.K, "RMP_Vehicle_" + t);

            IloNumVar[] artificial = new IloNumVar[m];
            for (int k = 0; k < m; k++) {
                IloColumn col = cplex.column(obj, artificialPenalty).and(cplex.column(coverEq[k], 1.0));
                artificial[k] = cplex.numVar(col, 0.0, Double.MAX_VALUE, "a_" + t + "_" + k);
            }

            HashSet<String> existingRouteKeys = new HashSet<String>();
            ArrayList<IloNumVar> routeVars = new ArrayList<IloNumVar>();

            for (int k = 0; k < m; k++) {
                int g = activeCustomers[k];
                double cost = ins.c[0][g] + ins.c[g][ins.n + 1];
                RouteColumn singleton = new RouteColumn(new int[]{g}, cost, qLocal[k]);
                addRouteColumn(cplex, obj, coverEq, vehLimit, localIndexByGlobal, singleton, routeVars, existingRouteKeys,
                        "init_" + t + "_" + k);
            }

            int generatedColumns = 0;
            double[] lastDualUiLocal = null;
            double lastDualU0 = Double.NaN;
            String status = "Unknown";
            boolean optimal = false;

            while (true) {
                boolean solved = cplex.solve();
                status = cplex.getStatus().toString();
                optimal = solved && status.startsWith("Optimal");
                if (!solved) {
                    return new Result(false, false, false, false, status, Double.NaN,
                            null, null, Double.NaN, generatedColumns, m);
                }

                lastDualUiLocal = new double[m];
                for (int k = 0; k < m; k++) {
                    lastDualUiLocal[k] = cplex.getDual(coverEq[k]);
                }
                lastDualU0 = cplex.getDual(vehLimit);

                PricingEspprcSolver.Result pricing = pricingSolver.findNegativeRoutes(
                        ins,
                        activeCustomers,
                        qLocal,
                        lastDualUiLocal,
                        lastDualU0,
                        existingRouteKeys,
                        maxColumnsPerPricing(m)
                );

                if (!pricing.foundNegativeColumn) {
                    break;
                }

                ArrayList<RouteColumn> newRoutes = pricing.routes;
                if (newRoutes == null || newRoutes.isEmpty()) {
                    // Safety fallback: preserve old behavior if pricing only returns the best route.
                    addRouteColumn(
                            cplex, obj, coverEq, vehLimit, localIndexByGlobal, pricing.route, routeVars, existingRouteKeys,
                            "cg_" + t + "_" + generatedColumns
                    );
                    generatedColumns++;
                } else {
                    for (int r = 0; r < newRoutes.size(); r++) {
                        addRouteColumn(
                                cplex, obj, coverEq, vehLimit, localIndexByGlobal, newRoutes.get(r), routeVars, existingRouteKeys,
                                "cg_" + t + "_" + generatedColumns
                        );
                        generatedColumns++;
                    }
                }
            }

            boolean artificialClean = true;
            for (int k = 0; k < m; k++) {
                if (cplex.getValue(artificial[k]) > ART_EPS) {
                    artificialClean = false;
                    break;
                }
            }

            double[] dualUiByGlobal = new double[ins.n + 1];
            if (lastDualUiLocal != null) {
                for (int k = 0; k < m; k++) {
                    dualUiByGlobal[activeCustomers[k]] = lastDualUiLocal[k];
                }
            }
            double[][] dualUiByPrev = buildDualByPrev(ins, t, dualUiByGlobal, prevVisitBySupplier);

            return new Result(
                    true,
                    optimal,
                    true,
                    artificialClean,
                    artificialClean ? status : (status + "_ArtificialPositive"),
                    cplex.getObjValue(),
                    dualUiByGlobal,
                    dualUiByPrev,
                    lastDualU0,
                    generatedColumns,
                    m
            );
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve route-master LP at t=" + t, e);
        }
    }

    private Result solveBySubsetColumnsLp(
            Instance ins,
            int t,
            int[] activeCustomers,
            double[] qLocal,
            int m,
            int[] prevVisitBySupplier
    ) {
        if (subsetRouteOracle == null || subsetRouteOracleInstance != ins) {
            subsetRouteOracle = new ExactSubsetRouteCostOracle(ins);
            subsetRouteOracleInstance = ins;
        }

        double totalLoad = 0.0;
        for (int k = 0; k < m; k++) {
            totalLoad += qLocal[k];
        }
        if (totalLoad > ins.K * ins.Q + 1e-9) {
            return new Result(false, false, false, false, "InfeasibleByTotalCapacity", Double.NaN,
                    null, null, Double.NaN, 0, m);
        }

        try (IloCplex cplex = new IloCplex()) {
            configureLp(cplex);

            IloObjective obj = cplex.addMinimize();
            IloRange[] coverEq = new IloRange[m];
            for (int k = 0; k < m; k++) {
                coverEq[k] = cplex.addEq(cplex.linearNumExpr(), 1.0, "RMP_Subset_Cover_" + t + "_" + k);
            }
            IloRange vehLimit = cplex.addLe(cplex.linearNumExpr(), ins.K, "RMP_Subset_Vehicle_" + t);

            int localMaskCount = 1 << m;
            double[] subsetLoadLocal = new double[localMaskCount];
            int[] subsetGlobalMask = new int[localMaskCount];

            for (int mask = 1; mask < localMaskCount; mask++) {
                int bit = mask & -mask;
                int local = Integer.numberOfTrailingZeros(bit);
                subsetLoadLocal[mask] = subsetLoadLocal[mask ^ bit] + qLocal[local];
                subsetGlobalMask[mask] = subsetGlobalMask[mask ^ bit] | (1 << (activeCustomers[local] - 1));
            }

            int generatedColumns = 0;
            for (int mask = 1; mask < localMaskCount; mask++) {
                if (subsetLoadLocal[mask] > ins.Q + 1e-9) {
                    continue;
                }
                double routeCost = subsetRouteOracle.routeCostByGlobalMask[subsetGlobalMask[mask]];
                if (!(routeCost < Double.POSITIVE_INFINITY / 2.0)) {
                    continue;
                }

                IloColumn col = cplex.column(obj, routeCost).and(cplex.column(vehLimit, 1.0));
                int bits = mask;
                while (bits != 0) {
                    int bit = bits & -bits;
                    int local = Integer.numberOfTrailingZeros(bit);
                    col = col.and(cplex.column(coverEq[local], 1.0));
                    bits ^= bit;
                }
                cplex.numVar(col, 0.0, Double.MAX_VALUE, "xi_subset_" + t + "_" + mask);
                generatedColumns++;
            }

            boolean solved = cplex.solve();
            String status = cplex.getStatus().toString();
            boolean optimal = solved && status.startsWith("Optimal");
            if (!solved) {
                return new Result(false, false, false, false, status, Double.NaN,
                        null, null, Double.NaN, generatedColumns, m);
            }

            double[] dualUiByGlobal = new double[ins.n + 1];
            double[] rawDualUiLocal = new double[m];
            for (int k = 0; k < m; k++) {
                rawDualUiLocal[k] = cplex.getDual(coverEq[k]);
            }
            double rawDualU0 = cplex.getDual(vehLimit);

            DualNormalization dualNorm = normalizeSubsetLpDual(
                    rawDualUiLocal,
                    rawDualU0,
                    subsetLoadLocal,
                    subsetGlobalMask,
                    subsetRouteOracle.routeCostByGlobalMask,
                    cplex.getObjValue(),
                    ins.K,
                    ins,
                    m
            );
            if (!dualNorm.valid) {
                return new Result(false, false, false, false,
                        status + "_InvalidDualNormalization", Double.NaN,
                        null, null, Double.NaN, generatedColumns, m);
            }
            for (int k = 0; k < m; k++) {
                dualUiByGlobal[activeCustomers[k]] = dualNorm.dualUiLocal[k];
            }
            double dualU0 = dualNorm.dualU0;
            double[][] dualUiByPrev = buildDualByPrev(ins, t, dualUiByGlobal, prevVisitBySupplier);

            return new Result(
                    true,
                    optimal,
                    true,
                    true,
                    status + "_SubsetColumnsLP",
                    cplex.getObjValue(),
                    dualUiByGlobal,
                    dualUiByPrev,
                    dualU0,
                    generatedColumns,
                    m
            );
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve subset-columns route-master LP at t=" + t, e);
        }
    }

    private static double[][] buildDualByPrev(
            Instance ins,
            int t,
            double[] dualUiByGlobal,
            int[] prevVisitBySupplier
    ) {
        double[][] byPrev = new double[ins.n + 1][ins.l + 2];
        if (prevVisitBySupplier == null) {
            return byPrev;
        }
        for (int i = 1; i <= ins.n; i++) {
            if (dualUiByGlobal == null || i >= dualUiByGlobal.length) {
                continue;
            }
            int v = (i < prevVisitBySupplier.length) ? prevVisitBySupplier[i] : -1;
            if (v >= 0 && v <= t - 1 && v < byPrev[i].length) {
                byPrev[i][v] = dualUiByGlobal[i];
            }
        }
        return byPrev;
    }

    private static DualNormalization normalizeSubsetLpDual(
            double[] rawDualUiLocal,
            double rawDualU0,
            double[] subsetLoadLocal,
            int[] subsetGlobalMask,
            double[] routeCostByGlobalMask,
            double lpObjective,
            int vehicleLimitK,
            Instance ins,
            int m
    ) {
        double[] candUi = new double[m];
        if (isDualFeasibleCandidate(
                rawDualUiLocal, rawDualU0, subsetLoadLocal, subsetGlobalMask, routeCostByGlobalMask,
                lpObjective, vehicleLimitK, ins, m
        )) {
            System.arraycopy(rawDualUiLocal, 0, candUi, 0, m);
            return new DualNormalization(true, candUi, rawDualU0);
        }
        for (int k = 0; k < m; k++) {
            candUi[k] = -rawDualUiLocal[k];
        }
        double flippedU0 = -rawDualU0;
        if (isDualFeasibleCandidate(
                candUi, flippedU0, subsetLoadLocal, subsetGlobalMask, routeCostByGlobalMask,
                lpObjective, vehicleLimitK, ins, m
        )) {
            return new DualNormalization(true, candUi, flippedU0);
        }
        return new DualNormalization(false, null, Double.NaN);
    }

    private static boolean isDualFeasibleCandidate(
            double[] dualUiLocal,
            double dualU0,
            double[] subsetLoadLocal,
            int[] subsetGlobalMask,
            double[] routeCostByGlobalMask,
            double lpObjective,
            int vehicleLimitK,
            Instance ins,
            int m
    ) {
        if (dualU0 > 1e-7) {
            return false;
        }
        int localMaskCount = 1 << m;
        for (int mask = 1; mask < localMaskCount; mask++) {
            if (subsetLoadLocal[mask] > ins.Q + 1e-9) {
                continue;
            }
            double lhs = dualU0;
            int bits = mask;
            while (bits != 0) {
                int bit = bits & -bits;
                int local = Integer.numberOfTrailingZeros(bit);
                lhs += dualUiLocal[local];
                bits ^= bit;
            }
            double cost = routeCostByGlobalMask[subsetGlobalMask[mask]];
            if (lhs > cost + 1e-6) {
                return false;
            }
        }
        double dualObj = vehicleLimitK * dualU0;
        for (int k = 0; k < m; k++) {
            dualObj += dualUiLocal[k];
        }
        double tol = 1e-5 * Math.max(1.0, Math.abs(lpObjective));
        if (Math.abs(dualObj - lpObjective) > tol) {
            return false;
        }
        return true;
    }

    private static final class DualNormalization {
        final boolean valid;
        final double[] dualUiLocal;
        final double dualU0;

        DualNormalization(boolean valid, double[] dualUiLocal, double dualU0) {
            this.valid = valid;
            this.dualUiLocal = dualUiLocal;
            this.dualU0 = dualU0;
        }
    }

    private static int[] collectActiveCustomers(Instance ins, double[] qBar, int[] zBar) {
        ArrayList<Integer> list = new ArrayList<Integer>();
        for (int i = 1; i <= ins.n; i++) {
            if (zBar[i] != 0) {
                list.add(i);
            }
        }
        int[] arr = new int[list.size()];
        for (int k = 0; k < list.size(); k++) {
            arr[k] = list.get(k);
        }
        return arr;
    }

    private static double computeArtificialPenalty(Instance ins, int customerCount) {
        double maxArc = 0.0;
        for (int i = 0; i < ins.nodeCount; i++) {
            for (int j = 0; j < ins.nodeCount; j++) {
                if (i != j && ins.c[i][j] > maxArc) {
                    maxArc = ins.c[i][j];
                }
            }
        }
        return 1_000_000.0 + maxArc * (customerCount + 2);
    }

    private int maxColumnsPerPricing(int customerCount) {
        if (customerCount >= 13) {
            return maxColumnsPerPricingLarge;
        }
        if (customerCount >= 9) {
            return maxColumnsPerPricingMedium;
        }
        return maxColumnsPerPricingSmall;
    }

    private static void addRouteColumn(
            IloCplex cplex,
            IloObjective obj,
            IloRange[] coverEq,
            IloRange vehLimit,
            int[] localIndexByGlobal,
            RouteColumn route,
            ArrayList<IloNumVar> routeVars,
            HashSet<String> existingRouteKeys,
            String varName
    ) throws IloException {
        String key = route.key();
        if (existingRouteKeys.contains(key)) {
            return;
        }
        existingRouteKeys.add(key);

        IloColumn col = cplex.column(obj, route.cost).and(cplex.column(vehLimit, 1.0));
        boolean[] coveredLocal = new boolean[coverEq.length];
        for (int pos = 0; pos < route.globalCustomersInOrder.length; pos++) {
            int g = route.globalCustomersInOrder[pos];
            int local = (g >= 0 && g < localIndexByGlobal.length) ? localIndexByGlobal[g] : -1;
            if (local >= 0 && !coveredLocal[local]) {
                col = col.and(cplex.column(coverEq[local], 1.0));
                coveredLocal[local] = true;
            }
        }
        IloNumVar var = cplex.numVar(col, 0.0, Double.MAX_VALUE, "xi_" + varName);
        routeVars.add(var);
    }

    private static void configureLp(IloCplex cplex) throws IloException {
        cplex.setParam(IloCplex.Param.TimeLimit, CplexConfig.TIME_LIMIT_SEC);
        cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Dual);
        cplex.setParam(IloCplex.Param.Threads, 1);
        if (!LOG_TO_CONSOLE) {
            cplex.setOut(null);
            cplex.setWarning(null);
        }
    }

    private static final class ExactSubsetRouteCostOracle {
        private static final double INF = 1e100;

        final double[] routeCostByGlobalMask;

        ExactSubsetRouteCostOracle(Instance ins) {
            if (ins.n <= 0 || ins.n > 30) {
                throw new IllegalArgumentException("ExactSubsetRouteCostOracle requires 1 <= n <= 30");
            }
            int n = ins.n;
            int maskCount = 1 << n;
            this.routeCostByGlobalMask = new double[maskCount];
            double[][] pathToLast = new double[maskCount][n];
            for (int mask = 0; mask < maskCount; mask++) {
                Arrays.fill(pathToLast[mask], INF);
            }
            routeCostByGlobalMask[0] = 0.0;

            for (int j = 0; j < n; j++) {
                int mask = 1 << j;
                int gj = j + 1;
                pathToLast[mask][j] = ins.c[0][gj];
                routeCostByGlobalMask[mask] = ins.c[0][gj] + ins.c[gj][ins.n + 1];
            }

            for (int mask = 1; mask < maskCount; mask++) {
                if ((mask & (mask - 1)) == 0) {
                    continue;
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
                routeCostByGlobalMask[mask] = bestRoute;
            }
        }
    }
}
