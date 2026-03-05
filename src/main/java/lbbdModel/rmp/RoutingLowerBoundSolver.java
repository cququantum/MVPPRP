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
import java.util.HashSet;

/**
 * Compute LB_R in lbbd.tex Section F by solving RF LP relaxation via column generation.
 */
public final class RoutingLowerBoundSolver {
    private static final boolean LOG_TO_CONSOLE = false;
    private static final double EPS = 1e-9;
    private static final double ART_EPS = 1e-7;
    private static final int MAX_COLUMNS_PER_PRICING = 16;

    public static final class Result {
        public final boolean feasible;
        public final boolean optimal;
        public final String status;
        public final double lbR;
        public final int iterations;
        public final int generatedColumns;

        public Result(
                boolean feasible,
                boolean optimal,
                String status,
                double lbR,
                int iterations,
                int generatedColumns
        ) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.status = status;
            this.lbR = lbR;
            this.iterations = iterations;
            this.generatedColumns = generatedColumns;
        }
    }

    private final Instance ins;

    public RoutingLowerBoundSolver(Instance ins) {
        this.ins = ins;
    }

    public Result solve() {
        int n = ins.n;
        double[] qMin = new double[n + 1];
        double[] serviceTimes = new double[n + 1];
        double totalLoad = 0.0;

        for (int i = 1; i <= n; i++) {
            double q = Double.POSITIVE_INFINITY;
            double horizonLoad = ins.Ii0[i];
            for (int t = 1; t <= ins.l; t++) {
                q = Math.min(q, ins.s[i][t]);
                horizonLoad += ins.s[i][t];
            }
            if (Double.isNaN(q) || Double.isInfinite(q)) {
                q = 0.0;
            }
            qMin[i] = Math.max(0.0, q);
            serviceTimes[i] = Math.max(0.0, Math.ceil(horizonLoad / ins.Q - EPS));
            totalLoad += horizonLoad;
            if (serviceTimes[i] > 0.0 && qMin[i] > ins.Q + EPS) {
                return new Result(false, false, "InfeasibleBySingleNodeMinDemand", Double.NaN, 0, 0);
            }
        }

        double mL = Math.max(0.0, Math.ceil(totalLoad / ins.Q - EPS));
        double mU = ins.K * ins.l;
        if (mL - mU > EPS) {
            return new Result(false, false, "InfeasibleByRouteCountBounds", Double.NaN, 0, 0);
        }

        boolean allZero = true;
        for (int i = 1; i <= n; i++) {
            if (serviceTimes[i] > EPS) {
                allZero = false;
                break;
            }
        }
        if (allZero && mL <= EPS) {
            return new Result(true, true, "TrivialZeroLBR", 0.0, 0, 0);
        }

        int[] customers = new int[n];
        double[] qLocal = new double[n];
        for (int k = 0; k < n; k++) {
            customers[k] = k + 1;
            qLocal[k] = qMin[k + 1];
        }

        PricingEspprcSolver pricingSolver = new PricingEspprcSolver();
        try (IloCplex cplex = new IloCplex()) {
            configureLp(cplex);
            IloObjective obj = cplex.addMinimize();
            double artPenalty = computeArtificialPenalty(ins, n);

            IloRange[] coverGe = new IloRange[n];
            for (int i = 1; i <= n; i++) {
                coverGe[i - 1] = cplex.addGe(cplex.linearNumExpr(), serviceTimes[i], "RF_Cover_" + i);
            }
            IloRange routeCountMin = cplex.addGe(cplex.linearNumExpr(), mL, "RF_RouteMin");
            IloRange routeCountMax = cplex.addLe(cplex.linearNumExpr(), mU, "RF_RouteMax");

            IloNumVar[] coverArtificial = new IloNumVar[n];
            for (int i = 0; i < n; i++) {
                IloColumn col = cplex.column(obj, artPenalty).and(cplex.column(coverGe[i], 1.0));
                coverArtificial[i] = cplex.numVar(col, 0.0, Double.MAX_VALUE, "rf_art_cover_" + (i + 1));
            }

            HashSet<String> existingRouteKeys = new HashSet<String>();

            for (int i = 1; i <= n; i++) {
                if (qMin[i] > ins.Q + EPS) {
                    continue;
                }
                RouteColumn singleton = new RouteColumn(
                        new int[]{i},
                        ins.c[0][i] + ins.c[i][ins.n + 1],
                        qMin[i]
                );
                addRouteColumn(cplex, obj, coverGe, routeCountMin, routeCountMax, singleton, existingRouteKeys, "rf_init_" + i);
            }

            int iteration = 0;
            int generatedColumns = 0;
            String status = "Unknown";
            boolean optimal = false;

            while (true) {
                iteration++;
                boolean solved = cplex.solve();
                status = cplex.getStatus().toString();
                optimal = solved && status.startsWith("Optimal");
                if (!solved) {
                    return new Result(false, false, status, Double.NaN, iteration, generatedColumns);
                }

                double[] dualCover = new double[n];
                for (int k = 0; k < n; k++) {
                    dualCover[k] = cplex.getDual(coverGe[k]);
                }
                double dualRoute = cplex.getDual(routeCountMin) + cplex.getDual(routeCountMax);

                PricingEspprcSolver.Result pricing = pricingSolver.findNegativeRoutes(
                        ins,
                        customers,
                        qLocal,
                        dualCover,
                        dualRoute,
                        existingRouteKeys,
                        MAX_COLUMNS_PER_PRICING
                );

                if (!pricing.foundNegativeColumn) {
                    break;
                }

                ArrayList<RouteColumn> newRoutes = pricing.routes;
                if (newRoutes == null || newRoutes.isEmpty()) {
                    addRouteColumn(cplex, obj, coverGe, routeCountMin, routeCountMax, pricing.route, existingRouteKeys,
                            "rf_cg_" + iteration + "_" + generatedColumns);
                    generatedColumns++;
                } else {
                    for (int r = 0; r < newRoutes.size(); r++) {
                        addRouteColumn(cplex, obj, coverGe, routeCountMin, routeCountMax, newRoutes.get(r), existingRouteKeys,
                                "rf_cg_" + iteration + "_" + generatedColumns);
                        generatedColumns++;
                    }
                }
            }

            boolean artificialClean = true;
            for (int i = 0; i < n; i++) {
                if (cplex.getValue(coverArtificial[i]) > ART_EPS) {
                    artificialClean = false;
                    break;
                }
            }
            if (!artificialClean) {
                return new Result(false, false, status + "_ArtificialPositive", Double.NaN, iteration, generatedColumns);
            }

            double lbR = cplex.getObjValue();
            return new Result(true, optimal, status, lbR, iteration, generatedColumns);
        } catch (IloException e) {
            throw new RuntimeException("Failed to compute LB_R", e);
        }
    }

    private static void addRouteColumn(
            IloCplex cplex,
            IloObjective obj,
            IloRange[] coverGe,
            IloRange routeCountMin,
            IloRange routeCountMax,
            RouteColumn route,
            HashSet<String> existingRouteKeys,
            String name
    ) throws IloException {
        if (route == null) {
            return;
        }
        String key = route.key();
        if (existingRouteKeys.contains(key)) {
            return;
        }
        existingRouteKeys.add(key);

        IloColumn col = cplex.column(obj, route.cost)
                .and(cplex.column(routeCountMin, 1.0))
                .and(cplex.column(routeCountMax, 1.0));
        for (int idx = 0; idx < route.globalCustomersInOrder.length; idx++) {
            int g = route.globalCustomersInOrder[idx];
            if (g >= 1 && g <= coverGe.length) {
                col = col.and(cplex.column(coverGe[g - 1], 1.0));
            }
        }
        cplex.numVar(col, 0.0, Double.MAX_VALUE, "zeta_" + name);
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

    private static double computeArtificialPenalty(Instance ins, int customerCount) {
        double maxArc = 0.0;
        for (int i = 0; i < ins.nodeCount; i++) {
            for (int j = 0; j < ins.nodeCount; j++) {
                if (i == j) {
                    continue;
                }
                if (ins.c[i][j] > maxArc) {
                    maxArc = ins.c[i][j];
                }
            }
        }
        return 1_000_000.0 + maxArc * (customerCount + 2);
    }
}
