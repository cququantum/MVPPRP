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
 * Initial dual-cut solver T1 from speed.tex Section B.2.
 * For each period t, solve a scenario set-covering LP by column generation:
 *   min sum_r b_r * xi_r
 *   s.t. sum_r a_{s,r} * xi_r >= 1, for every scenario s=(i,v)
 *        sum_r xi_r <= K
 *        xi_r >= 0
 */
public final class T1InitialCutSolver {
    private static final boolean LOG_TO_CONSOLE = false;
    private static final double ART_EPS = 1e-7;
    private static final long NANOS_PER_SECOND = 1_000_000_000L;

    private final Instance ins;
    private final int maxColumnsPerPricingSmall;
    private final int maxColumnsPerPricingMedium;
    private final int maxColumnsPerPricingLarge;

    public static final class T1Result {
        public final boolean feasible;
        public final boolean optimal;
        public final boolean pricingProvedOptimal;
        public final boolean artificialClean;
        public final double lpObjective;
        public final double[][] dualW; // [n+1][l+2], indexed by [i][v]
        public final double dualU0;
        public final int generatedColumns;

        public T1Result(
                boolean feasible,
                boolean optimal,
                boolean pricingProvedOptimal,
                boolean artificialClean,
                double lpObjective,
                double[][] dualW,
                double dualU0,
                int generatedColumns
        ) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.pricingProvedOptimal = pricingProvedOptimal;
            this.artificialClean = artificialClean;
            this.lpObjective = lpObjective;
            this.dualW = dualW;
            this.dualU0 = dualU0;
            this.generatedColumns = generatedColumns;
        }
    }

    public T1InitialCutSolver(Instance ins) {
        this.ins = ins;
        if (ins.n <= 12) {
            this.maxColumnsPerPricingSmall = 8;
            this.maxColumnsPerPricingMedium = 12;
            this.maxColumnsPerPricingLarge = 16;
        } else if (ins.n <= 20) {
            this.maxColumnsPerPricingSmall = 10;
            this.maxColumnsPerPricingMedium = 15;
            this.maxColumnsPerPricingLarge = 20;
        } else {
            this.maxColumnsPerPricingSmall = 12;
            this.maxColumnsPerPricingMedium = 18;
            this.maxColumnsPerPricingLarge = 24;
        }
    }

    public T1Result solve(int t) {
        return solve(t, CplexConfig.TIME_LIMIT_SEC);
    }

    public T1Result solve(int t, double timeLimitSec) {
        if (timeLimitSec <= 0.0) {
            return new T1Result(false, false, false, false, Double.NaN, null, Double.NaN, 0);
        }
        long deadlineNs = System.nanoTime() + Math.max(1L, Math.round(timeLimitSec * NANOS_PER_SECOND));
        if (t < 1 || t > ins.l) {
            throw new IllegalArgumentException("t must be in 1..l");
        }

        ArrayList<Integer> supplierList = new ArrayList<Integer>();
        ArrayList<Integer> prevList = new ArrayList<Integer>();
        ArrayList<Double> demandList = new ArrayList<Double>();
        for (int i = 1; i <= ins.n; i++) {
            for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                double demand = ins.g(i, v, t);
                if (demand > ins.Q + 1e-9) {
                    continue;
                }
                supplierList.add(i);
                prevList.add(v);
                demandList.add(demand);
            }
        }

        int scenarioCount = supplierList.size();
        if (scenarioCount == 0) {
            return new T1Result(
                    true,
                    true,
                    true,
                    true,
                    0.0,
                    new double[ins.n + 1][ins.l + 2],
                    0.0,
                    0
            );
        }

        int[] scenarioSupplier = new int[scenarioCount];
        int[] scenarioPrev = new int[scenarioCount];
        double[] scenarioDemand = new double[scenarioCount];
        for (int s = 0; s < scenarioCount; s++) {
            scenarioSupplier[s] = supplierList.get(s);
            scenarioPrev[s] = prevList.get(s);
            scenarioDemand[s] = demandList.get(s);
        }

        final double artificialPenalty = computeArtificialPenalty(ins, scenarioCount);
        PricingEspprcSolver pricingSolver = new PricingEspprcSolver();

        try (IloCplex cplex = new IloCplex()) {
            configureLp(cplex, timeLimitSec);

            IloObjective obj = cplex.addMinimize();
            IloRange[] coverGe = new IloRange[scenarioCount];
            for (int s = 0; s < scenarioCount; s++) {
                coverGe[s] = cplex.addGe(cplex.linearNumExpr(), 1.0, "T1_Cover_" + t + "_s" + s);
            }
            IloRange vehLimit = cplex.addLe(cplex.linearNumExpr(), ins.K, "T1_Vehicle_" + t);

            IloNumVar[] artificial = new IloNumVar[scenarioCount];
            for (int s = 0; s < scenarioCount; s++) {
                IloColumn col = cplex.column(obj, artificialPenalty).and(cplex.column(coverGe[s], 1.0));
                artificial[s] = cplex.numVar(col, 0.0, Double.MAX_VALUE, "t1_a_" + t + "_s" + s);
            }

            HashSet<String> existingRouteKeys = new HashSet<String>();
            for (int s = 0; s < scenarioCount; s++) {
                int supplier = scenarioSupplier[s];
                double cost = ins.c[0][supplier] + ins.c[supplier][ins.n + 1];
                ScenarioRouteColumn singleton = new ScenarioRouteColumn(new int[]{s}, cost, scenarioDemand[s]);
                addScenarioRouteColumn(
                        cplex,
                        obj,
                        coverGe,
                        vehLimit,
                        singleton,
                        existingRouteKeys,
                        "t1_init_" + t + "_s" + s
                );
            }

            int generatedColumns = 0;
            double[] lastDualByScenario = null;
            double lastDualU0 = Double.NaN;
            boolean optimal = false;

            while (true) {
                if (isPastDeadline(deadlineNs)) {
                    return new T1Result(false, false, false, false, Double.NaN, null, Double.NaN, generatedColumns);
                }
                cplex.setParam(IloCplex.Param.TimeLimit, normalizedTimeLimitSec(remainingSeconds(deadlineNs)));
                boolean solved = cplex.solve();
                String status = cplex.getStatus().toString();
                optimal = solved && status.startsWith("Optimal");
                if (!solved) {
                    return new T1Result(false, false, false, false, Double.NaN, null, Double.NaN, generatedColumns);
                }

                lastDualByScenario = new double[scenarioCount];
                for (int s = 0; s < scenarioCount; s++) {
                    lastDualByScenario[s] = cplex.getDual(coverGe[s]);
                }
                lastDualU0 = cplex.getDual(vehLimit);

                PricingEspprcSolver.ScenarioResult pricing = pricingSolver.findNegativeScenarioRoutes(
                        ins,
                        scenarioSupplier,
                        scenarioDemand,
                        lastDualByScenario,
                        lastDualU0,
                        existingRouteKeys,
                        maxColumnsPerPricing(Math.min(scenarioCount, ins.n))
                );

                if (!pricing.foundNegativeColumn) {
                    break;
                }

                ArrayList<ScenarioRouteColumn> newRoutes = pricing.routes;
                if (newRoutes == null || newRoutes.isEmpty()) {
                    addScenarioRouteColumn(
                            cplex,
                            obj,
                            coverGe,
                            vehLimit,
                            pricing.route,
                            existingRouteKeys,
                            "t1_cg_" + t + "_" + generatedColumns
                    );
                    generatedColumns++;
                } else {
                    for (int r = 0; r < newRoutes.size(); r++) {
                        addScenarioRouteColumn(
                                cplex,
                                obj,
                                coverGe,
                                vehLimit,
                                newRoutes.get(r),
                                existingRouteKeys,
                                "t1_cg_" + t + "_" + generatedColumns
                        );
                        generatedColumns++;
                    }
                }
            }

            boolean artificialClean = true;
            for (int s = 0; s < scenarioCount; s++) {
                if (cplex.getValue(artificial[s]) > ART_EPS) {
                    artificialClean = false;
                    break;
                }
            }

            double[][] dualW = new double[ins.n + 1][ins.l + 2];
            if (lastDualByScenario != null) {
                for (int s = 0; s < scenarioCount; s++) {
                    dualW[scenarioSupplier[s]][scenarioPrev[s]] = lastDualByScenario[s];
                }
            }

            return new T1Result(
                    true,
                    optimal,
                    true,
                    artificialClean,
                    cplex.getObjValue(),
                    dualW,
                    lastDualU0,
                    generatedColumns
            );
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve T1 initial LP at t=" + t, e);
        }
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

    private static void addScenarioRouteColumn(
            IloCplex cplex,
            IloObjective obj,
            IloRange[] coverGe,
            IloRange vehLimit,
            ScenarioRouteColumn route,
            HashSet<String> existingRouteKeys,
            String varName
    ) throws IloException {
        if (route == null) {
            return;
        }
        String key = route.key();
        if (existingRouteKeys.contains(key)) {
            return;
        }
        existingRouteKeys.add(key);

        IloColumn col = cplex.column(obj, route.cost).and(cplex.column(vehLimit, 1.0));
        boolean[] covered = new boolean[coverGe.length];
        for (int p = 0; p < route.scenarioIndicesInOrder.length; p++) {
            int s = route.scenarioIndicesInOrder[p];
            if (s >= 0 && s < coverGe.length && !covered[s]) {
                col = col.and(cplex.column(coverGe[s], 1.0));
                covered[s] = true;
            }
        }
        cplex.numVar(col, 0.0, Double.MAX_VALUE, "xi_t1_" + varName);
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

    private static boolean isPastDeadline(long deadlineNs) {
        return System.nanoTime() >= deadlineNs;
    }

    private static double remainingSeconds(long deadlineNs) {
        long remainingNs = deadlineNs - System.nanoTime();
        if (remainingNs <= 0L) {
            return 0.0;
        }
        return remainingNs / (double) NANOS_PER_SECOND;
    }

    private static double normalizedTimeLimitSec(double timeLimitSec) {
        return Math.max(1e-3, timeLimitSec);
    }

    private static void configureLp(IloCplex cplex, double timeLimitSec) throws IloException {
        cplex.setParam(IloCplex.Param.TimeLimit, normalizedTimeLimitSec(timeLimitSec));
        cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Dual);
        cplex.setParam(IloCplex.Param.Threads, 1);
        if (!LOG_TO_CONSOLE) {
            cplex.setOut(null);
            cplex.setWarning(null);
        }
    }
}
