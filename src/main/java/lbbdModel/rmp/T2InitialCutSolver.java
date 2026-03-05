package lbbdModel.rmp;

import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import instance.Instance;
import model.CplexConfig;

import java.util.ArrayList;
import java.util.HashSet;

/**
 * Initial dual-cut solver T2 from speed.tex Section B.3.
 * It solves a global LP relaxation (production + inventory + lambda flow),
 * while routing in each period is represented by set-partitioning columns.
 */
public final class T2InitialCutSolver {
    private static final boolean LOG_TO_CONSOLE = false;
    private static final double ART_EPS = 1e-7;

    private final Instance ins;
    private final int maxColumnsPerPricingSmall;
    private final int maxColumnsPerPricingMedium;
    private final int maxColumnsPerPricingLarge;

    public static final class T2Result {
        public final boolean feasible;
        public final boolean optimal;
        public final boolean pricingProvedOptimal;
        public final boolean artificialClean;
        public final double lpObjective;
        public final double[][][] dualWByPeriod; // [l+1][n+1][l+2]
        public final double[] dualU0ByPeriod;    // [l+1]
        public final int totalGeneratedColumns;

        public T2Result(
                boolean feasible,
                boolean optimal,
                boolean pricingProvedOptimal,
                boolean artificialClean,
                double lpObjective,
                double[][][] dualWByPeriod,
                double[] dualU0ByPeriod,
                int totalGeneratedColumns
        ) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.pricingProvedOptimal = pricingProvedOptimal;
            this.artificialClean = artificialClean;
            this.lpObjective = lpObjective;
            this.dualWByPeriod = dualWByPeriod;
            this.dualU0ByPeriod = dualU0ByPeriod;
            this.totalGeneratedColumns = totalGeneratedColumns;
        }
    }

    public T2InitialCutSolver(Instance ins) {
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

    public T2Result solve() {
        int n = ins.n;
        int l = ins.l;
        PricingEspprcSolver pricingSolver = new PricingEspprcSolver();

        @SuppressWarnings("unchecked")
        HashSet<String>[] existingRouteKeysByPeriod = new HashSet[l + 1];
        int[][] scenarioSupplierByPeriod = new int[l + 1][];
        int[][] scenarioPrevByPeriod = new int[l + 1][];
        double[][] scenarioDemandByPeriod = new double[l + 1][];
        IloRange[][] coverEqByPeriod = new IloRange[l + 1][];
        IloNumVar[][] artificialByPeriod = new IloNumVar[l + 1][];
        IloRange[] vehicleLimit = new IloRange[l + 1];

        try (IloCplex cplex = new IloCplex()) {
            configureLp(cplex);

            IloLinearNumExpr objExpr = cplex.linearNumExpr();

            IloNumVar[] y = new IloNumVar[l + 1];
            IloNumVar[] p = new IloNumVar[l + 1];
            IloNumVar[] p0 = new IloNumVar[l + 1];
            IloNumVar[] i0 = new IloNumVar[l + 1];
            IloNumVar[][][] lambda = new IloNumVar[n + 1][l + 2][l + 2];

            for (int t = 1; t <= l; t++) {
                y[t] = cplex.numVar(0.0, 1.0, "t2_y_" + t);
                p[t] = cplex.numVar(0.0, Double.MAX_VALUE, "t2_p_" + t);
                p0[t] = cplex.numVar(0.0, Double.MAX_VALUE, "t2_P0_" + t);
                i0[t] = cplex.numVar(0.0, Double.MAX_VALUE, "t2_I0_" + t);
                objExpr.addTerm(ins.u, p[t]);
                objExpr.addTerm(ins.f, y[t]);
                objExpr.addTerm(ins.h0, i0[t]);
                objExpr.addTerm(ins.hp, p0[t]);
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l + 1; t++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (ins.g(i, v, t) > ins.Q + 1e-9) {
                            continue;
                        }
                        IloNumVar var = cplex.numVar(0.0, 1.0, "t2_lambda_" + i + "_" + v + "_" + t);
                        lambda[i][v][t] = var;
                        objExpr.addTerm(ins.e(i, v, t), var);
                    }
                }
            }
            IloObjective obj = cplex.addMinimize(objExpr);

            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr finishedBalance = cplex.linearNumExpr();
                double rhsFinished = ins.dt[t];
                if (t == 1) {
                    rhsFinished -= ins.P00;
                } else {
                    finishedBalance.addTerm(1.0, p0[t - 1]);
                }
                finishedBalance.addTerm(1.0, p[t]);
                finishedBalance.addTerm(-1.0, p0[t]);
                cplex.addEq(finishedBalance, rhsFinished, "T2_FinishedBalance_" + t);

                IloLinearNumExpr prodCap = cplex.linearNumExpr();
                prodCap.addTerm(1.0, p[t]);
                prodCap.addTerm(-ins.C, y[t]);
                cplex.addLe(prodCap, 0.0, "T2_ProdCap_" + t);

                cplex.addLe(p0[t], ins.Lp, "T2_FinishedCap_" + t);
            }

            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr factoryBalance = cplex.linearNumExpr();
                double rhsFactory = 0.0;
                if (t == 1) {
                    rhsFactory = -ins.I00;
                } else {
                    factoryBalance.addTerm(1.0, i0[t - 1]);
                }
                for (int i = 1; i <= n; i++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (lambda[i][v][t] != null) {
                            factoryBalance.addTerm(ins.g(i, v, t), lambda[i][v][t]);
                        }
                    }
                }
                factoryBalance.addTerm(-ins.k, p[t]);
                factoryBalance.addTerm(-1.0, i0[t]);
                cplex.addEq(factoryBalance, rhsFactory, "T2_FactoryBalance_" + t);

                cplex.addLe(i0[t], ins.L0, "T2_FactoryCap_" + t);
            }

            for (int i = 1; i <= n; i++) {
                IloLinearNumExpr startFlow = cplex.linearNumExpr();
                for (int t = 1; t <= ins.mu[i][0]; t++) {
                    if (lambda[i][0][t] != null) {
                        startFlow.addTerm(1.0, lambda[i][0][t]);
                    }
                }
                cplex.addEq(startFlow, 1.0, "T2_StartFlow_" + i);
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l; t++) {
                    IloLinearNumExpr balance = cplex.linearNumExpr();
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (lambda[i][v][t] != null) {
                            balance.addTerm(1.0, lambda[i][v][t]);
                        }
                    }
                    for (int tau = t + 1; tau <= ins.mu[i][t]; tau++) {
                        if (lambda[i][t][tau] != null) {
                            balance.addTerm(-1.0, lambda[i][t][tau]);
                        }
                    }
                    cplex.addEq(balance, 0.0, "T2_LambdaFlow_" + i + "_" + t);
                }
            }

            for (int i = 1; i <= n; i++) {
                IloLinearNumExpr terminalFlow = cplex.linearNumExpr();
                for (int t = ins.pi[i][l + 1]; t <= l; t++) {
                    if (lambda[i][t][l + 1] != null) {
                        terminalFlow.addTerm(1.0, lambda[i][t][l + 1]);
                    }
                }
                cplex.addEq(terminalFlow, 1.0, "T2_TerminalFlow_" + i);
            }

            final double artificialPenalty = computeArtificialPenalty(ins, ins.n);
            for (int t = 1; t <= l; t++) {
                ArrayList<Integer> supplierList = new ArrayList<Integer>();
                ArrayList<Integer> prevList = new ArrayList<Integer>();
                ArrayList<Double> demandList = new ArrayList<Double>();
                for (int i = 1; i <= n; i++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (lambda[i][v][t] == null) {
                            continue;
                        }
                        supplierList.add(i);
                        prevList.add(v);
                        demandList.add(ins.g(i, v, t));
                    }
                }

                int scenarioCount = supplierList.size();
                scenarioSupplierByPeriod[t] = new int[scenarioCount];
                scenarioPrevByPeriod[t] = new int[scenarioCount];
                scenarioDemandByPeriod[t] = new double[scenarioCount];
                for (int s = 0; s < scenarioCount; s++) {
                    scenarioSupplierByPeriod[t][s] = supplierList.get(s);
                    scenarioPrevByPeriod[t][s] = prevList.get(s);
                    scenarioDemandByPeriod[t][s] = demandList.get(s);
                }

                coverEqByPeriod[t] = new IloRange[scenarioCount];
                for (int s = 0; s < scenarioCount; s++) {
                    int supplier = scenarioSupplierByPeriod[t][s];
                    int prev = scenarioPrevByPeriod[t][s];
                    IloLinearNumExpr lhs = cplex.linearNumExpr();
                    lhs.addTerm(-1.0, lambda[supplier][prev][t]);
                    coverEqByPeriod[t][s] = cplex.addEq(lhs, 0.0, "T2_Cover_" + t + "_s" + s);
                }

                vehicleLimit[t] = cplex.addLe(cplex.linearNumExpr(), ins.K, "T2_Vehicle_" + t);
                existingRouteKeysByPeriod[t] = new HashSet<String>();

                artificialByPeriod[t] = new IloNumVar[scenarioCount];
                for (int s = 0; s < scenarioCount; s++) {
                    IloColumn col = cplex.column(obj, artificialPenalty).and(cplex.column(coverEqByPeriod[t][s], 1.0));
                    artificialByPeriod[t][s] = cplex.numVar(col, 0.0, Double.MAX_VALUE, "t2_a_" + t + "_s" + s);
                }

                for (int s = 0; s < scenarioCount; s++) {
                    int supplier = scenarioSupplierByPeriod[t][s];
                    double cost = ins.c[0][supplier] + ins.c[supplier][ins.n + 1];
                    ScenarioRouteColumn singleton = new ScenarioRouteColumn(
                            new int[]{s},
                            cost,
                            scenarioDemandByPeriod[t][s]
                    );
                    addScenarioRouteColumn(
                            cplex,
                            obj,
                            coverEqByPeriod[t],
                            vehicleLimit[t],
                            singleton,
                            existingRouteKeysByPeriod[t],
                            "t2_init_" + t + "_s" + s
                    );
                }
            }

            int totalGeneratedColumns = 0;
            boolean optimal = false;
            while (true) {
                boolean solved = cplex.solve();
                String status = cplex.getStatus().toString();
                optimal = solved && status.startsWith("Optimal");
                if (!solved || !optimal) {
                    return new T2Result(false, false, false, false, Double.NaN, null, null, totalGeneratedColumns);
                }

                boolean foundAnyNewColumn = false;
                for (int t = 1; t <= l; t++) {
                    int scenarioCount = scenarioSupplierByPeriod[t].length;
                    if (scenarioCount == 0) {
                        continue;
                    }

                    double[] dualByScenario = new double[scenarioCount];
                    double dualU0;
                    try {
                        for (int s = 0; s < scenarioCount; s++) {
                            dualByScenario[s] = cplex.getDual(coverEqByPeriod[t][s]);
                        }
                        dualU0 = cplex.getDual(vehicleLimit[t]);
                    } catch (IloException dualErr) {
                        return new T2Result(false, false, false, false, Double.NaN, null, null, totalGeneratedColumns);
                    }

                    PricingEspprcSolver.ScenarioResult pricing = pricingSolver.findNegativeScenarioRoutes(
                            ins,
                            scenarioSupplierByPeriod[t],
                            scenarioDemandByPeriod[t],
                            dualByScenario,
                            dualU0,
                            existingRouteKeysByPeriod[t],
                            maxColumnsPerPricing(Math.min(scenarioCount, ins.n))
                    );
                    if (!pricing.foundNegativeColumn) {
                        continue;
                    }

                    ArrayList<ScenarioRouteColumn> newRoutes = pricing.routes;
                    if (newRoutes == null || newRoutes.isEmpty()) {
                        addScenarioRouteColumn(
                                cplex,
                                obj,
                                coverEqByPeriod[t],
                                vehicleLimit[t],
                                pricing.route,
                                existingRouteKeysByPeriod[t],
                                "t2_cg_" + t + "_" + totalGeneratedColumns
                        );
                        totalGeneratedColumns++;
                        foundAnyNewColumn = true;
                    } else {
                        for (int r = 0; r < newRoutes.size(); r++) {
                            addScenarioRouteColumn(
                                    cplex,
                                    obj,
                                    coverEqByPeriod[t],
                                    vehicleLimit[t],
                                    newRoutes.get(r),
                                    existingRouteKeysByPeriod[t],
                                    "t2_cg_" + t + "_" + totalGeneratedColumns
                            );
                            totalGeneratedColumns++;
                            foundAnyNewColumn = true;
                        }
                    }
                }

                if (!foundAnyNewColumn) {
                    break;
                }
            }

            boolean artificialClean = true;
            for (int t = 1; t <= l && artificialClean; t++) {
                for (int s = 0; s < artificialByPeriod[t].length; s++) {
                    if (cplex.getValue(artificialByPeriod[t][s]) > ART_EPS) {
                        artificialClean = false;
                        break;
                    }
                }
            }

            double[][][] dualWByPeriod = new double[l + 1][n + 1][l + 2];
            double[] dualU0ByPeriod = new double[l + 1];
            try {
                for (int t = 1; t <= l; t++) {
                    dualU0ByPeriod[t] = cplex.getDual(vehicleLimit[t]);
                    for (int s = 0; s < scenarioSupplierByPeriod[t].length; s++) {
                        int supplier = scenarioSupplierByPeriod[t][s];
                        int prev = scenarioPrevByPeriod[t][s];
                        dualWByPeriod[t][supplier][prev] = cplex.getDual(coverEqByPeriod[t][s]);
                    }
                }
            } catch (IloException dualErr) {
                return new T2Result(false, false, false, false, Double.NaN, null, null, totalGeneratedColumns);
            }

            return new T2Result(
                    true,
                    optimal,
                    true,
                    artificialClean,
                    cplex.getObjValue(),
                    dualWByPeriod,
                    dualU0ByPeriod,
                    totalGeneratedColumns
            );
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve T2 initial LP", e);
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
            IloRange[] coverEq,
            IloRange vehicleLimit,
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

        IloColumn col = cplex.column(obj, route.cost).and(cplex.column(vehicleLimit, 1.0));
        boolean[] covered = new boolean[coverEq.length];
        for (int p = 0; p < route.scenarioIndicesInOrder.length; p++) {
            int s = route.scenarioIndicesInOrder[p];
            if (s >= 0 && s < coverEq.length && !covered[s]) {
                col = col.and(cplex.column(coverEq[s], 1.0));
                covered[s] = true;
            }
        }
        cplex.numVar(col, 0.0, Double.MAX_VALUE, "xi_t2_" + varName);
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

    private static void configureLp(IloCplex cplex) throws IloException {
        cplex.setParam(IloCplex.Param.TimeLimit, CplexConfig.TIME_LIMIT_SEC);
        cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Dual);
        cplex.setParam(IloCplex.Param.Threads, 1);
        if (!LOG_TO_CONSOLE) {
            cplex.setOut(null);
            cplex.setWarning(null);
        }
    }
}
