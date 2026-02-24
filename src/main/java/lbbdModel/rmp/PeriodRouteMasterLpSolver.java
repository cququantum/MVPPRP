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
 * Solve the LP relaxation of the set-partitioning CVRP subproblem via column generation,
 * and return duals used in lbbd.tex Eq. (3-18).
 */
public final class PeriodRouteMasterLpSolver {
    private static final boolean LOG_TO_CONSOLE = false;
    private static final double RC_EPS = 1e-8;
    private static final double ART_EPS = 1e-7;

    public static final class Result {
        public final boolean feasible;
        public final boolean optimal;
        public final boolean pricingProvedOptimal;
        public final boolean artificialClean;
        public final String status;
        public final double lpObjective;
        public final double[] dualUiByGlobal; // size n+1
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
            this.dualU0 = dualU0;
            this.generatedColumns = generatedColumns;
            this.customerCount = customerCount;
        }
    }

    public Result solve(Instance ins, int t, double[] qBar, int[] zBar) {
        if (t < 1 || t > ins.l) {
            throw new IllegalArgumentException("t must be in 1..l");
        }
        int[] activeCustomers = collectActiveCustomers(ins, qBar, zBar);
        int m = activeCustomers.length;
        if (m == 0) {
            return new Result(true, true, true, true, "TrivialZeroRmp", 0.0, new double[ins.n + 1], 0.0, 0, 0);
        }

        double[] qLocal = new double[m];
        for (int k = 0; k < m; k++) {
            qLocal[k] = qBar[activeCustomers[k]];
            if (qLocal[k] > ins.Q + 1e-9) {
                return new Result(false, false, false, false, "InfeasibleBySingleCustomerCapacity", Double.NaN,
                        null, Double.NaN, 0, m);
            }
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
                addRouteColumn(cplex, obj, coverEq, vehLimit, activeCustomers, singleton, routeVars, existingRouteKeys,
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
                    return new Result(false, false, false, false, status, Double.NaN, null, Double.NaN, generatedColumns, m);
                }

                lastDualUiLocal = new double[m];
                for (int k = 0; k < m; k++) {
                    lastDualUiLocal[k] = cplex.getDual(coverEq[k]);
                }
                lastDualU0 = cplex.getDual(vehLimit);

                PricingEspprcSolver.Result pricing = pricingSolver.findBestNegativeRoute(
                        ins,
                        activeCustomers,
                        qLocal,
                        lastDualUiLocal,
                        lastDualU0,
                        existingRouteKeys
                );

                if (!pricing.foundNegativeColumn) {
                    break;
                }

                addRouteColumn(
                        cplex, obj, coverEq, vehLimit, activeCustomers, pricing.route, routeVars, existingRouteKeys,
                        "cg_" + t + "_" + generatedColumns
                );
                generatedColumns++;
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

            return new Result(
                    true,
                    optimal,
                    true,
                    artificialClean,
                    artificialClean ? status : (status + "_ArtificialPositive"),
                    cplex.getObjValue(),
                    dualUiByGlobal,
                    lastDualU0,
                    generatedColumns,
                    m
            );
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve route-master LP at t=" + t, e);
        }
    }

    private static int[] collectActiveCustomers(Instance ins, double[] qBar, int[] zBar) {
        ArrayList<Integer> list = new ArrayList<Integer>();
        for (int i = 1; i <= ins.n; i++) {
            if (zBar[i] != 0 && qBar[i] > 1e-9) {
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

    private static void addRouteColumn(
            IloCplex cplex,
            IloObjective obj,
            IloRange[] coverEq,
            IloRange vehLimit,
            int[] activeCustomers,
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
        boolean[] coveredLocal = new boolean[activeCustomers.length];
        for (int pos = 0; pos < route.globalCustomersInOrder.length; pos++) {
            int g = route.globalCustomersInOrder[pos];
            int local = indexOf(activeCustomers, g);
            if (local >= 0 && !coveredLocal[local]) {
                col = col.and(cplex.column(coverEq[local], 1.0));
                coveredLocal[local] = true;
            }
        }
        IloNumVar var = cplex.numVar(col, 0.0, Double.MAX_VALUE, "xi_" + varName);
        routeVars.add(var);
    }

    private static int indexOf(int[] arr, int x) {
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] == x) {
                return i;
            }
        }
        return -1;
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
