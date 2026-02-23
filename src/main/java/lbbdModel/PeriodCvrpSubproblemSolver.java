package lbbdModel;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import instance.Instance;
import model.CplexConfig;

public final class PeriodCvrpSubproblemSolver {
    private static final boolean LOG_TO_CONSOLE = false;

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
        if (t < 1 || t > ins.l) {
            throw new IllegalArgumentException("t must be in 1..l");
        }
        if (qBar == null || qBar.length < ins.n + 1) {
            throw new IllegalArgumentException("qBar length must be at least n+1");
        }
        if (zBar == null || zBar.length < ins.n + 1) {
            throw new IllegalArgumentException("zBar length must be at least n+1");
        }

        double totalPickup = 0.0;
        int visitCount = 0;
        for (int i = 1; i <= ins.n; i++) {
            totalPickup += qBar[i];
            if (zBar[i] != 0) {
                visitCount++;
            }
        }
        if (visitCount == 0 && totalPickup <= 1e-9) {
            return new Result(true, true, "TrivialZero", 0.0);
        }
        if (totalPickup > ins.Q * ins.K + 1e-9) {
            return new Result(false, false, "InfeasibleByTotalCapacity", Double.NaN);
        }

        try (IloCplex cplex = new IloCplex()) {
            configure(cplex);

            int n = ins.n;
            int nodeCount = ins.nodeCount;

            IloNumVar m = cplex.intVar(0, ins.K, "m_" + t);
            IloNumVar[] u = new IloNumVar[nodeCount];
            IloNumVar[][] x = new IloNumVar[nodeCount][nodeCount];

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

            IloLinearNumExpr vehCap = cplex.linearNumExpr();
            vehCap.addTerm(-ins.Q, m);
            cplex.addLe(vehCap, -totalPickup, "VehicleCap_" + t);

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
                cplex.addEq(incoming, zBar[j], "VisitLink_" + j + "_" + t);
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
                    double rhs = -ins.bigM;
                    if (isPickupNode(i, n)) {
                        rhs += qBar[i];
                    }
                    cplex.addGe(mtz, rhs, "MTZ_" + i + "_" + j + "_" + t);
                }
            }

            for (int i = 0; i < nodeCount; i++) {
                cplex.addLe(u[i], ins.Q, "LoadCap_" + i + "_" + t);
            }
            cplex.addEq(u[0], 0.0, "DepotLoad_" + t);

            boolean solved = cplex.solve();
            String status = cplex.getStatus().toString();
            boolean optimal = status.startsWith("Optimal");
            if (!solved) {
                return new Result(false, false, status, Double.NaN);
            }
            return new Result(true, optimal, status, cplex.getObjValue());
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve CVRP subproblem at t=" + t, e);
        }
    }

    public DualCutResult solveLpDualCut(Instance ins, int t, double[] qBar, int[] zBar) {
        if (t < 1 || t > ins.l) {
            throw new IllegalArgumentException("t must be in 1..l");
        }
        if (qBar == null || qBar.length < ins.n + 1) {
            throw new IllegalArgumentException("qBar length must be at least n+1");
        }
        if (zBar == null || zBar.length < ins.n + 1) {
            throw new IllegalArgumentException("zBar length must be at least n+1");
        }

        double totalPickup = 0.0;
        int visitCount = 0;
        for (int i = 1; i <= ins.n; i++) {
            totalPickup += qBar[i];
            if (zBar[i] != 0) {
                visitCount++;
            }
        }
        if (visitCount == 0 && totalPickup <= 1e-9) {
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
        if (totalPickup > ins.Q * ins.K + 1e-9) {
            return new DualCutResult(false, false, "InfeasibleByTotalCapacity", Double.NaN, Double.NaN, null, null);
        }

        try (IloCplex cplex = new IloCplex()) {
            configure(cplex);

            int n = ins.n;
            int nodeCount = ins.nodeCount;

            IloNumVar m = cplex.numVar(0.0, Double.MAX_VALUE, "m_lp_" + t);
            IloNumVar[] u = new IloNumVar[nodeCount];
            IloNumVar[][] x = new IloNumVar[nodeCount][nodeCount];

            for (int i = 0; i < nodeCount; i++) {
                u[i] = cplex.numVar(0.0, Double.MAX_VALUE, "u_lp_" + i + "_" + t);
            }
            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    x[i][j] = cplex.numVar(0.0, Double.MAX_VALUE, "x_lp_" + i + "_" + j + "_" + t);
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

            IloRange vehCap;
            IloRange[] visitLink = new IloRange[n + 1];
            IloRange[][] mtzRows = new IloRange[nodeCount][nodeCount];

            IloLinearNumExpr vehCapExpr = cplex.linearNumExpr();
            vehCapExpr.addTerm(-ins.Q, m);
            vehCap = cplex.addLe(vehCapExpr, -totalPickup, "VehicleCapLP_" + t);

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
                visitLink[j] = cplex.addEq(incoming, zBar[j], "VisitLinkLP_" + j + "_" + t);
            }
            cplex.addLe(m, ins.K, "VehicleCountLP_" + t);

            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    IloLinearNumExpr mtz = cplex.linearNumExpr();
                    mtz.addTerm(1.0, u[j]);
                    mtz.addTerm(-1.0, u[i]);
                    mtz.addTerm(-ins.bigM, x[i][j]);
                    double rhs = -ins.bigM;
                    if (isPickupNode(i, n)) {
                        rhs += qBar[i];
                    }
                    mtzRows[i][j] = cplex.addGe(mtz, rhs, "MTZLP_" + i + "_" + j + "_" + t);
                }
            }

            for (int i = 0; i < nodeCount; i++) {
                cplex.addLe(u[i], ins.Q, "LoadCapLP_" + i + "_" + t);
            }
            cplex.addEq(u[0], 0.0, "DepotLoadLP_" + t);

            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (i == j) {
                        continue;
                    }
                    cplex.addLe(x[i][j], 1.0, "XUpperLP_" + i + "_" + j + "_" + t);
                }
            }

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
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve LP CVRP dual-cut subproblem at t=" + t, e);
        }
    }

    private static boolean isPickupNode(int node, int n) {
        return node >= 1 && node <= n;
    }

    private static void configure(IloCplex cplex) throws IloException {
        cplex.setParam(IloCplex.Param.TimeLimit, CplexConfig.TIME_LIMIT_SEC);
        cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, CplexConfig.MIP_GAP);
        if (!LOG_TO_CONSOLE) {
            cplex.setOut(null);
            cplex.setWarning(null);
        }
    }
}
