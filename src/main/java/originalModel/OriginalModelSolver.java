package originalModel;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import instance.Instance;
import model.CplexConfig;
import model.SolveResult;

public final class OriginalModelSolver {

    public SolveResult solve(Instance ins) {
        long startNs = System.nanoTime();
        try (IloCplex cplex = new IloCplex()) {
            configure(cplex);

            int n = ins.n;
            int l = ins.l;
            int nodeCount = ins.nodeCount;

            IloNumVar[] y = new IloNumVar[l + 1];
            IloNumVar[] p = new IloNumVar[l + 1];
            IloNumVar[] p0 = new IloNumVar[l + 1];
            IloNumVar[] i0 = new IloNumVar[l + 1];
            IloNumVar[] m = new IloNumVar[l + 1];

            IloNumVar[][] inventory = new IloNumVar[n + 1][l + 1];
            IloNumVar[][] q = new IloNumVar[n + 1][l + 1];
            IloNumVar[][] z = new IloNumVar[n + 1][l + 1];

            IloNumVar[][][] x = new IloNumVar[nodeCount][nodeCount][l + 1];
            IloNumVar[][] u = new IloNumVar[nodeCount][l + 1];

            for (int t = 1; t <= l; t++) {
                y[t] = cplex.boolVar("y_" + t);
                p[t] = cplex.numVar(0.0, Double.MAX_VALUE, "p_" + t);
                p0[t] = cplex.numVar(0.0, Double.MAX_VALUE, "P0_" + t);
                i0[t] = cplex.numVar(0.0, Double.MAX_VALUE, "I0_" + t);
                m[t] = cplex.intVar(0, ins.K, "m_" + t);
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l; t++) {
                    inventory[i][t] = cplex.numVar(0.0, Double.MAX_VALUE, "I_" + i + "_" + t);
                    q[i][t] = cplex.numVar(0.0, Double.MAX_VALUE, "q_" + i + "_" + t);
                    z[i][t] = cplex.boolVar("z_" + i + "_" + t);
                }
            }

            for (int i = 0; i < nodeCount; i++) {
                for (int t = 1; t <= l; t++) {
                    u[i][t] = cplex.numVar(0.0, Double.MAX_VALUE, "u_" + i + "_" + t);
                }
            }

            for (int t = 1; t <= l; t++) {
                for (int i = 0; i < nodeCount; i++) {
                    for (int j = 0; j < nodeCount; j++) {
                        if (i == j) {
                            continue;
                        }
                        x[i][j][t] = cplex.boolVar("x_" + i + "_" + j + "_" + t);
                    }
                }
            }

            IloLinearNumExpr obj = cplex.linearNumExpr();
            for (int t = 1; t <= l; t++) {
                for (int i = 0; i < nodeCount; i++) {
                    for (int j = 0; j < nodeCount; j++) {
                        if (i == j) {
                            continue;
                        }
                        obj.addTerm(ins.c[i][j], x[i][j][t]);
                    }
                }
                obj.addTerm(ins.u, p[t]);
                obj.addTerm(ins.f, y[t]);
                obj.addTerm(ins.h0, i0[t]);
                obj.addTerm(ins.hp, p0[t]);
                for (int i = 1; i <= n; i++) {
                    obj.addTerm(ins.hi[i], inventory[i][t]);
                }
            }
            cplex.addMinimize(obj);

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
                cplex.addEq(finishedBalance, rhsFinished, "FinishedBalance_" + t);

                IloLinearNumExpr prodCap = cplex.linearNumExpr();
                prodCap.addTerm(1.0, p[t]);
                prodCap.addTerm(-ins.C, y[t]);
                cplex.addLe(prodCap, 0.0, "ProdCap_" + t);

                cplex.addLe(p0[t], ins.Lp, "FinishedCap_" + t);
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l; t++) {
                    IloLinearNumExpr supplierBalance = cplex.linearNumExpr();
                    double rhsSupplier = -ins.s[i][t];
                    if (t == 1) {
                        rhsSupplier -= ins.Ii0[i];
                    } else {
                        supplierBalance.addTerm(1.0, inventory[i][t - 1]);
                    }
                    supplierBalance.addTerm(-1.0, q[i][t]);
                    supplierBalance.addTerm(-1.0, inventory[i][t]);
                    cplex.addEq(supplierBalance, rhsSupplier, "SupplierBalance_" + i + "_" + t);

                    IloLinearNumExpr pickupAll = cplex.linearNumExpr();
                    double rhsPickupAll = ins.s[i][t] - ins.bigM;
                    pickupAll.addTerm(1.0, q[i][t]);
                    if (t == 1) {
                        rhsPickupAll += ins.Ii0[i];
                    } else {
                        pickupAll.addTerm(-1.0, inventory[i][t - 1]);
                    }
                    pickupAll.addTerm(-ins.bigM, z[i][t]);
                    cplex.addGe(pickupAll, rhsPickupAll, "PickupAllLB_" + i + "_" + t);

                    IloLinearNumExpr rhsNoVisitNoPickup = cplex.linearNumExpr();
                    rhsNoVisitNoPickup.addTerm(ins.bigM, z[i][t]);
                    cplex.addLe(q[i][t], rhsNoVisitNoPickup, "NoVisitNoPickup_" + i + "_" + t);

                    IloLinearNumExpr visitZeroInventory = cplex.linearNumExpr();
                    visitZeroInventory.addTerm(1.0, inventory[i][t]);
                    visitZeroInventory.addTerm(ins.bigM, z[i][t]);
                    cplex.addLe(visitZeroInventory, ins.bigM, "VisitInventoryZero_" + i + "_" + t);

                    IloLinearNumExpr noVisitCarry = cplex.linearNumExpr();
                    double rhsNoVisitCarry = ins.s[i][t];
                    noVisitCarry.addTerm(1.0, inventory[i][t]);
                    if (t == 1) {
                        rhsNoVisitCarry += ins.Ii0[i];
                    } else {
                        noVisitCarry.addTerm(-1.0, inventory[i][t - 1]);
                    }
                    noVisitCarry.addTerm(ins.bigM, z[i][t]);
                    cplex.addGe(noVisitCarry, rhsNoVisitCarry, "NoVisitCarry_" + i + "_" + t);

                    cplex.addLe(inventory[i][t], ins.Li[i], "SupplierCap_" + i + "_" + t);
                }
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
                    factoryBalance.addTerm(1.0, q[i][t]);
                }
                factoryBalance.addTerm(-ins.k, p[t]);
                factoryBalance.addTerm(-1.0, i0[t]);
                cplex.addEq(factoryBalance, rhsFactory, "FactoryBalance_" + t);

                cplex.addLe(i0[t], ins.L0, "FactoryCap_" + t);
            }

            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr vehCap = cplex.linearNumExpr();
                for (int i = 1; i <= n; i++) {
                    vehCap.addTerm(1.0, q[i][t]);
                }
                vehCap.addTerm(-ins.Q, m[t]);
                cplex.addLe(vehCap, 0.0, "VehicleCap_" + t);
            }

            for (int t = 1; t <= l; t++) {
                for (int i = 1; i <= n; i++) {
                    IloLinearNumExpr out = cplex.linearNumExpr();
                    IloLinearNumExpr in = cplex.linearNumExpr();
                    for (int j = 0; j < nodeCount; j++) {
                        if (i == j) {
                            continue;
                        }
                        out.addTerm(1.0, x[i][j][t]);
                        in.addTerm(1.0, x[j][i][t]);
                    }
                    cplex.addEq(out, in, "FlowBalance_" + i + "_" + t);
                }
            }

            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr depart = cplex.linearNumExpr();
                for (int j = 1; j <= n; j++) {
                    depart.addTerm(1.0, x[0][j][t]);
                }
                cplex.addEq(depart, m[t], "DepartFactory_" + t);

                IloLinearNumExpr back = cplex.linearNumExpr();
                for (int i = 1; i <= n; i++) {
                    back.addTerm(1.0, x[i][n + 1][t]);
                }
                cplex.addEq(back, m[t], "ReturnFactory_" + t);
            }

            for (int t = 1; t <= l; t++) {
                for (int j = 1; j <= n; j++) {
                    IloLinearNumExpr incoming = cplex.linearNumExpr();
                    for (int i = 0; i < nodeCount; i++) {
                        if (i == j) {
                            continue;
                        }
                        incoming.addTerm(1.0, x[i][j][t]);
                    }
                    cplex.addEq(incoming, z[j][t], "VisitLink_" + j + "_" + t);
                }
            }

            for (int t = 1; t <= l; t++) {
                for (int i = 0; i < nodeCount; i++) {
                    for (int j = 0; j < nodeCount; j++) {
                        if (i == j) {
                            continue;
                        }
                        IloLinearNumExpr mtz = cplex.linearNumExpr();
                        mtz.addTerm(1.0, u[j][t]);
                        mtz.addTerm(-1.0, u[i][t]);
                        mtz.addTerm(-ins.bigM, x[i][j][t]);
                        if (isPickupNode(i, n)) {
                            mtz.addTerm(-1.0, q[i][t]);
                        }
                        cplex.addGe(mtz, -ins.bigM, "MTZ_" + i + "_" + j + "_" + t);
                    }
                }
            }

            for (int t = 1; t <= l; t++) {
                for (int i = 0; i < nodeCount; i++) {
                    cplex.addLe(u[i][t], ins.Q, "LoadCap_" + i + "_" + t);
                }
                cplex.addEq(u[0][t], 0.0, "DepotLoad_" + t);
                cplex.addLe(m[t], ins.K, "VehicleCount_" + t);
            }

            boolean solved = cplex.solve();
            return buildResult("OriginalModel", solved, cplex, startNs);
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve OriginalModel", e);
        }
    }

    private static boolean isPickupNode(int node, int n) {
        return node >= 1 && node <= n;
    }

    private static void configure(IloCplex cplex) throws IloException {
        cplex.setParam(IloCplex.Param.TimeLimit, CplexConfig.TIME_LIMIT_SEC);
        cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, CplexConfig.MIP_GAP);
        if (!CplexConfig.LOG_TO_CONSOLE) {
            cplex.setOut(null);
            cplex.setWarning(null);
        }
    }

    private static SolveResult buildResult(String modelName, boolean solved, IloCplex cplex, long startNs) throws IloException {
        double sec = (System.nanoTime() - startNs) / 1_000_000_000.0;
        String status = cplex.getStatus().toString();
        boolean optimal = status.startsWith("Optimal");

        double objective = Double.NaN;
        if (solved) {
            objective = cplex.getObjValue();
        }

        double bestBound = safeBestBound(cplex);
        double gap = safeGap(cplex);

        return new SolveResult(modelName, solved, optimal, status, objective, bestBound, gap, sec);
    }

    private static double safeBestBound(IloCplex cplex) {
        try {
            return cplex.getBestObjValue();
        } catch (IloException e) {
            return Double.NaN;
        }
    }

    private static double safeGap(IloCplex cplex) {
        try {
            return cplex.getMIPRelativeGap();
        } catch (IloException e) {
            return Double.NaN;
        }
    }
}
