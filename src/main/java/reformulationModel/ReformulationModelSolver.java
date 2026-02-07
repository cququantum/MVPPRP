package reformulationModel;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import instance.Instance;
import model.CplexConfig;
import model.SolveResult;

public final class ReformulationModelSolver {

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

            IloNumVar[][] z = new IloNumVar[n + 1][l + 1];
            IloNumVar[][] u = new IloNumVar[nodeCount][l + 1];
            IloNumVar[][][] x = new IloNumVar[nodeCount][nodeCount][l + 1];
            IloNumVar[][][] lambda = new IloNumVar[n + 1][l + 2][l + 2];

            for (int t = 1; t <= l; t++) {
                y[t] = cplex.boolVar("y_" + t);
                p[t] = cplex.numVar(0.0, Double.MAX_VALUE, "p_" + t);
                p0[t] = cplex.numVar(0.0, Double.MAX_VALUE, "P0_" + t);
                i0[t] = cplex.numVar(0.0, Double.MAX_VALUE, "I0_" + t);
                m[t] = cplex.intVar(0, ins.K, "m_" + t);
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l; t++) {
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

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l + 1; t++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        lambda[i][v][t] = cplex.boolVar("lambda_" + i + "_" + v + "_" + t);
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
            }
            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l + 1; t++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        obj.addTerm(holdingCostOnArc(ins, i, v, t), lambda[i][v][t]);
                    }
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
                        factoryBalance.addTerm(ins.g(i, v, t), lambda[i][v][t]);
                    }
                }
                factoryBalance.addTerm(-ins.k, p[t]);
                factoryBalance.addTerm(-1.0, i0[t]);
                cplex.addEq(factoryBalance, rhsFactory, "FactoryBalance_" + t);

                cplex.addLe(i0[t], ins.L0, "FactoryCap_" + t);
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l; t++) {
                    IloLinearNumExpr link = cplex.linearNumExpr();
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        link.addTerm(1.0, lambda[i][v][t]);
                    }
                    link.addTerm(-1.0, z[i][t]);
                    cplex.addEq(link, 0.0, "LinkZLambda_" + i + "_" + t);
                }
            }

            for (int i = 1; i <= n; i++) {
                IloLinearNumExpr startFlow = cplex.linearNumExpr();
                for (int t = 1; t <= ins.mu[i][0]; t++) {
                    if (lambda[i][0][t] != null) {
                        startFlow.addTerm(1.0, lambda[i][0][t]);
                    }
                }
                cplex.addEq(startFlow, 1.0, "StartFlow_" + i);
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l; t++) {
                    IloLinearNumExpr balance = cplex.linearNumExpr();
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        balance.addTerm(1.0, lambda[i][v][t]);
                    }
                    for (int tau = t + 1; tau <= ins.mu[i][t]; tau++) {
                        if (lambda[i][t][tau] != null) {
                            balance.addTerm(-1.0, lambda[i][t][tau]);
                        }
                    }
                    cplex.addEq(balance, 0.0, "LambdaFlow_" + i + "_" + t);
                }
            }

            for (int i = 1; i <= n; i++) {
                IloLinearNumExpr terminalFlow = cplex.linearNumExpr();
                for (int t = ins.pi[i][l + 1]; t <= l; t++) {
                    if (lambda[i][t][l + 1] != null) {
                        terminalFlow.addTerm(1.0, lambda[i][t][l + 1]);
                    }
                }
                cplex.addEq(terminalFlow, 1.0, "TerminalFlow_" + i);
            }

            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr vehCap = cplex.linearNumExpr();
                for (int i = 1; i <= n; i++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        vehCap.addTerm(ins.g(i, v, t), lambda[i][v][t]);
                    }
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
                cplex.addLe(m[t], ins.K, "VehicleCount_" + t);
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
                            for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                                mtz.addTerm(-ins.g(i, v, t), lambda[i][v][t]);
                            }
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
            }

            boolean solved = cplex.solve();
            return buildResult("ReformulationModel", solved, cplex, startNs);
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve ReformulationModel", e);
        }
    }

    private static boolean isPickupNode(int node, int n) {
        return node >= 1 && node <= n;
    }

    private static double holdingCostOnArc(Instance ins, int i, int v, int t) {
        double e = ins.e(i, v, t);
        if (v == 0) {
            // Align with corrected model.tex:
            // e_{i0t} should sum j=1..t-1, while Instance.e uses j=0..t-1.
            // Therefore subtract one h_i * I_{i0}.
            return e - ins.hi[i] * ins.Ii0[i];
        }
        return e;
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
