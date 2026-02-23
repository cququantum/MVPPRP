package lbbdModel;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import instance.Instance;
import model.CplexConfig;
import model.SolveResult;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;

public final class LbbdReformulationSolver {

    private static final double CUT_EPS = 1e-6;
    private static final double RESULT_TOL = 1e-4;
    private static final int MAX_ITERATIONS = 5000;
    private static final boolean LOG_TO_CONSOLE = false;

    public SolveResult solve(Instance ins) {
        return solve(ins, Double.NaN, RESULT_TOL);
    }

    public SolveResult solve(Instance ins, double targetObjective, double targetTol) {
        long startNs = System.nanoTime();
        PeriodCvrpSubproblemSolver subproblemSolver = new PeriodCvrpSubproblemSolver();
        HashMap<String, PeriodCvrpSubproblemSolver.Result> subproblemCache =
                new HashMap<String, PeriodCvrpSubproblemSolver.Result>();
        HashMap<String, PeriodCvrpSubproblemSolver.DualCutResult> lpCutCache =
                new HashMap<String, PeriodCvrpSubproblemSolver.DualCutResult>();

        try (IloCplex cplex = new IloCplex()) {
            configure(cplex);
            MasterModel master = new MasterModel(cplex, ins);

            int iteration = 0;
            int feasibilityCuts = 0;
            int optimalityCuts = 0;
            int lpDualCuts = 0;
            int targetPruneCuts = 0;

            double bestUpperBound = Double.POSITIVE_INFINITY;
            double bestUpperBoundGapRef = Double.NaN;
            double bestLowerBound = Double.NEGATIVE_INFINITY;

            boolean allMasterOptimal = true;
            boolean allSubproblemsOptimal = true;

            while (iteration < MAX_ITERATIONS) {
                iteration++;

                boolean solved = cplex.solve();
                String masterStatus = safeStatus(cplex);
                if (!solved) {
                    return buildFailedResult(startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts, masterStatus);
                }
                if (!masterStatus.startsWith("Optimal")) {
                    allMasterOptimal = false;
                }

                MasterPoint point = master.extractPoint(ins);
                bestLowerBound = point.bmpObjective;

                double phi = 0.0;
                double[] phiByPeriod = new double[ins.l + 1];
                int infeasiblePeriod = -1;
                for (int t = 1; t <= ins.l; t++) {
                    String key = buildSubproblemCacheKey(point, t, ins.n);
                    PeriodCvrpSubproblemSolver.Result sub = subproblemCache.get(key);
                    if (sub == null) {
                        sub = subproblemSolver.solve(ins, t, point.qBar[t], point.zBar[t]);
                        if (sub.optimal || !sub.feasible) {
                            subproblemCache.put(key, sub);
                        }
                    }
                    if (!sub.feasible) {
                        infeasiblePeriod = t;
                        break;
                    }
                    if (!sub.optimal) {
                        allSubproblemsOptimal = false;
                    }
                    phiByPeriod[t] = sub.objective;
                    phi += sub.objective;

                    PeriodCvrpSubproblemSolver.DualCutResult lpCut = lpCutCache.get(key);
                    if (lpCut == null) {
                        lpCut = subproblemSolver.solveLpDualCut(ins, t, point.qBar[t], point.zBar[t]);
                        if (lpCut.feasible && lpCut.optimal) {
                            lpCutCache.put(key, lpCut);
                            master.addPeriodLpDualCut(t, lpCut);
                            lpDualCuts++;
                        }
                    }
                }

                if (infeasiblePeriod >= 0) {
                    System.out.println("[LBBD] iter=" + iteration
                            + " infeasiblePeriod=" + infeasiblePeriod
                            + " masterObj=" + fmt(point.bmpObjective)
                            + " cacheSize=" + subproblemCache.size()
                            + " lpCacheSize=" + lpCutCache.size()
                            + " lpDualCuts=" + lpDualCuts
                            + " targetPruneCuts=" + targetPruneCuts);
                    boolean addedStrong = master.addStrongFeasibilityCut(point, infeasiblePeriod);
                    if (!addedStrong) {
                        master.addNoGoodCut(point);
                    }
                    feasibilityCuts++;
                    continue;
                }

                double exactObjective = point.nonRouteCost + phi;
                if (exactObjective < bestUpperBound) {
                    bestUpperBound = exactObjective;
                    bestUpperBoundGapRef = point.bmpObjective;
                }

                double delta = phi - point.omegaSum;
                if (shouldLogIteration(iteration)) {
                    System.out.println("[LBBD] iter=" + iteration
                            + " masterObj=" + fmt(point.bmpObjective)
                            + " exactObj=" + fmt(exactObjective)
                            + " omegaSum=" + fmt(point.omegaSum)
                            + " phi=" + fmt(phi)
                            + " delta(phi-omega)=" + fmt(delta)
                            + " cacheSize=" + subproblemCache.size()
                            + " lpCacheSize=" + lpCutCache.size()
                            + " lpDualCuts=" + lpDualCuts
                            + " targetPruneCuts=" + targetPruneCuts);
                }

                if (!Double.isNaN(targetObjective) && Math.abs(exactObjective - targetObjective) <= targetTol) {
                    String status = "LBBD_TargetMatched(iter=" + iteration
                            + ",target=" + fmt(targetObjective)
                            + ",feasCuts=" + feasibilityCuts
                            + ",optCuts=" + optimalityCuts
                            + ",lpDualCuts=" + lpDualCuts
                            + ",targetPruneCuts=" + targetPruneCuts + ")";
                    double sec = elapsedSec(startNs);
                    return new SolveResult(
                            "LBBDReformulation",
                            true,
                            false,
                            status,
                            exactObjective,
                            point.bmpObjective,
                            relativeGap(exactObjective, point.bmpObjective),
                            sec
                    );
                }
                if (!Double.isNaN(targetObjective) && exactObjective < targetObjective - targetTol) {
                    throw new IllegalStateException(
                            "LBBD found objective better than reform target: exact="
                                    + exactObjective + ", target=" + targetObjective
                    );
                }

                boolean violatedOptimality = phi > point.omegaSum + CUT_EPS;
                if (violatedOptimality) {
                    master.addPeriodOptimalityCuts(point, phiByPeriod);
                    master.addOptimalityCut(point, phi);
                    optimalityCuts++;
                }

                if (!Double.isNaN(targetObjective) && exactObjective > targetObjective + targetTol) {
                    master.addNoGoodCut(point);
                    targetPruneCuts++;
                    continue;
                }

                if (violatedOptimality) {
                    continue;
                }

                double finalGap = relativeGap(exactObjective, point.bmpObjective);
                boolean optimal = allMasterOptimal && allSubproblemsOptimal && finalGap <= RESULT_TOL;
                String status = "LBBD_Converged(iter=" + iteration
                        + ",feasCuts=" + feasibilityCuts
                        + ",optCuts=" + optimalityCuts
                        + ",lpDualCuts=" + lpDualCuts
                        + ",targetPruneCuts=" + targetPruneCuts
                        + ",masterStatus=" + masterStatus + ")";

                double sec = elapsedSec(startNs);
                return new SolveResult(
                        "LBBDReformulation",
                        true,
                        optimal,
                        status,
                        exactObjective,
                        point.bmpObjective,
                        finalGap,
                        sec
                );
            }

            double sec = elapsedSec(startNs);
            boolean feasible = !Double.isInfinite(bestUpperBound) && !Double.isNaN(bestUpperBound);
            String status = "LBBD_MaxIterations(iter=" + iteration
                    + ",feasCuts=" + feasibilityCuts
                    + ",optCuts=" + optimalityCuts
                    + ",lpDualCuts=" + lpDualCuts
                    + ",targetPruneCuts=" + targetPruneCuts + ")";
            double objective = feasible ? bestUpperBound : Double.NaN;
            double bestBound = Double.isInfinite(bestLowerBound) ? Double.NaN : bestLowerBound;
            double gap = (feasible && !Double.isNaN(bestBound))
                    ? relativeGap(objective, bestBound)
                    : Double.NaN;
            if (feasible && !Double.isNaN(bestUpperBoundGapRef)) {
                gap = relativeGap(bestUpperBound, bestUpperBoundGapRef);
            }

            return new SolveResult("LBBDReformulation", feasible, false, status, objective, bestBound, gap, sec);
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve LBBDReformulation", e);
        }
    }

    private static boolean shouldLogIteration(int iteration) {
        return iteration <= 10 || iteration % 10 == 0;
    }

    private static String buildSubproblemCacheKey(MasterPoint point, int t, int n) {
        StringBuilder sb = new StringBuilder();
        sb.append(t).append(':');
        for (int i = 1; i <= n; i++) {
            if (i > 1) {
                sb.append(',');
            }
            sb.append(point.prevVisit[i][t]);
        }
        return sb.toString();
    }

    private static SolveResult buildFailedResult(
            long startNs,
            int iteration,
            int feasibilityCuts,
            int optimalityCuts,
            int lpDualCuts,
            String masterStatus
    ) {
        double sec = elapsedSec(startNs);
        String status = "LBBD_MasterSolveFailed(iter=" + iteration
                + ",feasCuts=" + feasibilityCuts
                + ",optCuts=" + optimalityCuts
                + ",lpDualCuts=" + lpDualCuts
                + ",masterStatus=" + masterStatus + ")";
        return new SolveResult("LBBDReformulation", false, false, status, Double.NaN, Double.NaN, Double.NaN, sec);
    }

    private static double elapsedSec(long startNs) {
        return (System.nanoTime() - startNs) / 1_000_000_000.0;
    }

    private static String fmt(double v) {
        return String.format(Locale.US, "%.6f", v);
    }

    private static double relativeGap(double upper, double lower) {
        if (Double.isNaN(upper) || Double.isNaN(lower)) {
            return Double.NaN;
        }
        double denom = Math.max(1.0, Math.abs(upper));
        return Math.abs(upper - lower) / denom;
    }

    private static String safeStatus(IloCplex cplex) {
        try {
            return cplex.getStatus().toString();
        } catch (IloException e) {
            return "Unknown";
        }
    }

    private static boolean isPickupNode(int node, int n) {
        return node >= 1 && node <= n;
    }

    private static double holdingCostOnArc(Instance ins, int i, int v, int t) {
        double e = ins.e(i, v, t);
        if (v == 0) {
            return e - ins.hi[i] * ins.Ii0[i];
        }
        return e;
    }

    private static void configure(IloCplex cplex) throws IloException {
        cplex.setParam(IloCplex.Param.TimeLimit, CplexConfig.TIME_LIMIT_SEC);
        cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, CplexConfig.MIP_GAP);
        if (!LOG_TO_CONSOLE) {
            cplex.setOut(null);
            cplex.setWarning(null);
        }
    }

    private static final class LambdaRef {
        final int i;
        final int v;
        final int t;
        final IloNumVar var;

        LambdaRef(int i, int v, int t, IloNumVar var) {
            this.i = i;
            this.v = v;
            this.t = t;
            this.var = var;
        }
    }

    private static final class MasterPoint {
        final double bmpObjective;
        final double nonRouteCost;
        final double omegaSum;
        final double[] omega;
        final double[][] qBar; // [t][i]
        final int[][] zBar;    // [t][i]
        final int[][] prevVisit; // [i][t]
        final boolean[][][] lambdaOne; // [i][v][t], t in 1..l only
        final int routeLambdaOneCount;

        MasterPoint(
                double bmpObjective,
                double nonRouteCost,
                double omegaSum,
                double[] omega,
                double[][] qBar,
                int[][] zBar,
                int[][] prevVisit,
                boolean[][][] lambdaOne,
                int routeLambdaOneCount
        ) {
            this.bmpObjective = bmpObjective;
            this.nonRouteCost = nonRouteCost;
            this.omegaSum = omegaSum;
            this.omega = omega;
            this.qBar = qBar;
            this.zBar = zBar;
            this.prevVisit = prevVisit;
            this.lambdaOne = lambdaOne;
            this.routeLambdaOneCount = routeLambdaOneCount;
        }
    }

    private static final class MasterModel {
        private final IloCplex cplex;
        private final Instance ins;

        private final IloNumVar[] omega;
        private final IloNumVar[] m;
        private IloNumVar[][] z;
        private final IloNumVar[][][] lambda;
        private final ArrayList<LambdaRef> routingLambdaRefs;

        private int feasCutCount = 0;
        private int noGoodCutCount = 0;
        private int optCutCount = 0;

        MasterModel(IloCplex cplex, Instance ins) throws IloException {
            this.cplex = cplex;
            this.ins = ins;
            this.omega = new IloNumVar[ins.l + 1];
            this.m = new IloNumVar[ins.l + 1];
            this.lambda = new IloNumVar[ins.n + 1][ins.l + 2][ins.l + 2];
            this.routingLambdaRefs = new ArrayList<LambdaRef>();
            build();
        }

        private void build() throws IloException {
            int n = ins.n;
            int l = ins.l;
            double[] minIncomingToSupplier = new double[n + 1];
            double[] minOutgoingFromSupplier = new double[n + 1];
            for (int j = 1; j <= n; j++) {
                double minIn = Double.POSITIVE_INFINITY;
                double minOut = Double.POSITIVE_INFINITY;
                for (int i = 0; i < ins.nodeCount; i++) {
                    if (i == j) {
                        continue;
                    }
                    if (ins.c[i][j] < minIn) {
                        minIn = ins.c[i][j];
                    }
                    if (ins.c[j][i] < minOut) {
                        minOut = ins.c[j][i];
                    }
                }
                minIncomingToSupplier[j] = minIn;
                minOutgoingFromSupplier[j] = minOut;
            }
            double minRoundTrip = Double.POSITIVE_INFINITY;
            for (int j = 1; j <= n; j++) {
                double rt = ins.c[0][j] + ins.c[j][n + 1];
                if (rt < minRoundTrip) {
                    minRoundTrip = rt;
                }
            }
            double loadBasedOmegaCoeff = (minRoundTrip > 0.0 && ins.Q > 0.0) ? (minRoundTrip / ins.Q) : 0.0;

            IloNumVar[] y = new IloNumVar[l + 1];
            IloNumVar[] p = new IloNumVar[l + 1];
            IloNumVar[] p0 = new IloNumVar[l + 1];
            IloNumVar[] i0 = new IloNumVar[l + 1];
            z = new IloNumVar[n + 1][l + 1];

            for (int t = 1; t <= l; t++) {
                y[t] = cplex.boolVar("y_" + t);
                p[t] = cplex.numVar(0.0, Double.MAX_VALUE, "p_" + t);
                p0[t] = cplex.numVar(0.0, Double.MAX_VALUE, "P0_" + t);
                i0[t] = cplex.numVar(0.0, Double.MAX_VALUE, "I0_" + t);
                m[t] = cplex.intVar(0, ins.K, "m_" + t);
                omega[t] = cplex.numVar(0.0, Double.MAX_VALUE, "omega_" + t);
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l; t++) {
                    z[i][t] = cplex.boolVar("z_" + i + "_" + t);
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l + 1; t++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        IloNumVar var = cplex.boolVar("lambda_" + i + "_" + v + "_" + t);
                        lambda[i][v][t] = var;
                        if (t <= l) {
                            routingLambdaRefs.add(new LambdaRef(i, v, t, var));
                        }
                    }
                }
            }

            IloLinearNumExpr obj = cplex.linearNumExpr();
            for (int t = 1; t <= l; t++) {
                obj.addTerm(1.0, omega[t]);
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

            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr vehCap = cplex.linearNumExpr();
                for (int i = 1; i <= n; i++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        vehCap.addTerm(ins.g(i, v, t), lambda[i][v][t]);
                    }
                }
                vehCap.addTerm(-ins.Q, m[t]);
                cplex.addLe(vehCap, 0.0, "VehicleCap_" + t);
                cplex.addLe(m[t], ins.K, "VehicleCount_" + t);

                IloLinearNumExpr lbIn = cplex.linearNumExpr();
                IloLinearNumExpr lbOut = cplex.linearNumExpr();
                IloLinearNumExpr lbLoad = cplex.linearNumExpr();
                IloLinearNumExpr lbM = cplex.linearNumExpr();
                lbIn.addTerm(1.0, omega[t]);
                lbOut.addTerm(1.0, omega[t]);
                lbLoad.addTerm(1.0, omega[t]);
                lbM.addTerm(1.0, omega[t]);
                lbM.addTerm(-minRoundTrip, m[t]);

                for (int i = 1; i <= n; i++) {
                    lbIn.addTerm(-minIncomingToSupplier[i], z[i][t]);
                    lbOut.addTerm(-minOutgoingFromSupplier[i], z[i][t]);
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        lbLoad.addTerm(-loadBasedOmegaCoeff * ins.g(i, v, t), lambda[i][v][t]);
                    }
                }
                cplex.addGe(lbIn, 0.0, "OmegaLB_In_" + t);
                cplex.addGe(lbOut, 0.0, "OmegaLB_Out_" + t);
                cplex.addGe(lbLoad, 0.0, "OmegaLB_Load_" + t);
                cplex.addGe(lbM, 0.0, "OmegaLB_M_" + t);
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
        }

        MasterPoint extractPoint(Instance ins) throws IloException {
            int n = ins.n;
            int l = ins.l;

            double[] omegaVal = new double[l + 1];
            double omegaSum = 0.0;
            for (int t = 1; t <= l; t++) {
                omegaVal[t] = cplex.getValue(omega[t]);
                omegaSum += omegaVal[t];
            }

            boolean[][][] lambdaOne = new boolean[n + 1][l + 1][l + 1];
            int[][] prevVisit = new int[n + 1][l + 1];
            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l; t++) {
                    prevVisit[i][t] = -1;
                }
            }
            double[][] qBar = new double[l + 1][n + 1];
            int[][] zBar = new int[l + 1][n + 1];

            int routeLambdaOneCount = 0;
            for (int idx = 0; idx < routingLambdaRefs.size(); idx++) {
                LambdaRef ref = routingLambdaRefs.get(idx);
                double value = cplex.getValue(ref.var);
                if (value > 0.5) {
                    if (!lambdaOne[ref.i][ref.v][ref.t]) {
                        lambdaOne[ref.i][ref.v][ref.t] = true;
                        routeLambdaOneCount++;
                    }
                    if (prevVisit[ref.i][ref.t] >= 0 && prevVisit[ref.i][ref.t] != ref.v) {
                        throw new IllegalStateException(
                                "Multiple active lambda predecessors for i=" + ref.i + ", t=" + ref.t
                        );
                    }
                    prevVisit[ref.i][ref.t] = ref.v;
                    zBar[ref.t][ref.i] = 1;
                    qBar[ref.t][ref.i] += ins.g(ref.i, ref.v, ref.t);
                }
            }

            double bmpObjective = cplex.getObjValue();
            double nonRouteCost = bmpObjective - omegaSum;

            return new MasterPoint(
                    bmpObjective,
                    nonRouteCost,
                    omegaSum,
                    omegaVal,
                    qBar,
                    zBar,
                    prevVisit,
                    lambdaOne,
                    routeLambdaOneCount
            );
        }

        boolean addStrongFeasibilityCut(MasterPoint point, int t) throws IloException {
            IloLinearNumExpr expr = cplex.linearNumExpr();
            int visited = 0;

            for (int i = 1; i <= ins.n; i++) {
                int v = point.prevVisit[i][t];
                if (v >= 0) {
                    expr.addTerm(-1.0, lambda[i][v][t]);
                    visited++;
                }
            }
            if (visited == 0) {
                return false;
            }

            feasCutCount++;
            cplex.addGe(expr, 1.0 - visited, "LBBD_StrongFeasCut_" + feasCutCount + "_t" + t);
            return true;
        }

        void addNoGoodCut(MasterPoint point) throws IloException {
            IloLinearNumExpr expr = cplex.linearNumExpr();
            for (int idx = 0; idx < routingLambdaRefs.size(); idx++) {
                LambdaRef ref = routingLambdaRefs.get(idx);
                if (point.lambdaOne[ref.i][ref.v][ref.t]) {
                    expr.addTerm(-1.0, ref.var);
                } else {
                    expr.addTerm(1.0, ref.var);
                }
            }
            double rhs = 1.0 - point.routeLambdaOneCount;
            noGoodCutCount++;
            cplex.addGe(expr, rhs, "LBBD_NoGoodCut_" + noGoodCutCount);
        }

        void addOptimalityCut(MasterPoint point, double phi) throws IloException {
            IloLinearNumExpr expr = cplex.linearNumExpr();
            for (int t = 1; t <= ins.l; t++) {
                expr.addTerm(1.0, omega[t]);
            }
            for (int idx = 0; idx < routingLambdaRefs.size(); idx++) {
                LambdaRef ref = routingLambdaRefs.get(idx);
                if (point.lambdaOne[ref.i][ref.v][ref.t]) {
                    expr.addTerm(-phi, ref.var);
                } else {
                    expr.addTerm(phi, ref.var);
                }
            }
            double rhs = phi - phi * point.routeLambdaOneCount;
            optCutCount++;
            cplex.addGe(expr, rhs, "LBBD_OptCut_" + optCutCount);
        }

        void addPeriodOptimalityCuts(MasterPoint point, double[] phiByPeriod) throws IloException {
            for (int t = 1; t <= ins.l; t++) {
                double phiT = phiByPeriod[t];
                if (phiT <= CUT_EPS) {
                    continue;
                }

                IloLinearNumExpr expr = cplex.linearNumExpr();
                expr.addTerm(1.0, omega[t]);

                int onesCount = 0;
                for (int idx = 0; idx < routingLambdaRefs.size(); idx++) {
                    LambdaRef ref = routingLambdaRefs.get(idx);
                    if (ref.t != t) {
                        continue;
                    }
                    if (point.lambdaOne[ref.i][ref.v][ref.t]) {
                        expr.addTerm(-phiT, ref.var);
                        onesCount++;
                    } else {
                        expr.addTerm(phiT, ref.var);
                    }
                }

                double rhs = phiT - phiT * onesCount;
                optCutCount++;
                cplex.addGe(expr, rhs, "LBBD_PeriodOptCut_" + optCutCount + "_t" + t);
            }
        }

        void addPeriodLpDualCut(int t, PeriodCvrpSubproblemSolver.DualCutResult cut) throws IloException {
            IloLinearNumExpr expr = cplex.linearNumExpr();
            expr.addTerm(1.0, omega[t]);

            for (int i = 1; i <= ins.n; i++) {
                if (Math.abs(cut.coeffZ[i]) > 1e-9) {
                    expr.addTerm(-cut.coeffZ[i], z[i][t]);
                }
                if (Math.abs(cut.coeffQ[i]) > 1e-9) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        expr.addTerm(-cut.coeffQ[i] * ins.g(i, v, t), lambda[i][v][t]);
                    }
                }
            }

            optCutCount++;
            cplex.addGe(expr, cut.constant, "LBBD_LpDualCut_" + optCutCount + "_t" + t);
        }

    }
}
