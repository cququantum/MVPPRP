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
    private static final int EXACT_DP_MAX_N = 15;

    private final BranchAndPriceCvrpSolver branchAndPriceSolver = new BranchAndPriceCvrpSolver();
    // Legacy exact solvers are kept only as internal validation or fallback building blocks.
    private final HashMap<Integer, ExactPeriodModel> exactModels = new HashMap<Integer, ExactPeriodModel>();
    private Instance subsetDpInstance = null;
    private ExactSubsetDpOracle subsetDpOracle = null;

    public static final class Result {
        public final boolean feasible;
        public final boolean optimal;
        public final boolean provenInfeasible;
        public final String status;
        public final double objective;

        public Result(boolean feasible, boolean optimal, boolean provenInfeasible, String status, double objective) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.provenInfeasible = provenInfeasible;
            this.status = status;
            this.objective = objective;
        }
    }

    public Result solve(Instance ins, int t, double[] qBar, int[] zBar) {
        return solve(ins, t, qBar, zBar, CplexConfig.TIME_LIMIT_SEC);
    }

    public Result solve(Instance ins, int t, double[] qBar, int[] zBar, double timeLimitSec) {
        validateInput(ins, t, qBar, zBar);
        if (timeLimitSec <= 0.0) {
            return new Result(false, false, false, "Subproblem_TimeLimit", Double.NaN);
        }

        double totalPickup = totalPickup(ins, qBar);
        int visitCount = visitCount(ins, zBar);
        if (visitCount == 0 && totalPickup <= EPS) {
            return new Result(true, true, false, "TrivialZero", 0.0);
        }
        if (totalPickup > ins.Q * ins.K + EPS) {
            return new Result(false, false, true, "InfeasibleByTotalCapacity", Double.NaN);
        }
        BranchAndPriceCvrpSolver.Result bp = branchAndPriceSolver.solve(ins, t, qBar, zBar, timeLimitSec);
        return new Result(bp.feasible, bp.optimal, bp.provenInfeasible, bp.status, bp.objective);
    }

    @Override
    public void close() {
        for (ExactPeriodModel model : exactModels.values()) {
            model.close();
        }
        exactModels.clear();
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
                    return new Result(false, false, true, "InfeasibleByPositivePickupWithoutVisit", Double.NaN);
                }
            }
            if (activeMask == 0) {
                return new Result(true, true, false, "TrivialZeroDP", 0.0);
            }

            buildSubsetLoads(qBar);

            // If any visited customer alone exceeds capacity, infeasible.
            int bits = activeMask;
            while (bits != 0) {
                int bit = bits & -bits;
                int idx = Integer.numberOfTrailingZeros(bit);
                if (qBar[idx + 1] > ins.Q + EPS) {
                    return new Result(false, false, true, "InfeasibleBySingleCustomerCapacity", Double.NaN);
                }
                bits ^= bit;
            }

            clearMemo();
            double best = solvePartition(activeMask, ins.K);
            if (best >= INF / 2) {
                return new Result(false, false, true, "InfeasibleByRoutePartitionDP", Double.NaN);
            }
            return new Result(true, true, false, "OptimalDP", best);
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
                    if (i == j || !isValidArc(i, j, n)) {
                        continue;
                    }
                    x[i][j] = cplex.boolVar("x_" + i + "_" + j + "_" + t);
                }
            }

            IloLinearNumExpr obj = cplex.linearNumExpr();
            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (x[i][j] == null) {
                        continue;
                    }
                    obj.addTerm(ins.c[i][j], x[i][j]);
                }
            }
            cplex.addMinimize(obj);

            IloLinearNumExpr vehCapExpr = cplex.linearNumExpr();
            vehCapExpr.addTerm(-ins.Q, m);
            vehCap = cplex.addLe(vehCapExpr, 0.0, "TotalPickupBound_" + t);

            for (int i = 1; i <= n; i++) {
                IloLinearNumExpr out = cplex.linearNumExpr();
                IloLinearNumExpr in = cplex.linearNumExpr();
                for (int j = 0; j < nodeCount; j++) {
                    if (x[i][j] != null) {
                        out.addTerm(1.0, x[i][j]);
                    }
                    if (x[j][i] != null) {
                        in.addTerm(1.0, x[j][i]);
                    }
                }
                cplex.addEq(out, in, "FlowBalance_" + i + "_" + t);
            }

            IloLinearNumExpr depart = cplex.linearNumExpr();
            for (int j = 1; j <= n; j++) {
                if (x[0][j] != null) {
                    depart.addTerm(1.0, x[0][j]);
                }
            }
            cplex.addEq(depart, m, "DepartFactory_" + t);

            IloLinearNumExpr back = cplex.linearNumExpr();
            for (int i = 1; i <= n; i++) {
                if (x[i][n + 1] != null) {
                    back.addTerm(1.0, x[i][n + 1]);
                }
            }
            cplex.addEq(back, m, "ReturnFactory_" + t);

            for (int j = 1; j <= n; j++) {
                IloLinearNumExpr incoming = cplex.linearNumExpr();
                for (int i = 0; i < nodeCount; i++) {
                    if (x[i][j] != null) {
                        incoming.addTerm(1.0, x[i][j]);
                    }
                }
                visitLink[j] = cplex.addEq(incoming, 0.0, "VisitLink_" + j + "_" + t);
            }
            cplex.addLe(m, ins.K, "VehicleCount_" + t);

            for (int i = 0; i < nodeCount; i++) {
                for (int j = 0; j < nodeCount; j++) {
                    if (x[i][j] == null) {
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
                return new Result(false, false, false, status, Double.NaN);
            }
            return new Result(true, optimal, false, status, cplex.getObjValue());
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
                    if (mtzRows[i][j] == null) {
                        continue;
                    }
                    mtzRows[i][j].setLB(rhs);
                }
            }
        }

        private static boolean isValidArc(int i, int j, int n) {
            int returnDepot = n + 1;
            if (i == j) return false;
            if (i == returnDepot) return false;
            if (i == 0 && (j < 1 || j > n)) return false;
            if (j == 0) return false;
            if (j == returnDepot && (i < 1 || i > n)) return false;
            return true;
        }

        void close() {
            cplex.end();
        }
    }

}
