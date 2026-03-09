package lbbdModel;

import instance.Instance;
import model.CplexConfig;

public final class PeriodCvrpSubproblemSolver implements AutoCloseable {
    private static final double EPS = 1e-9;

    private final BranchAndPriceCvrpSolver branchAndPriceSolver = new BranchAndPriceCvrpSolver();

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
        // No-op: this wrapper currently delegates all work to BranchAndPriceCvrpSolver.
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
}
