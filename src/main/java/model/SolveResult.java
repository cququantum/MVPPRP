package model;

import java.util.Locale;

public final class SolveResult {
    public final String modelName;
    public final boolean feasible;
    public final boolean optimal;
    public final String status;
    public final double objective;
    public final double bestBound;
    public final double mipGap;
    public final double solveTimeSec;

    public SolveResult(
            String modelName,
            boolean feasible,
            boolean optimal,
            String status,
            double objective,
            double bestBound,
            double mipGap,
            double solveTimeSec
    ) {
        this.modelName = modelName;
        this.feasible = feasible;
        this.optimal = optimal;
        this.status = status;
        this.objective = objective;
        this.bestBound = bestBound;
        this.mipGap = mipGap;
        this.solveTimeSec = solveTimeSec;
    }

    public String toSummaryLine() {
        StringBuilder sb = new StringBuilder();
        sb.append(modelName)
          .append(" | status=").append(status)
          .append(" | feasible=").append(feasible)
          .append(" | optimal=").append(optimal);

        if (!Double.isNaN(objective)) {
            sb.append(" | obj=").append(format(objective));
        }
        if (!Double.isNaN(bestBound)) {
            sb.append(" | bestBound=").append(format(bestBound));
        }
        if (!Double.isNaN(mipGap)) {
            sb.append(" | gap=").append(format(mipGap));
        }
        sb.append(" | timeSec=").append(format(solveTimeSec));
        return sb.toString();
    }

    private static String format(double v) {
        return String.format(Locale.US, "%.6f", v);
    }
}
