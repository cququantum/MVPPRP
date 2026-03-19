package lbbdModel;

import model.SolveResult;

public final class LbbdDetailedResult {
    public final SolveResult summary;
    public final int iterations;
    public final int feasibilityCuts;
    public final int optimalityCuts;
    public final int lpDualCuts;
    public final double[] phiByPeriod;
    public final double[][] qBar;
    public final int[][] zBar;
    public final int[][] prevVisit;
    public final double[] y;
    public final double[] p;
    public final double[] i0;
    public final double[] p0;
    public final double routeCost;
    public final double inventoryCost;
    public final double productionCost;

    public LbbdDetailedResult(
            SolveResult summary,
            int iterations,
            int feasibilityCuts,
            int optimalityCuts,
            int lpDualCuts,
            double[] phiByPeriod,
            double[][] qBar,
            int[][] zBar,
            int[][] prevVisit,
            double[] y,
            double[] p,
            double[] i0,
            double[] p0,
            double routeCost,
            double inventoryCost,
            double productionCost
    ) {
        this.summary = summary;
        this.iterations = iterations;
        this.feasibilityCuts = feasibilityCuts;
        this.optimalityCuts = optimalityCuts;
        this.lpDualCuts = lpDualCuts;
        this.phiByPeriod = copy(phiByPeriod);
        this.qBar = copy(qBar);
        this.zBar = copy(zBar);
        this.prevVisit = copy(prevVisit);
        this.y = copy(y);
        this.p = copy(p);
        this.i0 = copy(i0);
        this.p0 = copy(p0);
        this.routeCost = routeCost;
        this.inventoryCost = inventoryCost;
        this.productionCost = productionCost;
    }

    public LbbdDetailedResult copy() {
        return new LbbdDetailedResult(
                summary,
                iterations,
                feasibilityCuts,
                optimalityCuts,
                lpDualCuts,
                phiByPeriod,
                qBar,
                zBar,
                prevVisit,
                y,
                p,
                i0,
                p0,
                routeCost,
                inventoryCost,
                productionCost
        );
    }

    private static double[] copy(double[] src) {
        if (src == null) {
            return null;
        }
        double[] dst = new double[src.length];
        System.arraycopy(src, 0, dst, 0, src.length);
        return dst;
    }

    private static int[] copy(int[] src) {
        if (src == null) {
            return null;
        }
        int[] dst = new int[src.length];
        System.arraycopy(src, 0, dst, 0, src.length);
        return dst;
    }

    private static double[][] copy(double[][] src) {
        if (src == null) {
            return null;
        }
        double[][] dst = new double[src.length][];
        for (int i = 0; i < src.length; i++) {
            dst[i] = copy(src[i]);
        }
        return dst;
    }

    private static int[][] copy(int[][] src) {
        if (src == null) {
            return null;
        }
        int[][] dst = new int[src.length][];
        for (int i = 0; i < src.length; i++) {
            dst[i] = copy(src[i]);
        }
        return dst;
    }
}
