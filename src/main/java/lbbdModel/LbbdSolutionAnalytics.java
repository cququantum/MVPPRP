package lbbdModel;

import instance.Instance;

public final class LbbdSolutionAnalytics {
    private static final double EPS = 1e-9;

    private LbbdSolutionAnalytics() {
    }

    public static double[][] computeSupplierEndInventory(Instance ins, double[][] qBar, int[][] zBar) {
        double[][] supplierEndInventory = new double[ins.n + 1][ins.l + 1];
        double[] current = new double[ins.n + 1];
        for (int i = 1; i <= ins.n; i++) {
            current[i] = ins.Ii0[i];
        }

        for (int t = 1; t <= ins.l; t++) {
            for (int i = 1; i <= ins.n; i++) {
                double beforePickup = current[i] + ins.s[i][t];
                boolean visited = zBar != null
                        && t < zBar.length
                        && zBar[t] != null
                        && i < zBar[t].length
                        && zBar[t][i] != 0;
                if (visited) {
                    current[i] = 0.0;
                } else {
                    current[i] = beforePickup;
                }
                supplierEndInventory[i][t] = current[i];
            }
        }

        return supplierEndInventory;
    }

    public static double computeRouteCost(double[] phiByPeriod) {
        if (phiByPeriod == null) {
            return Double.NaN;
        }
        double total = 0.0;
        for (int t = 1; t < phiByPeriod.length; t++) {
            total += phiByPeriod[t];
        }
        return total;
    }

    public static double computeProductionCost(Instance ins, double[] y, double[] p) {
        if (y == null || p == null) {
            return Double.NaN;
        }
        double total = 0.0;
        for (int t = 1; t < y.length && t < p.length; t++) {
            total += ins.u * p[t] + ins.f * y[t];
        }
        return total;
    }

    public static double computeInventoryCost(
            Instance ins,
            double[] i0,
            double[] p0,
            double[][] qBar,
            int[][] zBar
    ) {
        if (i0 == null || p0 == null) {
            return Double.NaN;
        }
        double total = 0.0;
        for (int t = 1; t < i0.length && t < p0.length; t++) {
            total += ins.h0 * i0[t] + ins.hp * p0[t];
        }
        double[][] supplierEndInventory = computeSupplierEndInventory(ins, qBar, zBar);
        for (int i = 1; i <= ins.n; i++) {
            for (int t = 1; t <= ins.l; t++) {
                total += ins.hi[i] * supplierEndInventory[i][t];
            }
        }
        return total;
    }

    public static double computeAverageVisitFrequency(Instance ins, int[][] zBar) {
        if (zBar == null) {
            return Double.NaN;
        }
        double totalVisits = 0.0;
        for (int t = 1; t < zBar.length; t++) {
            if (zBar[t] == null) {
                continue;
            }
            for (int i = 1; i <= ins.n && i < zBar[t].length; i++) {
                totalVisits += zBar[t][i];
            }
        }
        return totalVisits / Math.max(1, ins.n);
    }

    public static double computeAveragePickupPerVisit(Instance ins, double[][] qBar, int[][] zBar) {
        if (qBar == null || zBar == null) {
            return Double.NaN;
        }
        double totalPickup = 0.0;
        double totalVisits = 0.0;
        for (int t = 1; t < qBar.length && t < zBar.length; t++) {
            if (qBar[t] != null) {
                for (int i = 1; i <= ins.n && i < qBar[t].length; i++) {
                    totalPickup += qBar[t][i];
                }
            }
            if (zBar[t] != null) {
                for (int i = 1; i <= ins.n && i < zBar[t].length; i++) {
                    totalVisits += zBar[t][i];
                }
            }
        }
        if (totalVisits <= EPS) {
            return 0.0;
        }
        return totalPickup / totalVisits;
    }

    public static int countProductionPeriods(double[] y) {
        if (y == null) {
            return 0;
        }
        int count = 0;
        for (int t = 1; t < y.length; t++) {
            if (y[t] > 0.5) {
                count++;
            }
        }
        return count;
    }

    public static double computeAverageProductionBatch(double[] p, double[] y) {
        if (p == null || y == null) {
            return Double.NaN;
        }
        int periods = countProductionPeriods(y);
        if (periods == 0) {
            return 0.0;
        }
        double total = 0.0;
        for (int t = 1; t < p.length; t++) {
            total += p[t];
        }
        return total / periods;
    }

    public static double computeRawInventoryAverage(Instance ins, double[] i0) {
        if (i0 == null) {
            return Double.NaN;
        }
        double total = 0.0;
        for (int t = 1; t < i0.length; t++) {
            total += i0[t];
        }
        return total / Math.max(1, ins.l);
    }

    public static double computeFinishedInventoryAverage(Instance ins, double[] p0) {
        if (p0 == null) {
            return Double.NaN;
        }
        double total = 0.0;
        for (int t = 1; t < p0.length; t++) {
            total += p0[t];
        }
        return total / Math.max(1, ins.l);
    }

    public static double computeSupplierInventoryAverage(Instance ins, double[][] supplierEndInventory) {
        if (supplierEndInventory == null) {
            return Double.NaN;
        }
        double total = 0.0;
        for (int i = 1; i <= ins.n; i++) {
            if (i >= supplierEndInventory.length || supplierEndInventory[i] == null) {
                continue;
            }
            for (int t = 1; t <= ins.l && t < supplierEndInventory[i].length; t++) {
                total += supplierEndInventory[i][t];
            }
        }
        return total / Math.max(1, ins.n * ins.l);
    }
}
