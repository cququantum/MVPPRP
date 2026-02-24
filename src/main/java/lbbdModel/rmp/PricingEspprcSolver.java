package lbbdModel.rmp;

import instance.Instance;

import java.util.HashSet;

/**
 * Exact elementary pricing via DFS enumeration with capacity pruning.
 * This is exact (not heuristic), but exponential in the number of active customers.
 */
public final class PricingEspprcSolver {
    private static final double RC_EPS = 1e-8;

    public static final class Result {
        public final boolean foundNegativeColumn;
        public final double bestReducedCost;
        public final RouteColumn route;

        public Result(boolean foundNegativeColumn, double bestReducedCost, RouteColumn route) {
            this.foundNegativeColumn = foundNegativeColumn;
            this.bestReducedCost = bestReducedCost;
            this.route = route;
        }
    }

    public Result findBestNegativeRoute(
            Instance ins,
            int[] activeGlobalCustomers,
            double[] qLocal,
            double[] dualUiLocal,
            double dualU0,
            HashSet<String> existingRouteKeys
    ) {
        if (activeGlobalCustomers.length == 0) {
            return new Result(false, 0.0, null);
        }
        Search search = new Search(ins, activeGlobalCustomers, qLocal, dualUiLocal, dualU0, existingRouteKeys);
        search.run();
        if (search.bestRoute != null && search.bestReducedCost < -RC_EPS) {
            return new Result(true, search.bestReducedCost, search.bestRoute);
        }
        return new Result(false, search.bestReducedCost, null);
    }

    private static final class Search {
        private final Instance ins;
        private final int[] customers;      // local -> global
        private final double[] q;
        private final double[] dual;
        private final double dualU0;
        private final HashSet<String> existingKeys;
        private final boolean[] used;
        private final int[] pathLocal;

        private final double[] minFromDepot;
        private final double[] minToDepot;

        private double bestReducedCost = Double.POSITIVE_INFINITY;
        private RouteColumn bestRoute = null;

        Search(
                Instance ins,
                int[] customers,
                double[] q,
                double[] dual,
                double dualU0,
                HashSet<String> existingKeys
        ) {
            this.ins = ins;
            this.customers = customers;
            this.q = q;
            this.dual = dual;
            this.dualU0 = dualU0;
            this.existingKeys = existingKeys;
            this.used = new boolean[customers.length];
            this.pathLocal = new int[customers.length];
            this.minFromDepot = new double[customers.length];
            this.minToDepot = new double[customers.length];
            for (int k = 0; k < customers.length; k++) {
                int g = customers[k];
                minFromDepot[k] = ins.c[0][g];
                minToDepot[k] = ins.c[g][ins.n + 1];
            }
        }

        void run() {
            for (int start = 0; start < customers.length; start++) {
                if (q[start] > ins.Q + 1e-9) {
                    continue;
                }
                used[start] = true;
                pathLocal[0] = start;
                int g = customers[start];
                double travel = minFromDepot[start];
                double dualGain = dual[start];
                double load = q[start];
                evaluateAndRecurse(1, g, travel, dualGain, load);
                used[start] = false;
            }
        }

        private void evaluateAndRecurse(int depth, int lastGlobalNode, double travel, double dualGain, double load) {
            double reducedCost = travel + ins.c[lastGlobalNode][ins.n + 1] - dualGain - dualU0;
            if (reducedCost < bestReducedCost - 1e-12) {
                String key = buildCurrentRouteKey(depth);
                if (!existingKeys.contains(key)) {
                    bestReducedCost = reducedCost;
                    bestRoute = buildRoute(depth, travel + ins.c[lastGlobalNode][ins.n + 1], load);
                }
            }

            // Simple pruning: even if we append more nodes, the route must eventually return to depot.
            // Lower bound is current reduced cost minus any future dual gains + nonnegative travel;
            // without a strong admissible bound, we only prune on capacity.
            for (int nxt = 0; nxt < customers.length; nxt++) {
                if (used[nxt]) {
                    continue;
                }
                double newLoad = load + q[nxt];
                if (newLoad > ins.Q + 1e-9) {
                    continue;
                }
                used[nxt] = true;
                pathLocal[depth] = nxt;
                int nextGlobal = customers[nxt];
                double newTravel = travel + ins.c[lastGlobalNode][nextGlobal];
                double newDualGain = dualGain + dual[nxt];
                evaluateAndRecurse(depth + 1, nextGlobal, newTravel, newDualGain, newLoad);
                used[nxt] = false;
            }
        }

        private String buildCurrentRouteKey(int depth) {
            StringBuilder sb = new StringBuilder();
            for (int idx = 0; idx < depth; idx++) {
                if (idx > 0) {
                    sb.append('-');
                }
                sb.append(customers[pathLocal[idx]]);
            }
            return sb.toString();
        }

        private RouteColumn buildRoute(int depth, double fullCost, double load) {
            int[] globals = new int[depth];
            for (int idx = 0; idx < depth; idx++) {
                globals[idx] = customers[pathLocal[idx]];
            }
            return new RouteColumn(globals, fullCost, load);
        }
    }
}
