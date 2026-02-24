package lbbdModel.rmp;

import instance.Instance;

import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Exact ESPPRC pricing by label expansion (label-setting style).
 *
 * Compared with plain DFS permutation enumeration, this implementation merges permutations
 * that reach the same exact state (last customer, visited set), which is an exact dominance
 * rule and substantially reduces repeated work.
 */
public final class PricingEspprcSolver {
    private static final double RC_EPS = 1e-8;
    private static final double DOM_EPS = 1e-12;

    public static final class Result {
        public final boolean foundNegativeColumn;
        public final double bestReducedCost;
        public final RouteColumn route;
        public final ArrayList<RouteColumn> routes;

        public Result(boolean foundNegativeColumn, double bestReducedCost, RouteColumn route, ArrayList<RouteColumn> routes) {
            this.foundNegativeColumn = foundNegativeColumn;
            this.bestReducedCost = bestReducedCost;
            this.route = route;
            this.routes = routes;
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
        return findNegativeRoutes(ins, activeGlobalCustomers, qLocal, dualUiLocal, dualU0, existingRouteKeys, 1);
    }

    public Result findNegativeRoutes(
            Instance ins,
            int[] activeGlobalCustomers,
            double[] qLocal,
            double[] dualUiLocal,
            double dualU0,
            HashSet<String> existingRouteKeys,
            int maxColumns
    ) {
        if (activeGlobalCustomers.length == 0) {
            return new Result(false, 0.0, null, new ArrayList<RouteColumn>());
        }
        Search search;
        if (activeGlobalCustomers.length <= 62) {
            search = new LongMaskSearch(ins, activeGlobalCustomers, qLocal, dualUiLocal, dualU0, existingRouteKeys, maxColumns);
        } else {
            search = new BitSetSearch(ins, activeGlobalCustomers, qLocal, dualUiLocal, dualU0, existingRouteKeys, maxColumns);
        }
        search.run();
        ArrayList<RouteColumn> routes = search.getCollectedRoutes();
        if (search.bestRoute != null && search.bestReducedCost < -RC_EPS) {
            return new Result(true, search.bestReducedCost, search.bestRoute, routes);
        }
        return new Result(false, search.bestReducedCost, null, routes);
    }

    private abstract static class Search {
        final Instance ins;
        final int[] customers;      // local -> global
        final double[] q;
        final double[] dual;
        final double dualU0;
        final HashSet<String> existingKeys;
        final int maxColumns;
        final ArrayList<CollectedColumn> collected;
        final HashSet<String> collectedKeys;

        double bestReducedCost = Double.POSITIVE_INFINITY;
        RouteColumn bestRoute = null;

        Search(
                Instance ins,
                int[] customers,
                double[] q,
                double[] dual,
                double dualU0,
                HashSet<String> existingKeys,
                int maxColumns
        ) {
            this.ins = ins;
            this.customers = customers;
            this.q = q;
            this.dual = dual;
            this.dualU0 = dualU0;
            this.existingKeys = existingKeys;
            this.maxColumns = Math.max(0, maxColumns);
            this.collected = new ArrayList<CollectedColumn>();
            this.collectedKeys = new HashSet<String>();
        }

        abstract void run();

        final ArrayList<RouteColumn> getCollectedRoutes() {
            ArrayList<RouteColumn> out = new ArrayList<RouteColumn>(collected.size());
            for (int i = 0; i < collected.size(); i++) {
                out.add(collected.get(i).route);
            }
            return out;
        }

        final void considerCompleteRoute(AbstractLabel label) {
            int lastGlobalNode = customers[label.lastLocal];
            double fullCost = label.travelCost + ins.c[lastGlobalNode][ins.n + 1];
            double reducedCost = fullCost - label.dualGain - dualU0;
            boolean canImproveBest = reducedCost < bestReducedCost - DOM_EPS;
            boolean canCollect = maxColumns > 0 && reducedCost < -RC_EPS
                    && (collected.size() < maxColumns
                    || reducedCost < collected.get(collected.size() - 1).reducedCost - DOM_EPS);
            if (!canImproveBest && !canCollect) {
                return;
            }
            int[] globals = reconstructGlobals(label);
            String key = buildRouteKey(globals);
            if (existingKeys.contains(key)) {
                return;
            }
            RouteColumn route = new RouteColumn(globals, fullCost, label.load);
            if (canImproveBest) {
                bestReducedCost = reducedCost;
                bestRoute = route;
            }
            if (maxColumns > 0 && reducedCost < -RC_EPS && !collectedKeys.contains(key)) {
                maybeCollect(new CollectedColumn(reducedCost, route, key));
            }
        }

        private int[] reconstructGlobals(AbstractLabel label) {
            int[] globals = new int[label.depth];
            AbstractLabel cur = label;
            for (int pos = label.depth - 1; pos >= 0; pos--) {
                globals[pos] = customers[cur.lastLocal];
                cur = cur.pred;
            }
            return globals;
        }

        private String buildRouteKey(int[] globals) {
            StringBuilder sb = new StringBuilder();
            for (int idx = 0; idx < globals.length; idx++) {
                if (idx > 0) {
                    sb.append('-');
                }
                sb.append(globals[idx]);
            }
            return sb.toString();
        }

        private void maybeCollect(CollectedColumn candidate) {
            if (collected.size() >= maxColumns) {
                CollectedColumn worst = collected.get(collected.size() - 1);
                if (candidate.reducedCost >= worst.reducedCost - DOM_EPS) {
                    return;
                }
                collectedKeys.remove(worst.key);
                collected.remove(collected.size() - 1);
            }
            int pos = collected.size();
            while (pos > 0 && candidate.reducedCost < collected.get(pos - 1).reducedCost - DOM_EPS) {
                pos--;
            }
            collected.add(pos, candidate);
            collectedKeys.add(candidate.key);
        }
    }

    private static final class CollectedColumn {
        final double reducedCost;
        final RouteColumn route;
        final String key;

        CollectedColumn(double reducedCost, RouteColumn route, String key) {
            this.reducedCost = reducedCost;
            this.route = route;
            this.key = key;
        }
    }

    private abstract static class AbstractLabel {
        final int lastLocal;
        final int depth;
        final double load;
        final double travelCost; // depot -> ... -> last (without return arc)
        final double dualGain;   // sum duals for visited customers
        final double partialReducedCost; // travelCost - dualGain
        final AbstractLabel pred;

        AbstractLabel(
                int lastLocal,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                AbstractLabel pred
        ) {
            this.lastLocal = lastLocal;
            this.depth = depth;
            this.load = load;
            this.travelCost = travelCost;
            this.dualGain = dualGain;
            this.partialReducedCost = travelCost - dualGain;
            this.pred = pred;
        }
    }

    private static final class LongMaskSearch extends Search {
        private final long fullMask;
        private final long[] localBit;
        private final double[] returnToDepot;
        private final int[] positiveDualOrder;
        private final ArrayDeque<LongLabel> queue;
        private final HashMap<Long, LongLabel>[] bestByMaskAtLast;

        LongMaskSearch(
                Instance ins,
                int[] customers,
                double[] q,
                double[] dual,
                double dualU0,
                HashSet<String> existingKeys,
                int maxColumns
        ) {
            super(ins, customers, q, dual, dualU0, existingKeys, maxColumns);
            this.fullMask = (1L << customers.length) - 1L;
            this.localBit = new long[customers.length];
            this.returnToDepot = new double[customers.length];
            for (int k = 0; k < customers.length; k++) {
                this.localBit[k] = 1L << k;
                this.returnToDepot[k] = ins.c[customers[k]][ins.n + 1];
            }
            this.positiveDualOrder = buildPositiveDualOrder(dual);
            this.queue = new ArrayDeque<LongLabel>();
            @SuppressWarnings("unchecked")
            HashMap<Long, LongLabel>[] maps = (HashMap<Long, LongLabel>[]) new HashMap[customers.length];
            for (int k = 0; k < customers.length; k++) {
                maps[k] = new HashMap<Long, LongLabel>();
            }
            this.bestByMaskAtLast = maps;
        }

        @Override
        void run() {
            for (int start = 0; start < customers.length; start++) {
                if (q[start] > ins.Q + 1e-9) {
                    continue;
                }
                long mask = (1L << start);
                int g = customers[start];
                LongLabel startLabel = new LongLabel(
                        start, 1, q[start], ins.c[0][g], dual[start], null, mask
                );
                if (register(startLabel)) {
                    queue.addLast(startLabel);
                }
            }

            while (!queue.isEmpty()) {
                LongLabel cur = queue.pollFirst();
                considerCompleteRoute(cur);
                if (cannotBeatBest(cur)) {
                    continue;
                }

                long remaining = fullMask & ~cur.visitedMask;
                while (remaining != 0L) {
                    long bit = remaining & -remaining;
                    int nxt = Long.numberOfTrailingZeros(bit);
                    remaining ^= bit;

                    double newLoad = cur.load + q[nxt];
                    if (newLoad > ins.Q + 1e-9) {
                        continue;
                    }

                    int curGlobal = customers[cur.lastLocal];
                    int nxtGlobal = customers[nxt];
                    LongLabel child = new LongLabel(
                            nxt,
                            cur.depth + 1,
                            newLoad,
                            cur.travelCost + ins.c[curGlobal][nxtGlobal],
                            cur.dualGain + dual[nxt],
                            cur,
                            cur.visitedMask | bit
                    );
                    if (register(child)) {
                        queue.addLast(child);
                    }
                }
            }
        }

        private boolean register(LongLabel label) {
            HashMap<Long, LongLabel> exactMap = bestByMaskAtLast[label.lastLocal];
            LongLabel incumbent = exactMap.get(label.visitedMask);
            if (incumbent != null) {
                if (incumbent.travelCost <= label.travelCost + DOM_EPS) {
                    return false;
                }
            }

            if (isDominatedBySubset(label, exactMap)) {
                return false;
            }

            exactMap.put(label.visitedMask, label);
            return true;
        }

        private boolean isDominatedBySubset(LongLabel label, HashMap<Long, LongLabel> exactMap) {
            if (label.depth <= 1) {
                return false;
            }
            long lastBit = localBit[label.lastLocal];
            long free = label.visitedMask ^ lastBit;
            long subFree = free;
            while (true) {
                long subsetMask = subFree | lastBit;
                LongLabel subset = exactMap.get(subsetMask);
                if (subset != null && subset.partialReducedCost <= label.partialReducedCost + DOM_EPS) {
                    return true;
                }
                if (subFree == 0L) {
                    break;
                }
                subFree = (subFree - 1L) & free;
            }
            return false;
        }

        private boolean cannotBeatBest(LongLabel label) {
            if (bestReducedCost == Double.POSITIVE_INFINITY) {
                return false;
            }
            double optimisticExtraDual = upperBoundExtraDualGain(label);
            double lowerBoundFinalReducedCost =
                    label.travelCost + returnToDepot[label.lastLocal] - label.dualGain - dualU0 - optimisticExtraDual;
            return lowerBoundFinalReducedCost >= bestReducedCost - DOM_EPS;
        }

        private double upperBoundExtraDualGain(LongLabel label) {
            if (positiveDualOrder.length == 0) {
                return 0.0;
            }
            double remCap = ins.Q - label.load;
            if (remCap <= 1e-9) {
                return 0.0;
            }
            double ub = 0.0;
            long visited = label.visitedMask;
            for (int p = 0; p < positiveDualOrder.length; p++) {
                int idx = positiveDualOrder[p];
                long bit = localBit[idx];
                if ((visited & bit) != 0L) {
                    continue;
                }
                double qq = q[idx];
                if (qq <= remCap + 1e-9) {
                    ub += dual[idx];
                    remCap -= qq;
                } else {
                    ub += dual[idx] * (remCap / qq);
                    break;
                }
                if (remCap <= 1e-9) {
                    break;
                }
            }
            return ub;
        }

        private int[] buildPositiveDualOrder(double[] dual) {
            int m = dual.length;
            int[] tmp = new int[m];
            int cnt = 0;
            for (int i = 0; i < m; i++) {
                if (dual[i] > 1e-12) {
                    tmp[cnt++] = i;
                }
            }
            if (cnt == 0) {
                return new int[0];
            }
            int[] arr = new int[cnt];
            System.arraycopy(tmp, 0, arr, 0, cnt);
            // insertion sort (m is small in pricing)
            for (int i = 1; i < cnt; i++) {
                int key = arr[i];
                double keyRatio = dual[key] / Math.max(1e-12, q[key]);
                int j = i - 1;
                while (j >= 0) {
                    double ratioJ = dual[arr[j]] / Math.max(1e-12, q[arr[j]]);
                    if (ratioJ >= keyRatio - 1e-15) {
                        break;
                    }
                    arr[j + 1] = arr[j];
                    j--;
                }
                arr[j + 1] = key;
            }
            return arr;
        }
    }

    private static final class LongLabel extends AbstractLabel {
        final long visitedMask;

        LongLabel(
                int lastLocal,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                AbstractLabel pred,
                long visitedMask
        ) {
            super(lastLocal, depth, load, travelCost, dualGain, pred);
            this.visitedMask = visitedMask;
        }
    }

    private static final class BitSetSearch extends Search {
        private final int customerCount;
        private final ArrayDeque<BitSetLabel> queue;
        private final HashMap<BitSetStateKey, BitSetLabel> bestByState;

        BitSetSearch(
                Instance ins,
                int[] customers,
                double[] q,
                double[] dual,
                double dualU0,
                HashSet<String> existingKeys,
                int maxColumns
        ) {
            super(ins, customers, q, dual, dualU0, existingKeys, maxColumns);
            this.customerCount = customers.length;
            this.queue = new ArrayDeque<BitSetLabel>();
            this.bestByState = new HashMap<BitSetStateKey, BitSetLabel>();
        }

        @Override
        void run() {
            for (int start = 0; start < customerCount; start++) {
                if (q[start] > ins.Q + 1e-9) {
                    continue;
                }
                BitSet visited = new BitSet(customerCount);
                visited.set(start);
                int g = customers[start];
                BitSetLabel startLabel = new BitSetLabel(
                        start, 1, q[start], ins.c[0][g], dual[start], null, visited
                );
                if (register(startLabel)) {
                    queue.addLast(startLabel);
                }
            }

            while (!queue.isEmpty()) {
                BitSetLabel cur = queue.pollFirst();
                considerCompleteRoute(cur);

                int nxt = cur.visited.nextClearBit(0);
                while (nxt >= 0 && nxt < customerCount) {
                    double newLoad = cur.load + q[nxt];
                    if (newLoad <= ins.Q + 1e-9) {
                        int curGlobal = customers[cur.lastLocal];
                        int nxtGlobal = customers[nxt];
                        BitSet nextVisited = (BitSet) cur.visited.clone();
                        nextVisited.set(nxt);
                        BitSetLabel child = new BitSetLabel(
                                nxt,
                                cur.depth + 1,
                                newLoad,
                                cur.travelCost + ins.c[curGlobal][nxtGlobal],
                                cur.dualGain + dual[nxt],
                                cur,
                                nextVisited
                        );
                        if (register(child)) {
                            queue.addLast(child);
                        }
                    }
                    nxt = cur.visited.nextClearBit(nxt + 1);
                }
            }
        }

        private boolean register(BitSetLabel label) {
            BitSetStateKey key = new BitSetStateKey(label.lastLocal, label.visited);
            BitSetLabel incumbent = bestByState.get(key);
            if (incumbent != null) {
                if (incumbent.travelCost <= label.travelCost + DOM_EPS) {
                    return false;
                }
            }
            bestByState.put(key, label);
            return true;
        }
    }

    private static final class BitSetLabel extends AbstractLabel {
        final BitSet visited;

        BitSetLabel(
                int lastLocal,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                AbstractLabel pred,
                BitSet visited
        ) {
            super(lastLocal, depth, load, travelCost, dualGain, pred);
            this.visited = visited;
        }
    }

    private static final class BitSetStateKey {
        final int lastLocal;
        final BitSet visited;
        final int hash;

        BitSetStateKey(int lastLocal, BitSet visited) {
            this.lastLocal = lastLocal;
            this.visited = visited;
            this.hash = 31 * lastLocal + visited.hashCode();
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (!(o instanceof BitSetStateKey)) {
                return false;
            }
            BitSetStateKey other = (BitSetStateKey) o;
            return this.lastLocal == other.lastLocal && this.visited.equals(other.visited);
        }

        @Override
        public int hashCode() {
            return hash;
        }
    }
}
