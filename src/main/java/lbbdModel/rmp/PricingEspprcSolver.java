package lbbdModel.rmp;

import instance.Instance;

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
        Search search;
        if (activeGlobalCustomers.length <= 62) {
            search = new LongMaskSearch(ins, activeGlobalCustomers, qLocal, dualUiLocal, dualU0, existingRouteKeys);
        } else {
            search = new BitSetSearch(ins, activeGlobalCustomers, qLocal, dualUiLocal, dualU0, existingRouteKeys);
        }
        search.run();
        if (search.bestRoute != null && search.bestReducedCost < -RC_EPS) {
            return new Result(true, search.bestReducedCost, search.bestRoute);
        }
        return new Result(false, search.bestReducedCost, null);
    }

    private abstract static class Search {
        final Instance ins;
        final int[] customers;      // local -> global
        final double[] q;
        final double[] dual;
        final double dualU0;
        final HashSet<String> existingKeys;

        double bestReducedCost = Double.POSITIVE_INFINITY;
        RouteColumn bestRoute = null;

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
        }

        abstract void run();

        final void considerCompleteRoute(AbstractLabel label) {
            int lastGlobalNode = customers[label.lastLocal];
            double fullCost = label.travelCost + ins.c[lastGlobalNode][ins.n + 1];
            double reducedCost = fullCost - label.dualGain - dualU0;
            if (reducedCost >= bestReducedCost - DOM_EPS) {
                return;
            }
            int[] globals = reconstructGlobals(label);
            String key = buildRouteKey(globals);
            if (existingKeys.contains(key)) {
                return;
            }
            bestReducedCost = reducedCost;
            bestRoute = new RouteColumn(globals, fullCost, label.load);
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
    }

    private abstract static class AbstractLabel {
        final int lastLocal;
        final int depth;
        final double load;
        final double travelCost; // depot -> ... -> last (without return arc)
        final double dualGain;   // sum duals for visited customers
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
            this.pred = pred;
        }
    }

    private static final class LongMaskSearch extends Search {
        private final long fullMask;
        private final ArrayDeque<LongLabel> queue;
        private final HashMap<LongStateKey, LongLabel> bestByState;

        LongMaskSearch(
                Instance ins,
                int[] customers,
                double[] q,
                double[] dual,
                double dualU0,
                HashSet<String> existingKeys
        ) {
            super(ins, customers, q, dual, dualU0, existingKeys);
            this.fullMask = (1L << customers.length) - 1L;
            this.queue = new ArrayDeque<LongLabel>();
            this.bestByState = new HashMap<LongStateKey, LongLabel>();
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
            LongStateKey key = new LongStateKey(label.lastLocal, label.visitedMask);
            LongLabel incumbent = bestByState.get(key);
            if (incumbent != null) {
                if (incumbent.travelCost <= label.travelCost + DOM_EPS) {
                    return false;
                }
            }
            bestByState.put(key, label);
            return true;
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

    private static final class LongStateKey {
        final int lastLocal;
        final long visitedMask;
        final int hash;

        LongStateKey(int lastLocal, long visitedMask) {
            this.lastLocal = lastLocal;
            this.visitedMask = visitedMask;
            int h = 31 * lastLocal + (int) (visitedMask ^ (visitedMask >>> 32));
            this.hash = h;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (!(o instanceof LongStateKey)) {
                return false;
            }
            LongStateKey other = (LongStateKey) o;
            return this.lastLocal == other.lastLocal && this.visitedMask == other.visitedMask;
        }

        @Override
        public int hashCode() {
            return hash;
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
                HashSet<String> existingKeys
        ) {
            super(ins, customers, q, dual, dualU0, existingKeys);
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
