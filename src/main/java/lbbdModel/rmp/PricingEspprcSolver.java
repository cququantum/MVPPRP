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
    private static final int SUPERSET_CLEANUP_MAX_EXTRA_BITS = 7;

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

    public static final class ScenarioResult {
        public final boolean foundNegativeColumn;
        public final double bestReducedCost;
        public final ScenarioRouteColumn route;
        public final ArrayList<ScenarioRouteColumn> routes;

        public ScenarioResult(
                boolean foundNegativeColumn,
                double bestReducedCost,
                ScenarioRouteColumn route,
                ArrayList<ScenarioRouteColumn> routes
        ) {
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

    public ScenarioResult findNegativeScenarioRoutes(
            Instance ins,
            int[] supplierByScenario,
            double[] demandByScenario,
            double[] dualByScenario,
            double dualU0,
            HashSet<String> existingRouteKeys,
            int maxColumns
    ) {
        int m = supplierByScenario.length;
        if (m == 0) {
            return new ScenarioResult(false, 0.0, null, new ArrayList<ScenarioRouteColumn>());
        }
        if (demandByScenario.length != m || dualByScenario.length != m) {
            throw new IllegalArgumentException("Scenario arrays length mismatch");
        }

        ScenarioSearch search;
        if (m <= 62 && ins.n <= 62) {
            search = new ScenarioLongSearch(
                    ins,
                    supplierByScenario,
                    demandByScenario,
                    dualByScenario,
                    dualU0,
                    existingRouteKeys,
                    maxColumns
            );
        } else {
            search = new ScenarioBitSetSearch(
                    ins,
                    supplierByScenario,
                    demandByScenario,
                    dualByScenario,
                    dualU0,
                    existingRouteKeys,
                    maxColumns
            );
        }
        search.run();
        ArrayList<ScenarioRouteColumn> routes = search.getCollectedRoutes();
        if (search.bestRoute != null && search.bestReducedCost < -RC_EPS) {
            return new ScenarioResult(true, search.bestReducedCost, search.bestRoute, routes);
        }
        return new ScenarioResult(false, search.bestReducedCost, null, routes);
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
                if (!isCurrentBest(cur)) {
                    continue;
                }
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

            if (isDominatedByProperSubset(label, exactMap)) {
                return false;
            }

            exactMap.put(label.visitedMask, label);
            removeDominatedSupersetsNearTail(label, exactMap);
            return true;
        }

        private boolean isCurrentBest(LongLabel label) {
            return bestByMaskAtLast[label.lastLocal].get(label.visitedMask) == label;
        }

        private boolean isDominatedByProperSubset(LongLabel label, HashMap<Long, LongLabel> exactMap) {
            if (label.depth <= 1) {
                return false;
            }
            long lastBit = localBit[label.lastLocal];
            long free = label.visitedMask ^ lastBit;
            long subFree = free;
            while (subFree != 0L) {
                long subsetMask = subFree | lastBit;
                LongLabel subset = exactMap.get(subsetMask);
                if (subset != null && subset.partialReducedCost <= label.partialReducedCost + DOM_EPS) {
                    return true;
                }
                subFree = (subFree - 1L) & free;
            }
            return false;
        }

        private void removeDominatedSupersetsNearTail(LongLabel label, HashMap<Long, LongLabel> exactMap) {
            if (label.depth <= 1) {
                return;
            }
            long mask = label.visitedMask;
            long extraPool = fullMask & ~mask;
            int extraCount = Long.bitCount(extraPool);
            if (extraCount == 0 || extraCount > SUPERSET_CLEANUP_MAX_EXTRA_BITS) {
                return;
            }
            double rc = label.partialReducedCost;
            long add = extraPool;
            while (add != 0L) {
                LongLabel other = exactMap.get(mask | add);
                if (other != null && rc <= other.partialReducedCost + DOM_EPS) {
                    exactMap.remove(mask | add);
                }
                add = (add - 1L) & extraPool;
            }
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
        private final double[] returnToDepot;
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
            this.returnToDepot = new double[customerCount];
            for (int k = 0; k < customerCount; k++) {
                this.returnToDepot[k] = ins.c[customers[k]][ins.n + 1];
            }
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
                if (!isCurrentBest(cur)) {
                    continue;
                }
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

        private boolean isCurrentBest(BitSetLabel label) {
            BitSetStateKey key = new BitSetStateKey(label.lastLocal, label.visited);
            return bestByState.get(key) == label;
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

    private abstract static class ScenarioSearch {
        final Instance ins;
        final int[] supplierByScenario;
        final double[] demandByScenario;
        final double[] dualByScenario;
        final double dualU0;
        final HashSet<String> existingRouteKeys;
        final int maxColumns;
        final ArrayList<ScenarioCollectedColumn> collected;
        final HashSet<String> collectedKeys;

        double bestReducedCost = Double.POSITIVE_INFINITY;
        ScenarioRouteColumn bestRoute = null;

        ScenarioSearch(
                Instance ins,
                int[] supplierByScenario,
                double[] demandByScenario,
                double[] dualByScenario,
                double dualU0,
                HashSet<String> existingRouteKeys,
                int maxColumns
        ) {
            this.ins = ins;
            this.supplierByScenario = supplierByScenario;
            this.demandByScenario = demandByScenario;
            this.dualByScenario = dualByScenario;
            this.dualU0 = dualU0;
            this.existingRouteKeys = existingRouteKeys;
            this.maxColumns = Math.max(0, maxColumns);
            this.collected = new ArrayList<ScenarioCollectedColumn>();
            this.collectedKeys = new HashSet<String>();
        }

        abstract void run();

        final ArrayList<ScenarioRouteColumn> getCollectedRoutes() {
            ArrayList<ScenarioRouteColumn> out = new ArrayList<ScenarioRouteColumn>(collected.size());
            for (int i = 0; i < collected.size(); i++) {
                out.add(collected.get(i).route);
            }
            return out;
        }

        final void considerCompleteRoute(AbstractScenarioLabel label) {
            int lastSupplier = supplierByScenario[label.lastScenario];
            double fullCost = label.travelCost + ins.c[lastSupplier][ins.n + 1];
            double reducedCost = fullCost - label.dualGain - dualU0;
            boolean canImproveBest = reducedCost < bestReducedCost - DOM_EPS;
            boolean canCollect = maxColumns > 0 && reducedCost < -RC_EPS
                    && (collected.size() < maxColumns
                    || reducedCost < collected.get(collected.size() - 1).reducedCost - DOM_EPS);
            if (!canImproveBest && !canCollect) {
                return;
            }
            int[] scenarioSeq = reconstructScenarioSeq(label);
            String key = buildScenarioKey(scenarioSeq);
            if (existingRouteKeys.contains(key)) {
                return;
            }
            ScenarioRouteColumn route = new ScenarioRouteColumn(scenarioSeq, fullCost, label.load);
            if (canImproveBest) {
                bestReducedCost = reducedCost;
                bestRoute = route;
            }
            if (maxColumns > 0 && reducedCost < -RC_EPS && !collectedKeys.contains(key)) {
                maybeCollect(new ScenarioCollectedColumn(reducedCost, route, key));
            }
        }

        private int[] reconstructScenarioSeq(AbstractScenarioLabel label) {
            int[] seq = new int[label.depth];
            AbstractScenarioLabel cur = label;
            for (int p = label.depth - 1; p >= 0; p--) {
                seq[p] = cur.lastScenario;
                cur = cur.pred;
            }
            return seq;
        }

        private static String buildScenarioKey(int[] seq) {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < seq.length; i++) {
                if (i > 0) {
                    sb.append('-');
                }
                sb.append('s').append(seq[i]);
            }
            return sb.toString();
        }

        private void maybeCollect(ScenarioCollectedColumn candidate) {
            if (collected.size() >= maxColumns) {
                ScenarioCollectedColumn worst = collected.get(collected.size() - 1);
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

    private static final class ScenarioCollectedColumn {
        final double reducedCost;
        final ScenarioRouteColumn route;
        final String key;

        ScenarioCollectedColumn(double reducedCost, ScenarioRouteColumn route, String key) {
            this.reducedCost = reducedCost;
            this.route = route;
            this.key = key;
        }
    }

    private abstract static class AbstractScenarioLabel {
        final int lastScenario;
        final int depth;
        final double load;
        final double travelCost;
        final double dualGain;
        final AbstractScenarioLabel pred;

        AbstractScenarioLabel(
                int lastScenario,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                AbstractScenarioLabel pred
        ) {
            this.lastScenario = lastScenario;
            this.depth = depth;
            this.load = load;
            this.travelCost = travelCost;
            this.dualGain = dualGain;
            this.pred = pred;
        }
    }

    private static final class ScenarioLongSearch extends ScenarioSearch {
        private final int scenarioCount;
        private final long fullMask;
        private final long[] scenarioBit;
        private final long[] supplierBitByScenario;
        private final double[] returnToDepot;
        private final ArrayDeque<ScenarioLongLabel> queue;
        private final HashMap<Long, ScenarioLongLabel>[] bestByMaskAtLast;

        ScenarioLongSearch(
                Instance ins,
                int[] supplierByScenario,
                double[] demandByScenario,
                double[] dualByScenario,
                double dualU0,
                HashSet<String> existingRouteKeys,
                int maxColumns
        ) {
            super(ins, supplierByScenario, demandByScenario, dualByScenario, dualU0, existingRouteKeys, maxColumns);
            this.scenarioCount = supplierByScenario.length;
            this.fullMask = (1L << scenarioCount) - 1L;
            this.scenarioBit = new long[scenarioCount];
            this.supplierBitByScenario = new long[scenarioCount];
            this.returnToDepot = new double[scenarioCount];
            for (int s = 0; s < scenarioCount; s++) {
                scenarioBit[s] = 1L << s;
                supplierBitByScenario[s] = 1L << (supplierByScenario[s] - 1);
                returnToDepot[s] = ins.c[supplierByScenario[s]][ins.n + 1];
            }
            this.queue = new ArrayDeque<ScenarioLongLabel>();
            @SuppressWarnings("unchecked")
            HashMap<Long, ScenarioLongLabel>[] arr = (HashMap<Long, ScenarioLongLabel>[]) new HashMap[scenarioCount];
            for (int s = 0; s < scenarioCount; s++) {
                arr[s] = new HashMap<Long, ScenarioLongLabel>();
            }
            this.bestByMaskAtLast = arr;
        }

        @Override
        void run() {
            for (int s = 0; s < scenarioCount; s++) {
                if (demandByScenario[s] > ins.Q + 1e-9) {
                    continue;
                }
                long visitedMask = scenarioBit[s];
                long supplierMask = supplierBitByScenario[s];
                int supplier = supplierByScenario[s];
                ScenarioLongLabel start = new ScenarioLongLabel(
                        s,
                        1,
                        demandByScenario[s],
                        ins.c[0][supplier],
                        dualByScenario[s],
                        null,
                        visitedMask,
                        supplierMask
                );
                if (register(start)) {
                    queue.addLast(start);
                }
            }

            while (!queue.isEmpty()) {
                ScenarioLongLabel cur = queue.pollFirst();
                if (!isCurrentBest(cur)) {
                    continue;
                }
                considerCompleteRoute(cur);

                long remaining = fullMask & ~cur.visitedMask;
                while (remaining != 0L) {
                    long bit = remaining & -remaining;
                    int nxt = Long.numberOfTrailingZeros(bit);
                    remaining ^= bit;

                    long nxtSupplierBit = supplierBitByScenario[nxt];
                    if ((cur.supplierMask & nxtSupplierBit) != 0L) {
                        continue;
                    }

                    double newLoad = cur.load + demandByScenario[nxt];
                    if (newLoad > ins.Q + 1e-9) {
                        continue;
                    }

                    int curSupplier = supplierByScenario[cur.lastScenario];
                    int nxtSupplier = supplierByScenario[nxt];
                    ScenarioLongLabel child = new ScenarioLongLabel(
                            nxt,
                            cur.depth + 1,
                            newLoad,
                            cur.travelCost + ins.c[curSupplier][nxtSupplier],
                            cur.dualGain + dualByScenario[nxt],
                            cur,
                            cur.visitedMask | bit,
                            cur.supplierMask | nxtSupplierBit
                    );
                    if (register(child)) {
                        queue.addLast(child);
                    }
                }
            }
        }

        private boolean register(ScenarioLongLabel label) {
            HashMap<Long, ScenarioLongLabel> map = bestByMaskAtLast[label.lastScenario];
            ScenarioLongLabel incumbent = map.get(label.visitedMask);
            if (incumbent != null && incumbent.travelCost <= label.travelCost + DOM_EPS) {
                return false;
            }
            map.put(label.visitedMask, label);
            return true;
        }

        private boolean isCurrentBest(ScenarioLongLabel label) {
            return bestByMaskAtLast[label.lastScenario].get(label.visitedMask) == label;
        }
    }

    private static final class ScenarioLongLabel extends AbstractScenarioLabel {
        final long visitedMask;
        final long supplierMask;

        ScenarioLongLabel(
                int lastScenario,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                AbstractScenarioLabel pred,
                long visitedMask,
                long supplierMask
        ) {
            super(lastScenario, depth, load, travelCost, dualGain, pred);
            this.visitedMask = visitedMask;
            this.supplierMask = supplierMask;
        }
    }

    private static final class ScenarioBitSetSearch extends ScenarioSearch {
        private final int scenarioCount;
        private final double[] returnToDepot;
        private final ArrayDeque<ScenarioBitSetLabel> queue;
        private final HashMap<ScenarioBitSetStateKey, ScenarioBitSetLabel> bestByState;

        ScenarioBitSetSearch(
                Instance ins,
                int[] supplierByScenario,
                double[] demandByScenario,
                double[] dualByScenario,
                double dualU0,
                HashSet<String> existingRouteKeys,
                int maxColumns
        ) {
            super(ins, supplierByScenario, demandByScenario, dualByScenario, dualU0, existingRouteKeys, maxColumns);
            this.scenarioCount = supplierByScenario.length;
            this.returnToDepot = new double[scenarioCount];
            for (int s = 0; s < scenarioCount; s++) {
                this.returnToDepot[s] = ins.c[supplierByScenario[s]][ins.n + 1];
            }
            this.queue = new ArrayDeque<ScenarioBitSetLabel>();
            this.bestByState = new HashMap<ScenarioBitSetStateKey, ScenarioBitSetLabel>();
        }

        @Override
        void run() {
            for (int s = 0; s < scenarioCount; s++) {
                if (demandByScenario[s] > ins.Q + 1e-9) {
                    continue;
                }
                BitSet visitedScenarios = new BitSet(scenarioCount);
                visitedScenarios.set(s);
                BitSet visitedSuppliers = new BitSet(ins.n + 1);
                visitedSuppliers.set(supplierByScenario[s]);
                int supplier = supplierByScenario[s];
                ScenarioBitSetLabel start = new ScenarioBitSetLabel(
                        s,
                        1,
                        demandByScenario[s],
                        ins.c[0][supplier],
                        dualByScenario[s],
                        null,
                        visitedScenarios,
                        visitedSuppliers
                );
                if (register(start)) {
                    queue.addLast(start);
                }
            }

            while (!queue.isEmpty()) {
                ScenarioBitSetLabel cur = queue.pollFirst();
                if (!isCurrentBest(cur)) {
                    continue;
                }
                considerCompleteRoute(cur);

                int nxt = cur.visitedScenarios.nextClearBit(0);
                while (nxt >= 0 && nxt < scenarioCount) {
                    int nxtSupplier = supplierByScenario[nxt];
                    if (cur.visitedSuppliers.get(nxtSupplier)) {
                        nxt = cur.visitedScenarios.nextClearBit(nxt + 1);
                        continue;
                    }

                    double newLoad = cur.load + demandByScenario[nxt];
                    if (newLoad > ins.Q + 1e-9) {
                        nxt = cur.visitedScenarios.nextClearBit(nxt + 1);
                        continue;
                    }

                    BitSet nextVisitedScenarios = (BitSet) cur.visitedScenarios.clone();
                    nextVisitedScenarios.set(nxt);
                    BitSet nextVisitedSuppliers = (BitSet) cur.visitedSuppliers.clone();
                    nextVisitedSuppliers.set(nxtSupplier);

                    int curSupplier = supplierByScenario[cur.lastScenario];
                    ScenarioBitSetLabel child = new ScenarioBitSetLabel(
                            nxt,
                            cur.depth + 1,
                            newLoad,
                            cur.travelCost + ins.c[curSupplier][nxtSupplier],
                            cur.dualGain + dualByScenario[nxt],
                            cur,
                            nextVisitedScenarios,
                            nextVisitedSuppliers
                    );
                    if (register(child)) {
                        queue.addLast(child);
                    }
                    nxt = cur.visitedScenarios.nextClearBit(nxt + 1);
                }
            }
        }

        private boolean register(ScenarioBitSetLabel label) {
            ScenarioBitSetStateKey key = new ScenarioBitSetStateKey(label.lastScenario, label.visitedScenarios);
            ScenarioBitSetLabel incumbent = bestByState.get(key);
            if (incumbent != null && incumbent.travelCost <= label.travelCost + DOM_EPS) {
                return false;
            }
            bestByState.put(key, label);
            return true;
        }

        private boolean isCurrentBest(ScenarioBitSetLabel label) {
            ScenarioBitSetStateKey key = new ScenarioBitSetStateKey(label.lastScenario, label.visitedScenarios);
            return bestByState.get(key) == label;
        }
    }

    private static final class ScenarioBitSetLabel extends AbstractScenarioLabel {
        final BitSet visitedScenarios;
        final BitSet visitedSuppliers;

        ScenarioBitSetLabel(
                int lastScenario,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                AbstractScenarioLabel pred,
                BitSet visitedScenarios,
                BitSet visitedSuppliers
        ) {
            super(lastScenario, depth, load, travelCost, dualGain, pred);
            this.visitedScenarios = visitedScenarios;
            this.visitedSuppliers = visitedSuppliers;
        }
    }

    private static final class ScenarioBitSetStateKey {
        final int lastScenario;
        final BitSet visitedScenarios;
        final int hash;

        ScenarioBitSetStateKey(int lastScenario, BitSet visitedScenarios) {
            this.lastScenario = lastScenario;
            this.visitedScenarios = visitedScenarios;
            this.hash = 31 * lastScenario + visitedScenarios.hashCode();
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (!(o instanceof ScenarioBitSetStateKey)) {
                return false;
            }
            ScenarioBitSetStateKey other = (ScenarioBitSetStateKey) o;
            return this.lastScenario == other.lastScenario
                    && this.visitedScenarios.equals(other.visitedScenarios);
        }

        @Override
        public int hashCode() {
            return hash;
        }
    }
}
