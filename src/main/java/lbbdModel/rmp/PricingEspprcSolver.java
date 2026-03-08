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
    private static final int SUBSET_DOMINANCE_MAX_EXTRA_BITS = 8;
    private static final int SUPERSET_CLEANUP_MAX_EXTRA_BITS = 7;
    private static final long BIDIRECTIONAL_RELAXED_TIME_BUDGET_NS = 300_000_000L;

    private static boolean isFinite(double value) {
        return !Double.isNaN(value) && !Double.isInfinite(value);
    }

    private static boolean shouldEnumerateSubsetDominance(long optionalBits) {
        return Long.bitCount(optionalBits) <= SUBSET_DOMINANCE_MAX_EXTRA_BITS;
    }

    /**
     * Optional branch-node route restrictions for branch-and-price.
     * For a customer pair (i,j):
     * - "separate" is modeled by forbidding routes that contain both.
     * - "together" is modeled by allowing only routes that contain both or neither.
     *
     * The arrays are indexed by the local customer index used in the current pricing call.
     */
    public static final class RoutePricingConstraints {
        final long[] requiredWithLocalLong;
        final long[] forbiddenWithLocalLong;
        final BitSet[] requiredWithLocalBits;
        final BitSet[] forbiddenWithLocalBits;
        final boolean[] forbiddenStartLocal;
        final boolean[] forbiddenReturnLocal;
        final boolean[][] forbiddenArcLocal;
        final double[] returnArcDualLocal;
        final double[][] arcDualLocal;
        final boolean hasArcDuals;

        private RoutePricingConstraints(
                long[] requiredWithLocalLong,
                long[] forbiddenWithLocalLong,
                BitSet[] requiredWithLocalBits,
                BitSet[] forbiddenWithLocalBits,
                boolean[] forbiddenStartLocal,
                boolean[] forbiddenReturnLocal,
                boolean[][] forbiddenArcLocal,
                double[] returnArcDualLocal,
                double[][] arcDualLocal
        ) {
            this.requiredWithLocalLong = requiredWithLocalLong;
            this.forbiddenWithLocalLong = forbiddenWithLocalLong;
            this.requiredWithLocalBits = requiredWithLocalBits;
            this.forbiddenWithLocalBits = forbiddenWithLocalBits;
            this.forbiddenStartLocal = forbiddenStartLocal;
            this.forbiddenReturnLocal = forbiddenReturnLocal;
            this.forbiddenArcLocal = forbiddenArcLocal;
            this.returnArcDualLocal = returnArcDualLocal;
            this.arcDualLocal = arcDualLocal;
            this.hasArcDuals = returnArcDualLocal != null || arcDualLocal != null;
        }

        public static RoutePricingConstraints fromLongMasks(long[] requiredWithLocalLong, long[] forbiddenWithLocalLong) {
            return new RoutePricingConstraints(requiredWithLocalLong, forbiddenWithLocalLong, null, null,
                    null, null, null, null, null);
        }

        public static RoutePricingConstraints fromBitSets(BitSet[] requiredWithLocalBits, BitSet[] forbiddenWithLocalBits) {
            return new RoutePricingConstraints(null, null, requiredWithLocalBits, forbiddenWithLocalBits,
                    null, null, null, null, null);
        }

        public static RoutePricingConstraints create(
                long[] requiredWithLocalLong,
                long[] forbiddenWithLocalLong,
                BitSet[] requiredWithLocalBits,
                BitSet[] forbiddenWithLocalBits,
                boolean[] forbiddenStartLocal,
                boolean[] forbiddenReturnLocal,
                boolean[][] forbiddenArcLocal,
                double[] returnArcDualLocal,
                double[][] arcDualLocal
        ) {
            return new RoutePricingConstraints(
                    requiredWithLocalLong,
                    forbiddenWithLocalLong,
                    requiredWithLocalBits,
                    forbiddenWithLocalBits,
                    forbiddenStartLocal,
                    forbiddenReturnLocal,
                    forbiddenArcLocal,
                    returnArcDualLocal,
                    arcDualLocal
            );
        }
    }

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
        return findNegativeRoutes(ins, activeGlobalCustomers, qLocal, dualUiLocal, dualU0, existingRouteKeys, maxColumns, null);
    }

    public Result findNegativeRoutes(
            Instance ins,
            int[] activeGlobalCustomers,
            double[] qLocal,
            double[] dualUiLocal,
            double dualU0,
            HashSet<String> existingRouteKeys,
            int maxColumns,
            RoutePricingConstraints constraints
    ) {
        if (activeGlobalCustomers.length == 0) {
            return new Result(false, 0.0, null, new ArrayList<RouteColumn>());
        }
        if (activeGlobalCustomers.length <= 62) {
            Search bidirectional = new BidirectionalLongSearch(
                    ins, activeGlobalCustomers, qLocal, dualUiLocal, dualU0, existingRouteKeys, maxColumns, constraints
            );
            bidirectional.run();
            ArrayList<RouteColumn> routes = bidirectional.getCollectedRoutes();
            if (bidirectional.bestRoute != null && bidirectional.bestReducedCost < -RC_EPS) {
                return new Result(true, bidirectional.bestReducedCost, bidirectional.bestRoute, routes);
            }
            // Keep the existing exact single-direction search as a no-negative-column confirmation fallback.
            Search fallback = new LongMaskSearch(
                    ins, activeGlobalCustomers, qLocal, dualUiLocal, dualU0, existingRouteKeys, maxColumns, constraints
            );
            fallback.run();
            ArrayList<RouteColumn> fallbackRoutes = fallback.getCollectedRoutes();
            if (fallback.bestRoute != null && fallback.bestReducedCost < -RC_EPS) {
                return new Result(true, fallback.bestReducedCost, fallback.bestRoute, fallbackRoutes);
            }
            return new Result(false, fallback.bestReducedCost, null, fallbackRoutes);
        }
        Search search = new BitSetSearch(ins, activeGlobalCustomers, qLocal, dualUiLocal, dualU0, existingRouteKeys, maxColumns, constraints);
        search.run();
        ArrayList<RouteColumn> routes = search.getCollectedRoutes();
        if (search.bestRoute != null && search.bestReducedCost < -RC_EPS) {
            return new Result(true, search.bestReducedCost, search.bestRoute, routes);
        }
        return new Result(false, search.bestReducedCost, null, routes);
    }

    public Result findNegativeRoutesRelaxedOnly(
            Instance ins,
            int[] activeGlobalCustomers,
            double[] qLocal,
            double[] dualUiLocal,
            double dualU0,
            HashSet<String> existingRouteKeys,
            int maxColumns,
            RoutePricingConstraints constraints
    ) {
        return findNegativeRoutesBidirectionalRelaxedOnly(
                ins,
                activeGlobalCustomers,
                qLocal,
                dualUiLocal,
                dualU0,
                existingRouteKeys,
                maxColumns,
                constraints
        );
    }

    private Result findNegativeRoutesBidirectionalRelaxedOnly(
            Instance ins,
            int[] activeGlobalCustomers,
            double[] qLocal,
            double[] dualUiLocal,
            double dualU0,
            HashSet<String> existingRouteKeys,
            int maxColumns,
            RoutePricingConstraints constraints
    ) {
        if (activeGlobalCustomers.length == 0) {
            return new Result(false, 0.0, null, new ArrayList<RouteColumn>());
        }
        if (activeGlobalCustomers.length > 62) {
            return new Result(false, 0.0, null, new ArrayList<RouteColumn>());
        }
        BidirectionalNgRelaxedLongSearch search = new BidirectionalNgRelaxedLongSearch(
                ins,
                activeGlobalCustomers,
                qLocal,
                dualUiLocal,
                dualU0,
                existingRouteKeys,
                maxColumns,
                constraints
        );
        search.run();
        ArrayList<RouteColumn> routes = search.getCollectedRoutes();
        if (!search.aborted && search.bestRoute != null && search.bestReducedCost < -RC_EPS) {
            return new Result(true, search.bestReducedCost, search.bestRoute, routes);
        }
        return new Result(false, search.bestReducedCost, null, routes);
    }

    public Result findNegativeRoutesRelaxedThenExact(
            Instance ins,
            int[] activeGlobalCustomers,
            double[] qLocal,
            double[] dualUiLocal,
            double dualU0,
            HashSet<String> existingRouteKeys,
            int maxColumns,
            RoutePricingConstraints constraints
    ) {
        Result relaxed = findNegativeRoutesBidirectionalRelaxedOnly(
                ins,
                activeGlobalCustomers,
                qLocal,
                dualUiLocal,
                dualU0,
                existingRouteKeys,
                maxColumns,
                constraints
        );
        if (relaxed.foundNegativeColumn) {
            return relaxed;
        }
        return findNegativeRoutes(
                ins,
                activeGlobalCustomers,
                qLocal,
                dualUiLocal,
                dualU0,
                existingRouteKeys,
                maxColumns,
                constraints
        );
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
        final RoutePricingConstraints constraints;
        final CompletionBound completionBound;

        double bestReducedCost = Double.POSITIVE_INFINITY;
        RouteColumn bestRoute = null;

        Search(
                Instance ins,
                int[] customers,
                double[] q,
                double[] dual,
                double dualU0,
                HashSet<String> existingKeys,
                int maxColumns,
                RoutePricingConstraints constraints
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
            this.constraints = constraints;
            this.completionBound = CompletionBound.tryBuild(ins, customers, q, dual, constraints);
        }

        abstract void run();
        abstract boolean isStartAllowed(int startLocal);
        abstract boolean canExtend(AbstractLabel label, int nextLocal);
        abstract boolean isCompleteRouteAllowed(AbstractLabel label);

        final boolean isStartArcAllowed(int startLocal) {
            return constraints == null
                    || constraints.forbiddenStartLocal == null
                    || !constraints.forbiddenStartLocal[startLocal];
        }

        final boolean isInternalArcAllowed(int fromLocal, int toLocal) {
            return constraints == null
                    || constraints.forbiddenArcLocal == null
                    || !constraints.forbiddenArcLocal[fromLocal][toLocal];
        }

        final boolean isReturnArcAllowed(int lastLocal) {
            return constraints == null
                    || constraints.forbiddenReturnLocal == null
                    || !constraints.forbiddenReturnLocal[lastLocal];
        }

        final double transitionArcDual(int fromLocal, int toLocal) {
            if (constraints == null || constraints.arcDualLocal == null) {
                return 0.0;
            }
            return constraints.arcDualLocal[fromLocal][toLocal];
        }

        final double returnArcDual(int lastLocal) {
            if (constraints == null || constraints.returnArcDualLocal == null) {
                return 0.0;
            }
            return constraints.returnArcDualLocal[lastLocal];
        }

        final boolean hasArcDuals() {
            return constraints != null && constraints.hasArcDuals;
        }

        final ArrayList<RouteColumn> getCollectedRoutes() {
            ArrayList<RouteColumn> out = new ArrayList<RouteColumn>(collected.size());
            for (int i = 0; i < collected.size(); i++) {
                out.add(collected.get(i).route);
            }
            return out;
        }

        final boolean canStillReachUsefulRoute(AbstractLabel label) {
            if (completionBound == null) {
                return true;
            }
            double completionLb = completionBound.lowerBound(label.lastLocal, label.load);
            if (Double.isNaN(completionLb)) {
                return true;
            }
            if (!isFinite(completionLb)) {
                return false;
            }
            double routeLowerBound = label.partialReducedCost + completionLb - dualU0;
            return routeLowerBound < lowerBoundThreshold() - DOM_EPS;
        }

        final boolean canStillReachUsefulBackward(int firstLocal, double usedLoad, double partialReducedCost) {
            if (completionBound == null) {
                return true;
            }
            double completionLb = completionBound.backwardLowerBound(firstLocal, usedLoad);
            if (Double.isNaN(completionLb)) {
                return true;
            }
            if (!isFinite(completionLb)) {
                return false;
            }
            double routeLowerBound = partialReducedCost + completionLb - dualU0;
            return routeLowerBound < lowerBoundThreshold() - DOM_EPS;
        }

        private double lowerBoundThreshold() {
            if (maxColumns > 0 && collected.size() >= maxColumns) {
                return collected.get(collected.size() - 1).reducedCost;
            }
            if (bestReducedCost < -RC_EPS) {
                return bestReducedCost;
            }
            return -RC_EPS;
        }

        final void considerCompleteRoute(AbstractLabel label) {
            if (!isCompleteRouteAllowed(label)) {
                return;
            }
            int lastGlobalNode = customers[label.lastLocal];
            double fullCost = label.travelCost + ins.c[lastGlobalNode][ins.n + 1];
            double reducedCost = fullCost - label.dualGain - label.arcDualGain - dualU0 - returnArcDual(label.lastLocal);
            int[] globals = reconstructGlobals(label);
            considerCandidateRoute(globals, fullCost, label.load, reducedCost);
        }

        final void considerCandidateRoute(int[] globals, double fullCost, double load, double reducedCost) {
            boolean canImproveBest = reducedCost < bestReducedCost - DOM_EPS;
            boolean canCollect = maxColumns > 0 && reducedCost < -RC_EPS
                    && (collected.size() < maxColumns
                    || reducedCost < collected.get(collected.size() - 1).reducedCost - DOM_EPS);
            if (!canImproveBest && !canCollect) {
                return;
            }
            RouteColumn route = new RouteColumn(globals, fullCost, load);
            String key = route.key();
            if (existingKeys.contains(key)) {
                return;
            }
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

    /**
     * Capacity-state completion bound in the same role as LRP's {@code CapTwoCycleFree}:
     * from a current last customer and used load, estimate the best possible reduced-cost
     * suffix back to the depot under the node's forbidden arcs and cut duals.
     *
     * The bound is only activated when vehicle capacity and customer demands are integral.
     * Otherwise pricing falls back to its existing exact/relaxed logic without this pruning.
     */
    private static final class CompletionBound {
        private static final double INF = 1e30;
        private static final double INTEGRAL_EPS = 1e-9;

        final int capacityUnits;
        final double[][] suffixLowerBoundByUsedLoad;
        final double[][] prefixLowerBoundByUsedLoad;

        private CompletionBound(
                int capacityUnits,
                double[][] suffixLowerBoundByUsedLoad,
                double[][] prefixLowerBoundByUsedLoad
        ) {
            this.capacityUnits = capacityUnits;
            this.suffixLowerBoundByUsedLoad = suffixLowerBoundByUsedLoad;
            this.prefixLowerBoundByUsedLoad = prefixLowerBoundByUsedLoad;
        }

        static CompletionBound tryBuild(
                Instance ins,
                int[] customers,
                double[] q,
                double[] dual,
                RoutePricingConstraints constraints
        ) {
            Integer qCapacity = toIntegralUnits(ins.Q);
            if (qCapacity == null || qCapacity < 0) {
                return null;
            }
            int[] qUnits = new int[q.length];
            for (int local = 0; local < q.length; local++) {
                Integer units = toIntegralUnits(q[local]);
                if (units == null || units <= 0 || units > qCapacity) {
                    return null;
                }
                qUnits[local] = units;
            }

            int m = customers.length;
            double[][] suffixBound = new double[m][qCapacity + 1];
            double[][] prefixBound = new double[m][qCapacity + 1];
            for (int last = 0; last < m; last++) {
                for (int used = 0; used <= qCapacity; used++) {
                    suffixBound[last][used] = INF;
                    prefixBound[last][used] = INF;
                }
            }

            for (int used = qCapacity; used >= 0; used--) {
                for (int last = 0; last < m; last++) {
                    double best = returnCost(ins, customers, last, constraints);
                    for (int next = 0; next < m; next++) {
                        if (next == last) {
                            continue;
                        }
                        if (constraints != null
                                && constraints.forbiddenArcLocal != null
                                && constraints.forbiddenArcLocal[last][next]) {
                            continue;
                        }
                        int nextUsed = used + qUnits[next];
                        if (nextUsed > qCapacity) {
                            continue;
                        }
                        double tail = suffixBound[next][nextUsed];
                        if (tail >= INF / 2.0) {
                            continue;
                        }
                        double cand = ins.c[customers[last]][customers[next]]
                                - dual[next]
                                - internalArcDual(last, next, constraints)
                                + tail;
                        if (cand < best) {
                            best = cand;
                        }
                    }
                    suffixBound[last][used] = best;
                }
            }

            for (int used = qCapacity; used >= 0; used--) {
                for (int first = 0; first < m; first++) {
                    double best = startCost(ins, customers, first, constraints);
                    for (int prev = 0; prev < m; prev++) {
                        if (prev == first) {
                            continue;
                        }
                        if (constraints != null
                                && constraints.forbiddenArcLocal != null
                                && constraints.forbiddenArcLocal[prev][first]) {
                            continue;
                        }
                        int nextUsed = used + qUnits[prev];
                        if (nextUsed > qCapacity) {
                            continue;
                        }
                        double head = prefixBound[prev][nextUsed];
                        if (head >= INF / 2.0) {
                            continue;
                        }
                        double cand = head
                                + ins.c[customers[prev]][customers[first]]
                                - dual[prev]
                                - internalArcDual(prev, first, constraints);
                        if (cand < best) {
                            best = cand;
                        }
                    }
                    prefixBound[first][used] = best;
                }
            }

            return new CompletionBound(qCapacity, suffixBound, prefixBound);
        }

        double lowerBound(int lastLocal, double usedLoad) {
            Integer usedUnits = toIntegralUnits(usedLoad);
            if (usedUnits == null || usedUnits < 0 || usedUnits > capacityUnits) {
                return Double.NaN;
            }
            return suffixLowerBoundByUsedLoad[lastLocal][usedUnits];
        }

        double backwardLowerBound(int firstLocal, double usedLoad) {
            Integer usedUnits = toIntegralUnits(usedLoad);
            if (usedUnits == null || usedUnits < 0 || usedUnits > capacityUnits) {
                return Double.NaN;
            }
            return prefixLowerBoundByUsedLoad[firstLocal][usedUnits];
        }

        private static double returnCost(
                Instance ins,
                int[] customers,
                int lastLocal,
                RoutePricingConstraints constraints
        ) {
            if (constraints != null
                    && constraints.forbiddenReturnLocal != null
                    && constraints.forbiddenReturnLocal[lastLocal]) {
                return INF;
            }
            double dual = 0.0;
            if (constraints != null && constraints.returnArcDualLocal != null) {
                dual = constraints.returnArcDualLocal[lastLocal];
            }
            return ins.c[customers[lastLocal]][ins.n + 1] - dual;
        }

        private static double startCost(
                Instance ins,
                int[] customers,
                int firstLocal,
                RoutePricingConstraints constraints
        ) {
            if (constraints != null
                    && constraints.forbiddenStartLocal != null
                    && constraints.forbiddenStartLocal[firstLocal]) {
                return INF;
            }
            return ins.c[0][customers[firstLocal]];
        }

        private static double internalArcDual(int fromLocal, int toLocal, RoutePricingConstraints constraints) {
            if (constraints == null || constraints.arcDualLocal == null) {
                return 0.0;
            }
            return constraints.arcDualLocal[fromLocal][toLocal];
        }

        private static Integer toIntegralUnits(double value) {
            long rounded = Math.round(value);
            if (Math.abs(value - rounded) > INTEGRAL_EPS) {
                return null;
            }
            if (rounded < Integer.MIN_VALUE || rounded > Integer.MAX_VALUE) {
                return null;
            }
            return (int) rounded;
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
        final double arcDualGain; // sum of cut duals on traversed arcs excluding return arc
        final double partialReducedCost; // travelCost - dualGain
        final AbstractLabel pred;

        AbstractLabel(
                int lastLocal,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                double arcDualGain,
                AbstractLabel pred
        ) {
            this.lastLocal = lastLocal;
            this.depth = depth;
            this.load = load;
            this.travelCost = travelCost;
            this.dualGain = dualGain;
            this.arcDualGain = arcDualGain;
            this.partialReducedCost = travelCost - dualGain - arcDualGain;
            this.pred = pred;
        }
    }

    private static final class BidirectionalLongSearch extends Search {
        private final int customerCount;
        private final long fullMask;
        private final long[] localBit;
        private final ArrayList<ForwardLabel>[] forwardUnextended;
        private final ArrayList<ForwardLabel>[] forwardExtended;
        private final ArrayList<BackwardLabel>[] backwardUnextended;
        private final ArrayList<BackwardLabel>[] backwardExtended;
        private final HashMap<Long, ForwardLabel>[] bestForwardByMaskAtLast;
        private final HashMap<Long, BackwardLabel>[] bestBackwardByMaskAtFirst;
        private final int[] joinedForwardCountByLast;
        private final int[] joinedBackwardCountByFirst;
        private final double stepSize;

        BidirectionalLongSearch(
                Instance ins,
                int[] customers,
                double[] q,
                double[] dual,
                double dualU0,
                HashSet<String> existingKeys,
                int maxColumns,
                RoutePricingConstraints constraints
        ) {
            super(ins, customers, q, dual, dualU0, existingKeys, maxColumns, constraints);
            this.customerCount = customers.length;
            this.fullMask = (1L << customers.length) - 1L;
            this.localBit = new long[customers.length];
            for (int local = 0; local < customers.length; local++) {
                this.localBit[local] = 1L << local;
            }
            @SuppressWarnings("unchecked")
            ArrayList<ForwardLabel>[] fUn = (ArrayList<ForwardLabel>[]) new ArrayList[customerCount];
            @SuppressWarnings("unchecked")
            ArrayList<ForwardLabel>[] fEx = (ArrayList<ForwardLabel>[]) new ArrayList[customerCount];
            @SuppressWarnings("unchecked")
            ArrayList<BackwardLabel>[] bUn = (ArrayList<BackwardLabel>[]) new ArrayList[customerCount];
            @SuppressWarnings("unchecked")
            ArrayList<BackwardLabel>[] bEx = (ArrayList<BackwardLabel>[]) new ArrayList[customerCount];
            @SuppressWarnings("unchecked")
            HashMap<Long, ForwardLabel>[] fBest = (HashMap<Long, ForwardLabel>[]) new HashMap[customerCount];
            @SuppressWarnings("unchecked")
            HashMap<Long, BackwardLabel>[] bBest = (HashMap<Long, BackwardLabel>[]) new HashMap[customerCount];
            for (int local = 0; local < customerCount; local++) {
                fUn[local] = new ArrayList<ForwardLabel>();
                fEx[local] = new ArrayList<ForwardLabel>();
                bUn[local] = new ArrayList<BackwardLabel>();
                bEx[local] = new ArrayList<BackwardLabel>();
                fBest[local] = new HashMap<Long, ForwardLabel>();
                bBest[local] = new HashMap<Long, BackwardLabel>();
            }
            this.forwardUnextended = fUn;
            this.forwardExtended = fEx;
            this.backwardUnextended = bUn;
            this.backwardExtended = bEx;
            this.bestForwardByMaskAtLast = fBest;
            this.bestBackwardByMaskAtFirst = bBest;
            this.joinedForwardCountByLast = new int[customerCount];
            this.joinedBackwardCountByFirst = new int[customerCount];
            this.stepSize = Math.max(1.0, ins.Q / 15.0);
        }

        @Override
        void run() {
            initializeBoundaryLabels();
            if (customerCount == 0) {
                return;
            }

            double forwardLimit = stepSize;
            double backwardLimit = stepSize;
            while (true) {
                extendForwardUpTo(forwardLimit);
                extendBackwardUpTo(backwardLimit);
                joinExtendedLabels();

                if (forwardLimit + backwardLimit >= ins.Q - 1e-9) {
                    extendForwardUpTo(ins.Q + 1.0);
                    extendBackwardUpTo(ins.Q + 1.0);
                    joinExtendedLabels();
                    return;
                }

                int remainingForward = countRemaining(forwardUnextended);
                int remainingBackward = countRemaining(backwardUnextended);
                if (remainingForward == 0 && remainingBackward == 0) {
                    return;
                }

                if (remainingForward <= remainingBackward) {
                    forwardLimit = Math.min(ins.Q, forwardLimit + stepSize);
                } else {
                    backwardLimit = Math.min(ins.Q, backwardLimit + stepSize);
                }
            }
        }

        @Override
        boolean isStartAllowed(int startLocal) {
            return isStartArcAllowed(startLocal);
        }

        @Override
        boolean canExtend(AbstractLabel label, int nextLocal) {
            if (!isInternalArcAllowed(label.lastLocal, nextLocal)) {
                return false;
            }
            if ((localBit[nextLocal] & ((ForwardLabel) label).visitedMask) != 0L) {
                return false;
            }
            return !violatesForbiddenPairMask(((ForwardLabel) label).visitedMask, nextLocal);
        }

        @Override
        boolean isCompleteRouteAllowed(AbstractLabel label) {
            return isReturnArcAllowed(label.lastLocal)
                    && isRouteMaskAllowed(((ForwardLabel) label).visitedMask);
        }

        private void initializeBoundaryLabels() {
            for (int local = 0; local < customerCount; local++) {
                if (q[local] > ins.Q + 1e-9) {
                    continue;
                }
                long mask = localBit[local];
                int global = customers[local];

                if (isStartAllowed(local)) {
                    ForwardLabel forward = new ForwardLabel(
                            local,
                            1,
                            q[local],
                            ins.c[0][global],
                            dual[local],
                            0.0,
                            null,
                            mask
                    );
                    if (registerForward(forward)) {
                        forwardUnextended[local].add(forward);
                    }
                }

                if (isReturnArcAllowed(local)) {
                    BackwardLabel backward = new BackwardLabel(
                            local,
                            1,
                            q[local],
                            ins.c[global][ins.n + 1],
                            dual[local],
                            returnArcDual(local),
                            null,
                            mask
                    );
                    if (registerBackward(backward)) {
                        backwardUnextended[local].add(backward);
                    }
                }
            }
        }

        private void extendForwardUpTo(double loadLimit) {
            while (true) {
                ForwardLabel label = pollForward(loadLimit);
                if (label == null) {
                    return;
                }
                if (!isCurrentForward(label)) {
                    continue;
                }
                forwardExtended[label.lastLocal].add(label);
                considerCompleteRoute(label);
                if (!canStillReachUsefulRoute(label)) {
                    continue;
                }

                long remaining = fullMask & ~label.visitedMask;
                while (remaining != 0L) {
                    long bit = remaining & -remaining;
                    int next = Long.numberOfTrailingZeros(bit);
                    remaining ^= bit;

                    if (!canExtend(label, next)) {
                        continue;
                    }
                    double newLoad = label.load + q[next];
                    if (newLoad > ins.Q + 1e-9) {
                        continue;
                    }
                    int fromGlobal = customers[label.lastLocal];
                    int toGlobal = customers[next];
                    ForwardLabel child = new ForwardLabel(
                            next,
                            label.depth + 1,
                            newLoad,
                            label.travelCost + ins.c[fromGlobal][toGlobal],
                            label.dualGain + dual[next],
                            label.arcDualGain + transitionArcDual(label.lastLocal, next),
                            label,
                            label.visitedMask | bit
                    );
                    if (registerForward(child)) {
                        forwardUnextended[next].add(child);
                    }
                }
            }
        }

        private void extendBackwardUpTo(double loadLimit) {
            while (true) {
                BackwardLabel label = pollBackward(loadLimit);
                if (label == null) {
                    return;
                }
                if (!isCurrentBackward(label)) {
                    continue;
                }
                backwardExtended[label.firstLocal].add(label);
                considerCompleteBackwardRoute(label);
                if (!canStillReachUsefulBackward(label.firstLocal, label.load, label.partialReducedCost)) {
                    continue;
                }

                long remaining = fullMask & ~label.visitedMask;
                while (remaining != 0L) {
                    long bit = remaining & -remaining;
                    int prev = Long.numberOfTrailingZeros(bit);
                    remaining ^= bit;

                    if (!canPrependBackward(label, prev)) {
                        continue;
                    }
                    double newLoad = label.load + q[prev];
                    if (newLoad > ins.Q + 1e-9) {
                        continue;
                    }
                    int prevGlobal = customers[prev];
                    int firstGlobal = customers[label.firstLocal];
                    BackwardLabel child = new BackwardLabel(
                            prev,
                            label.depth + 1,
                            newLoad,
                            ins.c[prevGlobal][firstGlobal] + label.travelCost,
                            label.dualGain + dual[prev],
                            label.arcDualGain + transitionArcDual(prev, label.firstLocal),
                            label,
                            label.visitedMask | bit
                    );
                    if (registerBackward(child)) {
                        backwardUnextended[prev].add(child);
                    }
                }
            }
        }

        private void joinExtendedLabels() {
            for (int forwardLast = 0; forwardLast < customerCount; forwardLast++) {
                ArrayList<ForwardLabel> fLabels = forwardExtended[forwardLast];
                if (fLabels.isEmpty()) {
                    continue;
                }
                int newForwardStart = joinedForwardCountByLast[forwardLast];
                for (int backwardFirst = 0; backwardFirst < customerCount; backwardFirst++) {
                    if (!isInternalArcAllowed(forwardLast, backwardFirst)) {
                        continue;
                    }
                    ArrayList<BackwardLabel> bLabels = backwardExtended[backwardFirst];
                    if (bLabels.isEmpty()) {
                        continue;
                    }
                    int newBackwardStart = joinedBackwardCountByFirst[backwardFirst];
                    if (newForwardStart < fLabels.size()) {
                        joinLabelRanges(fLabels, newForwardStart, fLabels.size(), bLabels, 0, bLabels.size());
                    }
                    if (newBackwardStart < bLabels.size() && newForwardStart > 0) {
                        joinLabelRanges(fLabels, 0, newForwardStart, bLabels, newBackwardStart, bLabels.size());
                    }
                }
            }
            for (int local = 0; local < customerCount; local++) {
                joinedForwardCountByLast[local] = forwardExtended[local].size();
                joinedBackwardCountByFirst[local] = backwardExtended[local].size();
            }
        }

        private void joinLabelRanges(
                ArrayList<ForwardLabel> fLabels,
                int fFrom,
                int fTo,
                ArrayList<BackwardLabel> bLabels,
                int bFrom,
                int bTo
        ) {
            for (int fi = fFrom; fi < fTo; fi++) {
                ForwardLabel fLabel = fLabels.get(fi);
                for (int bi = bFrom; bi < bTo; bi++) {
                    BackwardLabel bLabel = bLabels.get(bi);
                    if ((fLabel.visitedMask & bLabel.visitedMask) != 0L) {
                        continue;
                    }
                    double combinedLoad = fLabel.load + bLabel.load;
                    if (combinedLoad > ins.Q + 1e-9) {
                        continue;
                    }
                    long routeMask = fLabel.visitedMask | bLabel.visitedMask;
                    if (!isRouteMaskAllowed(routeMask)) {
                        continue;
                    }
                    double fullCost = fLabel.travelCost
                            + ins.c[customers[fLabel.lastLocal]][customers[bLabel.firstLocal]]
                            + bLabel.travelCost;
                    double reducedCost = fullCost
                            - fLabel.dualGain
                            - bLabel.dualGain
                            - fLabel.arcDualGain
                            - bLabel.arcDualGain
                            - transitionArcDual(fLabel.lastLocal, bLabel.firstLocal)
                            - dualU0;
                    int[] globals = combine(reconstructForwardGlobals(fLabel), reconstructBackwardGlobals(bLabel));
                    considerCandidateRoute(globals, fullCost, combinedLoad, reducedCost);
                }
            }
        }

        private void considerCompleteBackwardRoute(BackwardLabel label) {
            if (!isStartArcAllowed(label.firstLocal)) {
                return;
            }
            if (!isRouteMaskAllowed(label.visitedMask)) {
                return;
            }
            double fullCost = ins.c[0][customers[label.firstLocal]] + label.travelCost;
            double reducedCost = fullCost - label.dualGain - label.arcDualGain - dualU0;
            considerCandidateRoute(reconstructBackwardGlobals(label), fullCost, label.load, reducedCost);
        }

        private boolean canPrependBackward(BackwardLabel label, int prevLocal) {
            if (!isInternalArcAllowed(prevLocal, label.firstLocal)) {
                return false;
            }
            if ((label.visitedMask & localBit[prevLocal]) != 0L) {
                return false;
            }
            return !violatesForbiddenPairMask(label.visitedMask, prevLocal);
        }

        private boolean isRouteMaskAllowed(long visitedMask) {
            if (constraints == null) {
                return true;
            }
            if (constraints.forbiddenWithLocalLong != null) {
                long bits = visitedMask;
                while (bits != 0L) {
                    int local = Long.numberOfTrailingZeros(bits);
                    if ((visitedMask & constraints.forbiddenWithLocalLong[local]) != 0L) {
                        return false;
                    }
                    bits ^= (bits & -bits);
                }
            } else if (constraints.forbiddenWithLocalBits != null) {
                BitSet visited = bitSetFromMask(visitedMask);
                int bit = visited.nextSetBit(0);
                while (bit >= 0) {
                    BitSet forbidden = constraints.forbiddenWithLocalBits[bit];
                    if (forbidden != null) {
                        BitSet overlap = (BitSet) visited.clone();
                        overlap.and(forbidden);
                        if (!overlap.isEmpty()) {
                            return false;
                        }
                    }
                    bit = visited.nextSetBit(bit + 1);
                }
            }

            if (constraints.requiredWithLocalLong != null) {
                long bits = visitedMask;
                while (bits != 0L) {
                    int local = Long.numberOfTrailingZeros(bits);
                    if ((visitedMask & constraints.requiredWithLocalLong[local]) != constraints.requiredWithLocalLong[local]) {
                        return false;
                    }
                    bits ^= (bits & -bits);
                }
            } else if (constraints.requiredWithLocalBits != null) {
                BitSet visited = bitSetFromMask(visitedMask);
                int bit = visited.nextSetBit(0);
                while (bit >= 0) {
                    BitSet required = constraints.requiredWithLocalBits[bit];
                    if (required != null) {
                        BitSet missing = (BitSet) required.clone();
                        missing.andNot(visited);
                        if (!missing.isEmpty()) {
                            return false;
                        }
                    }
                    bit = visited.nextSetBit(bit + 1);
                }
            }
            return true;
        }

        private boolean violatesForbiddenPairMask(long visitedMask, int candidateLocal) {
            if (constraints == null) {
                return false;
            }
            if (constraints.forbiddenWithLocalLong != null) {
                return (visitedMask & constraints.forbiddenWithLocalLong[candidateLocal]) != 0L;
            }
            if (constraints.forbiddenWithLocalBits != null) {
                BitSet forbidden = constraints.forbiddenWithLocalBits[candidateLocal];
                if (forbidden == null || forbidden.isEmpty()) {
                    return false;
                }
                BitSet visited = bitSetFromMask(visitedMask);
                visited.and(forbidden);
                return !visited.isEmpty();
            }
            return false;
        }

        private BitSet bitSetFromMask(long mask) {
            BitSet bits = new BitSet(customerCount);
            long remaining = mask;
            while (remaining != 0L) {
                long bit = remaining & -remaining;
                bits.set(Long.numberOfTrailingZeros(bit));
                remaining ^= bit;
            }
            return bits;
        }

        private ForwardLabel pollForward(double loadLimit) {
            for (int local = 0; local < customerCount; local++) {
                ArrayList<ForwardLabel> labels = forwardUnextended[local];
                for (int idx = 0; idx < labels.size(); idx++) {
                    ForwardLabel label = labels.get(idx);
                    if (label.load <= loadLimit + 1e-9) {
                        labels.remove(idx);
                        return label;
                    }
                }
            }
            return null;
        }

        private BackwardLabel pollBackward(double loadLimit) {
            for (int local = 0; local < customerCount; local++) {
                ArrayList<BackwardLabel> labels = backwardUnextended[local];
                for (int idx = 0; idx < labels.size(); idx++) {
                    BackwardLabel label = labels.get(idx);
                    if (label.load <= loadLimit + 1e-9) {
                        labels.remove(idx);
                        return label;
                    }
                }
            }
            return null;
        }

        private int countRemaining(ArrayList<?>[] lists) {
            int total = 0;
            for (int local = 0; local < lists.length; local++) {
                total += lists[local].size();
            }
            return total;
        }

        private boolean registerForward(ForwardLabel label) {
            HashMap<Long, ForwardLabel> exactMap = bestForwardByMaskAtLast[label.lastLocal];
            ForwardLabel incumbent = exactMap.get(label.visitedMask);
            if (incumbent != null && incumbent.partialReducedCost <= label.partialReducedCost + DOM_EPS) {
                return false;
            }
            if (isForwardDominatedBySubset(label, exactMap)) {
                return false;
            }
            exactMap.put(label.visitedMask, label);
            removeDominatedForwardSupersets(label, exactMap);
            return true;
        }

        private boolean registerBackward(BackwardLabel label) {
            HashMap<Long, BackwardLabel> exactMap = bestBackwardByMaskAtFirst[label.firstLocal];
            BackwardLabel incumbent = exactMap.get(label.visitedMask);
            if (incumbent != null && incumbent.partialReducedCost <= label.partialReducedCost + DOM_EPS) {
                return false;
            }
            if (isBackwardDominatedBySubset(label, exactMap)) {
                return false;
            }
            exactMap.put(label.visitedMask, label);
            removeDominatedBackwardSupersets(label, exactMap);
            return true;
        }

        private boolean isCurrentForward(ForwardLabel label) {
            return bestForwardByMaskAtLast[label.lastLocal].get(label.visitedMask) == label;
        }

        private boolean isCurrentBackward(BackwardLabel label) {
            return bestBackwardByMaskAtFirst[label.firstLocal].get(label.visitedMask) == label;
        }

        private boolean isForwardDominatedBySubset(ForwardLabel label, HashMap<Long, ForwardLabel> exactMap) {
            if (label.depth <= 1) {
                return false;
            }
            long lastBit = localBit[label.lastLocal];
            long free = label.visitedMask ^ lastBit;
            if (!shouldEnumerateSubsetDominance(free)) {
                return false;
            }
            long subFree = free;
            while (subFree != 0L) {
                long subsetMask = subFree | lastBit;
                ForwardLabel subset = exactMap.get(subsetMask);
                if (subset != null && subset.partialReducedCost <= label.partialReducedCost + DOM_EPS) {
                    return true;
                }
                subFree = (subFree - 1L) & free;
            }
            return false;
        }

        private boolean isBackwardDominatedBySubset(BackwardLabel label, HashMap<Long, BackwardLabel> exactMap) {
            if (label.depth <= 1) {
                return false;
            }
            long firstBit = localBit[label.firstLocal];
            long free = label.visitedMask ^ firstBit;
            if (!shouldEnumerateSubsetDominance(free)) {
                return false;
            }
            long subFree = free;
            while (subFree != 0L) {
                long subsetMask = subFree | firstBit;
                BackwardLabel subset = exactMap.get(subsetMask);
                if (subset != null && subset.partialReducedCost <= label.partialReducedCost + DOM_EPS) {
                    return true;
                }
                subFree = (subFree - 1L) & free;
            }
            return false;
        }

        private void removeDominatedForwardSupersets(ForwardLabel label, HashMap<Long, ForwardLabel> exactMap) {
            removeDominatedSupersets(label.visitedMask, label.lastLocal, label.partialReducedCost, exactMap);
        }

        private void removeDominatedBackwardSupersets(BackwardLabel label, HashMap<Long, BackwardLabel> exactMap) {
            removeDominatedSupersets(label.visitedMask, label.firstLocal, label.partialReducedCost, exactMap);
        }

        private <T extends PartialRcState> void removeDominatedSupersets(
                long mask,
                int fixedLocal,
                double partialReducedCost,
                HashMap<Long, T> exactMap
        ) {
            long fixedBit = localBit[fixedLocal];
            long baseMask = mask | fixedBit;
            long extraPool = fullMask & ~baseMask;
            int extraCount = Long.bitCount(extraPool);
            if (extraCount == 0 || extraCount > SUPERSET_CLEANUP_MAX_EXTRA_BITS) {
                return;
            }
            long add = extraPool;
            while (add != 0L) {
                T other = exactMap.get(baseMask | add);
                if (other != null && partialReducedCost <= other.partialReducedCost() + DOM_EPS) {
                    exactMap.remove(baseMask | add);
                }
                add = (add - 1L) & extraPool;
            }
        }

        private int[] reconstructForwardGlobals(ForwardLabel label) {
            int[] globals = new int[label.depth];
            ForwardLabel cur = label;
            for (int pos = label.depth - 1; pos >= 0; pos--) {
                globals[pos] = customers[cur.lastLocal];
                cur = cur.pred;
            }
            return globals;
        }

        private int[] reconstructBackwardGlobals(BackwardLabel label) {
            int[] globals = new int[label.depth];
            BackwardLabel cur = label;
            for (int pos = 0; pos < label.depth; pos++) {
                globals[pos] = customers[cur.firstLocal];
                cur = cur.pred;
            }
            return globals;
        }

        private int[] combine(int[] left, int[] right) {
            int[] out = new int[left.length + right.length];
            System.arraycopy(left, 0, out, 0, left.length);
            System.arraycopy(right, 0, out, left.length, right.length);
            return out;
        }
    }

    private interface PartialRcState {
        double partialReducedCost();
    }

    private static final class ForwardLabel extends AbstractLabel implements PartialRcState {
        final long visitedMask;
        final ForwardLabel pred;

        ForwardLabel(
                int lastLocal,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                double arcDualGain,
                ForwardLabel pred,
                long visitedMask
        ) {
            super(lastLocal, depth, load, travelCost, dualGain, arcDualGain, pred);
            this.pred = pred;
            this.visitedMask = visitedMask;
        }

        @Override
        public double partialReducedCost() {
            return partialReducedCost;
        }
    }

    private static final class BackwardLabel implements PartialRcState {
        final int firstLocal;
        final int depth;
        final double load;
        final double travelCost; // first -> ... -> return depot
        final double dualGain;
        final double arcDualGain;
        final double partialReducedCost;
        final BackwardLabel pred; // next customer toward the return depot
        final long visitedMask;

        BackwardLabel(
                int firstLocal,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                double arcDualGain,
                BackwardLabel pred,
                long visitedMask
        ) {
            this.firstLocal = firstLocal;
            this.depth = depth;
            this.load = load;
            this.travelCost = travelCost;
            this.dualGain = dualGain;
            this.arcDualGain = arcDualGain;
            this.partialReducedCost = travelCost - dualGain - arcDualGain;
            this.pred = pred;
            this.visitedMask = visitedMask;
        }

        @Override
        public double partialReducedCost() {
            return partialReducedCost;
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
                int maxColumns,
                RoutePricingConstraints constraints
        ) {
            super(ins, customers, q, dual, dualU0, existingKeys, maxColumns, constraints);
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
                if (q[start] > ins.Q + 1e-9 || !isStartAllowed(start)) {
                    continue;
                }
                long mask = (1L << start);
                int g = customers[start];
                LongLabel startLabel = new LongLabel(
                        start, 1, q[start], ins.c[0][g], dual[start], 0.0, null, mask
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
                if (!canStillReachUsefulRoute(cur)) {
                    continue;
                }
                if (cannotBeatBest(cur)) {
                    continue;
                }

                long remaining = fullMask & ~cur.visitedMask;
                while (remaining != 0L) {
                    long bit = remaining & -remaining;
                    int nxt = Long.numberOfTrailingZeros(bit);
                    remaining ^= bit;

                    if (!canExtend(cur, nxt)) {
                        continue;
                    }
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
                            cur.arcDualGain + transitionArcDual(cur.lastLocal, nxt),
                            cur,
                            cur.visitedMask | bit
                    );
                    if (register(child)) {
                        queue.addLast(child);
                    }
                }
            }
        }

        @Override
        boolean isStartAllowed(int startLocal) {
            return isStartArcAllowed(startLocal);
        }

        @Override
        boolean canExtend(AbstractLabel label, int nextLocal) {
            if (!isInternalArcAllowed(label.lastLocal, nextLocal)) {
                return false;
            }
            if (constraints == null || constraints.forbiddenWithLocalLong == null) {
                return true;
            }
            LongLabel longLabel = (LongLabel) label;
            return (longLabel.visitedMask & constraints.forbiddenWithLocalLong[nextLocal]) == 0L;
        }

        @Override
        boolean isCompleteRouteAllowed(AbstractLabel label) {
            if (!isReturnArcAllowed(label.lastLocal)) {
                return false;
            }
            if (constraints == null || constraints.requiredWithLocalLong == null) {
                return true;
            }
            long visited = ((LongLabel) label).visitedMask;
            long bits = visited;
            while (bits != 0L) {
                int local = Long.numberOfTrailingZeros(bits);
                if ((visited & constraints.requiredWithLocalLong[local]) != constraints.requiredWithLocalLong[local]) {
                    return false;
                }
                bits ^= (bits & -bits);
            }
            return true;
        }

        private boolean register(LongLabel label) {
            HashMap<Long, LongLabel> exactMap = bestByMaskAtLast[label.lastLocal];
            LongLabel incumbent = exactMap.get(label.visitedMask);
            if (incumbent != null) {
                if (incumbent.partialReducedCost <= label.partialReducedCost + DOM_EPS) {
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
            if (!shouldEnumerateSubsetDominance(free)) {
                return false;
            }
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
            if (hasArcDuals()) {
                return false;
            }
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
                double arcDualGain,
                AbstractLabel pred,
                long visitedMask
        ) {
            super(lastLocal, depth, load, travelCost, dualGain, arcDualGain, pred);
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
                int maxColumns,
                RoutePricingConstraints constraints
        ) {
            super(ins, customers, q, dual, dualU0, existingKeys, maxColumns, constraints);
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
                if (q[start] > ins.Q + 1e-9 || !isStartAllowed(start)) {
                    continue;
                }
                BitSet visited = new BitSet(customerCount);
                visited.set(start);
                int g = customers[start];
                BitSetLabel startLabel = new BitSetLabel(
                        start, 1, q[start], ins.c[0][g], dual[start], 0.0, null, visited
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
                if (!canStillReachUsefulRoute(cur)) {
                    continue;
                }

                int nxt = cur.visited.nextClearBit(0);
                while (nxt >= 0 && nxt < customerCount) {
                    if (!canExtend(cur, nxt)) {
                        nxt = cur.visited.nextClearBit(nxt + 1);
                        continue;
                    }
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
                                cur.arcDualGain + transitionArcDual(cur.lastLocal, nxt),
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

        @Override
        boolean isStartAllowed(int startLocal) {
            return isStartArcAllowed(startLocal);
        }

        @Override
        boolean canExtend(AbstractLabel label, int nextLocal) {
            if (!isInternalArcAllowed(label.lastLocal, nextLocal)) {
                return false;
            }
            if (constraints == null || constraints.forbiddenWithLocalBits == null) {
                return true;
            }
            BitSetLabel bitSetLabel = (BitSetLabel) label;
            BitSet forbidden = constraints.forbiddenWithLocalBits[nextLocal];
            if (forbidden == null || forbidden.isEmpty()) {
                return true;
            }
            BitSet overlap = (BitSet) bitSetLabel.visited.clone();
            overlap.and(forbidden);
            return overlap.isEmpty();
        }

        @Override
        boolean isCompleteRouteAllowed(AbstractLabel label) {
            if (!isReturnArcAllowed(label.lastLocal)) {
                return false;
            }
            if (constraints == null || constraints.requiredWithLocalBits == null) {
                return true;
            }
            BitSet visited = ((BitSetLabel) label).visited;
            int bit = visited.nextSetBit(0);
            while (bit >= 0) {
                BitSet required = constraints.requiredWithLocalBits[bit];
                if (required != null && !required.isEmpty()) {
                    BitSet missing = (BitSet) required.clone();
                    missing.andNot(visited);
                    if (!missing.isEmpty()) {
                        return false;
                    }
                }
                bit = visited.nextSetBit(bit + 1);
            }
            return true;
        }

        private boolean register(BitSetLabel label) {
            BitSetStateKey key = new BitSetStateKey(label.lastLocal, label.visited);
            BitSetLabel incumbent = bestByState.get(key);
            if (incumbent != null) {
                if (incumbent.partialReducedCost <= label.partialReducedCost + DOM_EPS) {
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
                double arcDualGain,
                AbstractLabel pred,
                BitSet visited
        ) {
            super(lastLocal, depth, load, travelCost, dualGain, arcDualGain, pred);
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

    private static final class BidirectionalNgRelaxedLongSearch {
        private static final int DEFAULT_NG_SIZE = 5;

        final Instance ins;
        final int[] customers;
        final double[] q;
        final double[] dual;
        final double dualU0;
        final HashSet<String> existingKeys;
        final int maxColumns;
        final RoutePricingConstraints constraints;
        final long[] localBit;
        final ArrayList<CollectedColumn> collected;
        final HashSet<String> collectedKeys;
        final long[] ngMaskByLocal;
        final int maxDepth;
        final double stepSize;
        final CompletionBound completionBound;
        boolean aborted = false;

        double bestReducedCost = Double.POSITIVE_INFINITY;
        RouteColumn bestRoute = null;

        private double bestAnyReducedCost = Double.POSITIVE_INFINITY;
        private int[] bestAnyLocals = null;

        BidirectionalNgRelaxedLongSearch(
                Instance ins,
                int[] customers,
                double[] q,
                double[] dual,
                double dualU0,
                HashSet<String> existingKeys,
                int maxColumns,
                RoutePricingConstraints constraints
        ) {
            this.ins = ins;
            this.customers = customers;
            this.q = q;
            this.dual = dual;
            this.dualU0 = dualU0;
            this.existingKeys = existingKeys;
            this.maxColumns = Math.max(0, maxColumns);
            this.constraints = constraints;
            this.localBit = new long[customers.length];
            for (int local = 0; local < customers.length; local++) {
                this.localBit[local] = 1L << local;
            }
            this.collected = new ArrayList<CollectedColumn>();
            this.collectedKeys = new HashSet<String>();
            this.ngMaskByLocal = buildInitialNgMasks();
            this.maxDepth = Math.max(customers.length + 2, Math.min(customers.length * 2, customers.length + 10));
            this.stepSize = Math.max(1.0, ins.Q / 4.0);
            this.completionBound = CompletionBound.tryBuild(ins, customers, q, dual, constraints);
        }

        ArrayList<RouteColumn> getCollectedRoutes() {
            ArrayList<RouteColumn> out = new ArrayList<RouteColumn>(collected.size());
            for (int idx = 0; idx < collected.size(); idx++) {
                out.add(collected.get(idx).route);
            }
            return out;
        }

        void run() {
            long startNs = System.nanoTime();
            int rounds = 1;
            for (int round = 0; round < rounds; round++) {
                bestAnyReducedCost = Double.POSITIVE_INFINITY;
                bestAnyLocals = null;
                searchOnce(startNs);
                if (aborted || !collected.isEmpty()) {
                    return;
                }
                if (bestAnyLocals == null || bestAnyReducedCost >= -RC_EPS) {
                    return;
                }
                if (!updateNgSetsFromCycles(bestAnyLocals)) {
                    return;
                }
            }
        }

        private void searchOnce(long startNs) {
            @SuppressWarnings("unchecked")
            ArrayList<RelaxForwardLabel>[] forwardUnextended = (ArrayList<RelaxForwardLabel>[]) new ArrayList[customers.length];
            @SuppressWarnings("unchecked")
            ArrayList<RelaxForwardLabel>[] forwardExtended = (ArrayList<RelaxForwardLabel>[]) new ArrayList[customers.length];
            @SuppressWarnings("unchecked")
            ArrayList<RelaxBackwardLabel>[] backwardUnextended = (ArrayList<RelaxBackwardLabel>[]) new ArrayList[customers.length];
            @SuppressWarnings("unchecked")
            ArrayList<RelaxBackwardLabel>[] backwardExtended = (ArrayList<RelaxBackwardLabel>[]) new ArrayList[customers.length];
            for (int local = 0; local < customers.length; local++) {
                forwardUnextended[local] = new ArrayList<RelaxForwardLabel>();
                forwardExtended[local] = new ArrayList<RelaxForwardLabel>();
                backwardUnextended[local] = new ArrayList<RelaxBackwardLabel>();
                backwardExtended[local] = new ArrayList<RelaxBackwardLabel>();
            }

            initializeBoundaryRelaxedLabels(forwardUnextended, backwardUnextended);

            double forwardLimit = stepSize;
            double backwardLimit = stepSize;
            while (true) {
                if (isOverBudget(startNs)) {
                    aborted = true;
                    return;
                }
                extendForward(forwardUnextended, forwardExtended, forwardLimit, startNs);
                if (aborted) {
                    return;
                }
                extendBackward(backwardUnextended, backwardExtended, backwardLimit, startNs);
                if (aborted) {
                    return;
                }
                joinRelaxedLabels(forwardExtended, backwardExtended, startNs);
                if (!collected.isEmpty()) {
                    return;
                }
                if (forwardLimit + backwardLimit >= ins.Q - 1e-9) {
                    extendForward(forwardUnextended, forwardExtended, ins.Q + 1.0, startNs);
                    if (aborted) {
                        return;
                    }
                    extendBackward(backwardUnextended, backwardExtended, ins.Q + 1.0, startNs);
                    if (aborted) {
                        return;
                    }
                    joinRelaxedLabels(forwardExtended, backwardExtended, startNs);
                    return;
                }
                int remainingForward = countRemaining(forwardUnextended);
                int remainingBackward = countRemaining(backwardUnextended);
                if (remainingForward == 0 && remainingBackward == 0) {
                    return;
                }
                if (remainingForward <= remainingBackward) {
                    forwardLimit = Math.min(ins.Q, forwardLimit + stepSize);
                } else {
                    backwardLimit = Math.min(ins.Q, backwardLimit + stepSize);
                }
            }
        }

        private void initializeBoundaryRelaxedLabels(
                ArrayList<RelaxForwardLabel>[] forwardUnextended,
                ArrayList<RelaxBackwardLabel>[] backwardUnextended
        ) {
            for (int local = 0; local < customers.length; local++) {
                if (q[local] > ins.Q + 1e-9) {
                    continue;
                }
                long bit = localBit[local];
                double load = q[local];
                long memory = capacityForbiddenMask(load) | bit;
                int global = customers[local];

                if (constraints == null || constraints.forbiddenStartLocal == null || !constraints.forbiddenStartLocal[local]) {
                    forwardUnextended[local].add(new RelaxForwardLabel(
                            local,
                            1,
                            load,
                            ins.c[0][global],
                            dual[local],
                            0.0,
                            memory,
                            bit,
                            0L,
                            null
                    ));
                }

                if (constraints == null || constraints.forbiddenReturnLocal == null || !constraints.forbiddenReturnLocal[local]) {
                    backwardUnextended[local].add(new RelaxBackwardLabel(
                            local,
                            1,
                            load,
                            ins.c[global][ins.n + 1],
                            dual[local],
                            returnArcDual(local),
                            memory,
                            bit,
                            0L,
                            null
                    ));
                }
            }
        }

        private void extendForward(
                ArrayList<RelaxForwardLabel>[] forwardUnextended,
                ArrayList<RelaxForwardLabel>[] forwardExtended,
                double loadLimit,
                long startNs
        ) {
            while (true) {
                if (isOverBudget(startNs)) {
                    aborted = true;
                    return;
                }
                RelaxForwardLabel label = pollForward(forwardUnextended, loadLimit);
                if (label == null) {
                    return;
                }
                insertSortedByPartialReducedCost(forwardExtended[label.lastLocal], label);
                considerForwardCompleteRoute(label);
                if (maxColumns > 0 && collected.size() >= maxColumns) {
                    return;
                }
                if (!canStillReachUsefulForward(label)) {
                    continue;
                }
                if (label.depth >= maxDepth) {
                    continue;
                }
                for (int next = 0; next < customers.length; next++) {
                    if (!canExtendForwardRelaxed(label, next)) {
                        continue;
                    }
                    double newLoad = label.load + q[next];
                    if (newLoad > ins.Q + 1e-9) {
                        continue;
                    }
                    int fromGlobal = customers[label.lastLocal];
                    int toGlobal = customers[next];
                    long nextBit = localBit[next];
                    long nextAllVisited = label.allVisitedMask | nextBit;
                    long nextRepeated = label.repeatedMask | (label.allVisitedMask & nextBit);
                    long nextMemory = (label.memoryMask & ngMaskByLocal[next]) | nextBit | capacityForbiddenMask(newLoad);
                    RelaxForwardLabel child = new RelaxForwardLabel(
                            next,
                            label.depth + 1,
                            newLoad,
                            label.travelCost + ins.c[fromGlobal][toGlobal],
                            label.dualGain + dual[next],
                            label.arcDualGain + transitionArcDual(label.lastLocal, next),
                            nextMemory,
                            nextAllVisited,
                            nextRepeated,
                            label
                    );
                    if (registerRelaxed(forwardUnextended[next], child)) {
                        forwardUnextended[next].add(child);
                    }
                }
            }
        }

        private void extendBackward(
                ArrayList<RelaxBackwardLabel>[] backwardUnextended,
                ArrayList<RelaxBackwardLabel>[] backwardExtended,
                double loadLimit,
                long startNs
        ) {
            while (true) {
                if (isOverBudget(startNs)) {
                    aborted = true;
                    return;
                }
                RelaxBackwardLabel label = pollBackward(backwardUnextended, loadLimit);
                if (label == null) {
                    return;
                }
                insertSortedByPartialReducedCost(backwardExtended[label.firstLocal], label);
                considerBackwardCompleteRoute(label);
                if (maxColumns > 0 && collected.size() >= maxColumns) {
                    return;
                }
                if (!canStillReachUsefulBackward(label)) {
                    continue;
                }
                if (label.depth >= maxDepth) {
                    continue;
                }
                for (int prev = 0; prev < customers.length; prev++) {
                    if (!canPrependBackwardRelaxed(label, prev)) {
                        continue;
                    }
                    double newLoad = label.load + q[prev];
                    if (newLoad > ins.Q + 1e-9) {
                        continue;
                    }
                    int prevGlobal = customers[prev];
                    int firstGlobal = customers[label.firstLocal];
                    long prevBit = localBit[prev];
                    long nextAllVisited = label.allVisitedMask | prevBit;
                    long nextRepeated = label.repeatedMask | (label.allVisitedMask & prevBit);
                    long nextMemory = (label.memoryMask & ngMaskByLocal[prev]) | prevBit | capacityForbiddenMask(newLoad);
                    RelaxBackwardLabel child = new RelaxBackwardLabel(
                            prev,
                            label.depth + 1,
                            newLoad,
                            ins.c[prevGlobal][firstGlobal] + label.travelCost,
                            label.dualGain + dual[prev],
                            label.arcDualGain + transitionArcDual(prev, label.firstLocal),
                            nextMemory,
                            nextAllVisited,
                            nextRepeated,
                            label
                    );
                    if (registerRelaxed(backwardUnextended[prev], child)) {
                        backwardUnextended[prev].add(child);
                    }
                }
            }
        }

        private void joinRelaxedLabels(
                ArrayList<RelaxForwardLabel>[] forwardExtended,
                ArrayList<RelaxBackwardLabel>[] backwardExtended,
                long startNs
        ) {
            for (int forwardLast = 0; forwardLast < customers.length; forwardLast++) {
                if (isOverBudget(startNs)) {
                    aborted = true;
                    return;
                }
                ArrayList<RelaxForwardLabel> fLabels = forwardExtended[forwardLast];
                if (fLabels.isEmpty()) {
                    continue;
                }
                for (int backwardFirst = 0; backwardFirst < customers.length; backwardFirst++) {
                    if (isOverBudget(startNs)) {
                        aborted = true;
                        return;
                    }
                    if (constraints != null && constraints.forbiddenArcLocal != null
                            && constraints.forbiddenArcLocal[forwardLast][backwardFirst]) {
                        continue;
                    }
                    ArrayList<RelaxBackwardLabel> bLabels = backwardExtended[backwardFirst];
                    if (bLabels.isEmpty()) {
                        continue;
                    }
                    double bridgeReducedCost = ins.c[customers[forwardLast]][customers[backwardFirst]]
                            - transitionArcDual(forwardLast, backwardFirst)
                            - dualU0;
                    double bestBackwardPartial = bLabels.get(0).partialReducedCost;
                    for (int fi = 0; fi < fLabels.size(); fi++) {
                        RelaxForwardLabel fLabel = fLabels.get(fi);
                        if (fLabel.partialReducedCost + bestBackwardPartial + bridgeReducedCost
                                >= lowerBoundThreshold() - DOM_EPS) {
                            break;
                        }
                        for (int bi = 0; bi < bLabels.size(); bi++) {
                            RelaxBackwardLabel bLabel = bLabels.get(bi);
                            if (fLabel.partialReducedCost + bLabel.partialReducedCost + bridgeReducedCost
                                    >= lowerBoundThreshold() - DOM_EPS) {
                                break;
                            }
                            if ((fLabel.memoryMask & bLabel.allVisitedMask) != 0L
                                    || (bLabel.memoryMask & fLabel.allVisitedMask) != 0L) {
                                continue;
                            }
                            double load = fLabel.load + bLabel.load;
                            if (load > ins.Q + 1e-9) {
                                continue;
                            }
                            long allVisited = fLabel.allVisitedMask | bLabel.allVisitedMask;
                            if (!isRequiredSetSatisfied(allVisited) || violatesRouteLevelPairs(allVisited)) {
                                continue;
                            }
                            long repeated = fLabel.repeatedMask | bLabel.repeatedMask | (fLabel.allVisitedMask & bLabel.allVisitedMask);
                            double fullCost = fLabel.travelCost
                                    + ins.c[customers[fLabel.lastLocal]][customers[bLabel.firstLocal]]
                                    + bLabel.travelCost;
                            double reducedCost = fLabel.partialReducedCost
                                    + bLabel.partialReducedCost
                                    + bridgeReducedCost;
                            int[] locals = combineLocals(reconstructForwardLocals(fLabel), reconstructBackwardLocals(bLabel));
                            considerAnyCandidate(reducedCost, locals);
                            if (repeated != 0L || reducedCost >= -RC_EPS) {
                                continue;
                            }
                            considerElementaryCandidate(locals, fullCost, load, reducedCost);
                            if (collected.size() >= maxColumns && maxColumns > 0) {
                                return;
                            }
                        }
                    }
                }
            }
        }

        private <T extends PartialRcState> void insertSortedByPartialReducedCost(ArrayList<T> labels, T label) {
            int lo = 0;
            int hi = labels.size();
            double key = label.partialReducedCost();
            while (lo < hi) {
                int mid = (lo + hi) >>> 1;
                if (labels.get(mid).partialReducedCost() <= key + DOM_EPS) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            labels.add(lo, label);
        }

        private boolean isOverBudget(long startNs) {
            return System.nanoTime() - startNs > BIDIRECTIONAL_RELAXED_TIME_BUDGET_NS;
        }

        private boolean canStillReachUsefulForward(RelaxForwardLabel label) {
            if (completionBound == null) {
                return true;
            }
            double completionLb = completionBound.lowerBound(label.lastLocal, label.load);
            if (Double.isNaN(completionLb)) {
                return true;
            }
            if (!isFinite(completionLb)) {
                return false;
            }
            double routeLowerBound = label.partialReducedCost + completionLb - dualU0;
            return routeLowerBound < lowerBoundThreshold() - DOM_EPS;
        }

        private boolean canStillReachUsefulBackward(RelaxBackwardLabel label) {
            if (completionBound == null) {
                return true;
            }
            double completionLb = completionBound.backwardLowerBound(label.firstLocal, label.load);
            if (Double.isNaN(completionLb)) {
                return true;
            }
            if (!isFinite(completionLb)) {
                return false;
            }
            double routeLowerBound = label.partialReducedCost + completionLb - dualU0;
            return routeLowerBound < lowerBoundThreshold() - DOM_EPS;
        }

        private double lowerBoundThreshold() {
            if (maxColumns > 0 && collected.size() >= maxColumns) {
                return collected.get(collected.size() - 1).reducedCost;
            }
            if (bestReducedCost < -RC_EPS) {
                return bestReducedCost;
            }
            return -RC_EPS;
        }

        private void considerForwardCompleteRoute(RelaxForwardLabel label) {
            if (constraints != null && constraints.forbiddenReturnLocal != null
                    && constraints.forbiddenReturnLocal[label.lastLocal]) {
                return;
            }
            if (!isRequiredSetSatisfied(label.allVisitedMask) || violatesRouteLevelPairs(label.allVisitedMask)) {
                return;
            }
            double fullCost = label.travelCost + ins.c[customers[label.lastLocal]][ins.n + 1];
            double reducedCost = fullCost - label.dualGain - label.arcDualGain - dualU0 - returnArcDual(label.lastLocal);
            int[] locals = reconstructForwardLocals(label);
            considerAnyCandidate(reducedCost, locals);
            if (label.repeatedMask != 0L || reducedCost >= -RC_EPS) {
                return;
            }
            considerElementaryCandidate(locals, fullCost, label.load, reducedCost);
        }

        private void considerBackwardCompleteRoute(RelaxBackwardLabel label) {
            if (constraints != null && constraints.forbiddenStartLocal != null
                    && constraints.forbiddenStartLocal[label.firstLocal]) {
                return;
            }
            if (!isRequiredSetSatisfied(label.allVisitedMask) || violatesRouteLevelPairs(label.allVisitedMask)) {
                return;
            }
            double fullCost = ins.c[0][customers[label.firstLocal]] + label.travelCost;
            double reducedCost = fullCost - label.dualGain - label.arcDualGain - dualU0;
            int[] locals = reconstructBackwardLocals(label);
            considerAnyCandidate(reducedCost, locals);
            if (label.repeatedMask != 0L || reducedCost >= -RC_EPS) {
                return;
            }
            considerElementaryCandidate(locals, fullCost, label.load, reducedCost);
        }

        private void considerAnyCandidate(double reducedCost, int[] locals) {
            if (reducedCost < bestAnyReducedCost - DOM_EPS) {
                bestAnyReducedCost = reducedCost;
                bestAnyLocals = locals;
            }
        }

        private void considerElementaryCandidate(int[] locals, double fullCost, double load, double reducedCost) {
            int[] globals = new int[locals.length];
            for (int idx = 0; idx < locals.length; idx++) {
                globals[idx] = customers[locals[idx]];
            }
            RouteColumn route = new RouteColumn(globals, fullCost, load);
            String key = route.key();
            if (existingKeys.contains(key) || collectedKeys.contains(key)) {
                return;
            }
            if (reducedCost < bestReducedCost - DOM_EPS) {
                bestReducedCost = reducedCost;
                bestRoute = route;
            }
            maybeCollect(new CollectedColumn(reducedCost, route, key));
        }

        private boolean canExtendForwardRelaxed(RelaxForwardLabel label, int nextLocal) {
            if (nextLocal == label.lastLocal) {
                return false;
            }
            if ((label.memoryMask & localBit[nextLocal]) != 0L) {
                return false;
            }
            if (constraints != null && constraints.forbiddenArcLocal != null
                    && constraints.forbiddenArcLocal[label.lastLocal][nextLocal]) {
                return false;
            }
            if (constraints != null && constraints.forbiddenWithLocalLong != null
                    && (label.allVisitedMask & constraints.forbiddenWithLocalLong[nextLocal]) != 0L) {
                return false;
            }
            return true;
        }

        private boolean canPrependBackwardRelaxed(RelaxBackwardLabel label, int prevLocal) {
            if (prevLocal == label.firstLocal) {
                return false;
            }
            if ((label.memoryMask & localBit[prevLocal]) != 0L) {
                return false;
            }
            if (constraints != null && constraints.forbiddenArcLocal != null
                    && constraints.forbiddenArcLocal[prevLocal][label.firstLocal]) {
                return false;
            }
            if (constraints != null && constraints.forbiddenWithLocalLong != null
                    && (label.allVisitedMask & constraints.forbiddenWithLocalLong[prevLocal]) != 0L) {
                return false;
            }
            return true;
        }

        private boolean isRequiredSetSatisfied(long visitedMask) {
            if (constraints == null || constraints.requiredWithLocalLong == null) {
                return true;
            }
            long bits = visitedMask;
            while (bits != 0L) {
                int local = Long.numberOfTrailingZeros(bits);
                if ((visitedMask & constraints.requiredWithLocalLong[local]) != constraints.requiredWithLocalLong[local]) {
                    return false;
                }
                bits ^= (bits & -bits);
            }
            return true;
        }

        private boolean violatesRouteLevelPairs(long visitedMask) {
            if (constraints == null) {
                return false;
            }
            if (constraints.forbiddenWithLocalLong != null) {
                long bits = visitedMask;
                while (bits != 0L) {
                    int local = Long.numberOfTrailingZeros(bits);
                    if ((visitedMask & constraints.forbiddenWithLocalLong[local]) != 0L) {
                        return true;
                    }
                    bits ^= (bits & -bits);
                }
            }
            return false;
        }

        private boolean registerRelaxed(ArrayList<? extends RelaxLabelState> labelsBase, RelaxLabelState candidate) {
            @SuppressWarnings("unchecked")
            ArrayList<RelaxLabelState> labels = (ArrayList<RelaxLabelState>) labelsBase;
            for (int idx = 0; idx < labels.size(); idx++) {
                RelaxLabelState incumbent = labels.get(idx);
                if (dominates(incumbent, candidate)) {
                    return false;
                }
            }
            for (int idx = 0; idx < labels.size(); idx++) {
                if (dominates(candidate, labels.get(idx))) {
                    labels.remove(idx);
                    idx--;
                }
            }
            return true;
        }

        private boolean dominates(RelaxLabelState a, RelaxLabelState b) {
            if (a.load() > b.load() + 1e-9) {
                return false;
            }
            if (a.partialReducedCost() > b.partialReducedCost() + DOM_EPS) {
                return false;
            }
            if ((a.memoryMask() & b.memoryMask()) != a.memoryMask()) {
                return false;
            }
            return (a.repeatedMask() & ~b.repeatedMask()) == 0L;
        }

        private RelaxForwardLabel pollForward(ArrayList<RelaxForwardLabel>[] labels, double loadLimit) {
            for (int local = 0; local < labels.length; local++) {
                ArrayList<RelaxForwardLabel> list = labels[local];
                for (int idx = 0; idx < list.size(); idx++) {
                    RelaxForwardLabel label = list.get(idx);
                    if (label.load <= loadLimit + 1e-9) {
                        list.remove(idx);
                        return label;
                    }
                }
            }
            return null;
        }

        private RelaxBackwardLabel pollBackward(ArrayList<RelaxBackwardLabel>[] labels, double loadLimit) {
            for (int local = 0; local < labels.length; local++) {
                ArrayList<RelaxBackwardLabel> list = labels[local];
                for (int idx = 0; idx < list.size(); idx++) {
                    RelaxBackwardLabel label = list.get(idx);
                    if (label.load <= loadLimit + 1e-9) {
                        list.remove(idx);
                        return label;
                    }
                }
            }
            return null;
        }

        private int countRemaining(ArrayList<?>[] labels) {
            int total = 0;
            for (int local = 0; local < labels.length; local++) {
                total += labels[local].size();
            }
            return total;
        }

        private long[] buildInitialNgMasks() {
            int m = customers.length;
            int ngSize = Math.max(1, Math.min(recommendedNgSize(m), m));
            long[] masks = new long[m];
            for (int i = 0; i < m; i++) {
                masks[i] = localBit[i];
                boolean[] selected = new boolean[m];
                selected[i] = true;
                for (int count = 1; count < ngSize; count++) {
                    double best = Double.POSITIVE_INFINITY;
                    int bestIdx = -1;
                    int gi = customers[i];
                    for (int j = 0; j < m; j++) {
                        if (selected[j]) {
                            continue;
                        }
                        int gj = customers[j];
                        double dist = ins.c[gi][gj];
                        if (dist < best - 1e-12) {
                            best = dist;
                            bestIdx = j;
                        }
                    }
                    if (bestIdx < 0) {
                        break;
                    }
                    selected[bestIdx] = true;
                    masks[i] |= localBit[bestIdx];
                }
            }
            return masks;
        }

        private int recommendedNgSize(int customerCount) {
            return Math.max(DEFAULT_NG_SIZE, Math.min(customerCount / 3, 12));
        }

        private long capacityForbiddenMask(double load) {
            double rem = ins.Q - load;
            long mask = 0L;
            for (int local = 0; local < q.length; local++) {
                if (q[local] > rem + 1e-9) {
                    mask |= localBit[local];
                }
            }
            return mask;
        }

        private double transitionArcDual(int fromLocal, int toLocal) {
            if (constraints == null || constraints.arcDualLocal == null) {
                return 0.0;
            }
            return constraints.arcDualLocal[fromLocal][toLocal];
        }

        private double returnArcDual(int lastLocal) {
            if (constraints == null || constraints.returnArcDualLocal == null) {
                return 0.0;
            }
            return constraints.returnArcDualLocal[lastLocal];
        }

        private int[] reconstructForwardLocals(RelaxForwardLabel label) {
            int[] locals = new int[label.depth];
            RelaxForwardLabel cur = label;
            for (int pos = label.depth - 1; pos >= 0; pos--) {
                locals[pos] = cur.lastLocal;
                cur = cur.pred;
            }
            return locals;
        }

        private int[] reconstructBackwardLocals(RelaxBackwardLabel label) {
            int[] locals = new int[label.depth];
            RelaxBackwardLabel cur = label;
            for (int pos = 0; pos < label.depth; pos++) {
                locals[pos] = cur.firstLocal;
                cur = cur.pred;
            }
            return locals;
        }

        private int[] combineLocals(int[] left, int[] right) {
            int[] out = new int[left.length + right.length];
            System.arraycopy(left, 0, out, 0, left.length);
            System.arraycopy(right, 0, out, left.length, right.length);
            return out;
        }

        private boolean updateNgSetsFromCycles(int[] locals) {
            boolean changed = false;
            for (int i = 0; i < locals.length - 1; i++) {
                for (int j = i + 1; j < locals.length; j++) {
                    if (locals[j] != locals[i]) {
                        continue;
                    }
                    int start = locals[i];
                    for (int k = i + 1; k < j; k++) {
                        int node = locals[k];
                        long before = ngMaskByLocal[node];
                        ngMaskByLocal[node] |= localBit[start];
                        changed |= before != ngMaskByLocal[node];
                    }
                    break;
                }
            }
            return changed;
        }

        private void maybeCollect(CollectedColumn candidate) {
            if (candidate == null || collectedKeys.contains(candidate.key)) {
                return;
            }
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

    private interface RelaxLabelState extends PartialRcState {
        double load();
        long memoryMask();
        long repeatedMask();
    }

    private static final class RelaxForwardLabel implements RelaxLabelState {
        final int lastLocal;
        final int depth;
        final double load;
        final double travelCost;
        final double dualGain;
        final double arcDualGain;
        final double partialReducedCost;
        final long memoryMask;
        final long allVisitedMask;
        final long repeatedMask;
        final RelaxForwardLabel pred;

        RelaxForwardLabel(
                int lastLocal,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                double arcDualGain,
                long memoryMask,
                long allVisitedMask,
                long repeatedMask,
                RelaxForwardLabel pred
        ) {
            this.lastLocal = lastLocal;
            this.depth = depth;
            this.load = load;
            this.travelCost = travelCost;
            this.dualGain = dualGain;
            this.arcDualGain = arcDualGain;
            this.partialReducedCost = travelCost - dualGain - arcDualGain;
            this.memoryMask = memoryMask;
            this.allVisitedMask = allVisitedMask;
            this.repeatedMask = repeatedMask;
            this.pred = pred;
        }

        @Override
        public double partialReducedCost() {
            return partialReducedCost;
        }

        @Override
        public double load() {
            return load;
        }

        @Override
        public long memoryMask() {
            return memoryMask;
        }

        @Override
        public long repeatedMask() {
            return repeatedMask;
        }
    }

    private static final class RelaxBackwardLabel implements RelaxLabelState {
        final int firstLocal;
        final int depth;
        final double load;
        final double travelCost;
        final double dualGain;
        final double arcDualGain;
        final double partialReducedCost;
        final long memoryMask;
        final long allVisitedMask;
        final long repeatedMask;
        final RelaxBackwardLabel pred;

        RelaxBackwardLabel(
                int firstLocal,
                int depth,
                double load,
                double travelCost,
                double dualGain,
                double arcDualGain,
                long memoryMask,
                long allVisitedMask,
                long repeatedMask,
                RelaxBackwardLabel pred
        ) {
            this.firstLocal = firstLocal;
            this.depth = depth;
            this.load = load;
            this.travelCost = travelCost;
            this.dualGain = dualGain;
            this.arcDualGain = arcDualGain;
            this.partialReducedCost = travelCost - dualGain - arcDualGain;
            this.memoryMask = memoryMask;
            this.allVisitedMask = allVisitedMask;
            this.repeatedMask = repeatedMask;
            this.pred = pred;
        }

        @Override
        public double partialReducedCost() {
            return partialReducedCost;
        }

        @Override
        public double load() {
            return load;
        }

        @Override
        public long memoryMask() {
            return memoryMask;
        }

        @Override
        public long repeatedMask() {
            return repeatedMask;
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
        final double partialReducedCost;
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
            this.partialReducedCost = travelCost - dualGain;
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
            if (incumbent != null && incumbent.partialReducedCost <= label.partialReducedCost + DOM_EPS) {
                return false;
            }
            if (isDominatedByProperSubset(label, map)) {
                return false;
            }
            map.put(label.visitedMask, label);
            removeDominatedSupersets(label, map);
            return true;
        }

        private boolean isCurrentBest(ScenarioLongLabel label) {
            return bestByMaskAtLast[label.lastScenario].get(label.visitedMask) == label;
        }

        private boolean isDominatedByProperSubset(ScenarioLongLabel label, HashMap<Long, ScenarioLongLabel> map) {
            if (label.depth <= 1) {
                return false;
            }
            long lastBit = scenarioBit[label.lastScenario];
            long free = label.visitedMask ^ lastBit;
            if (!shouldEnumerateSubsetDominance(free)) {
                return false;
            }
            long subFree = free;
            while (subFree != 0L) {
                long subsetMask = subFree | lastBit;
                ScenarioLongLabel subset = map.get(subsetMask);
                if (subset != null && subset.partialReducedCost <= label.partialReducedCost + DOM_EPS) {
                    return true;
                }
                subFree = (subFree - 1L) & free;
            }
            return false;
        }

        private void removeDominatedSupersets(ScenarioLongLabel label, HashMap<Long, ScenarioLongLabel> map) {
            long extraPool = fullMask & ~label.visitedMask;
            int extraCount = Long.bitCount(extraPool);
            if (extraCount == 0 || extraCount > SUPERSET_CLEANUP_MAX_EXTRA_BITS) {
                return;
            }
            long add = extraPool;
            while (add != 0L) {
                long supersetMask = label.visitedMask | add;
                ScenarioLongLabel other = map.get(supersetMask);
                if (other != null && label.partialReducedCost <= other.partialReducedCost + DOM_EPS) {
                    map.remove(supersetMask);
                }
                add = (add - 1L) & extraPool;
            }
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
