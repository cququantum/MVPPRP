package lbbdModel;

import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import instance.Instance;
import lbbdModel.rmp.PricingEspprcSolver;
import lbbdModel.rmp.RouteColumn;
import model.CplexConfig;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashSet;
import java.util.PriorityQueue;

/**
 * Exact CVRP subproblem solver by branch-and-price.
 *
 * Master problem:
 * - set-partitioning over route columns
 * - cover each active customer exactly once
 * - use at most K vehicles
 *
 * Pricing:
 * - exact ESPPRC enumeration already implemented in {@link PricingEspprcSolver}
 *
 * Branching:
 * - Ryan-Foster style on customer pairs:
 *   together: customers i and j must belong to the same route
 *   separate: customers i and j cannot belong to the same route
 */
public final class BranchAndPriceCvrpSolver {
    private static final double EPS = 1e-7;
    private static final double ART_EPS = 1e-7;
    private static final boolean LOG_TO_CONSOLE = false;
    private static final long NANOS_PER_SECOND = 1_000_000_000L;

    public static final class Result {
        public final boolean feasible;
        public final boolean optimal;
        public final boolean provenInfeasible;
        public final String status;
        public final double objective;
        public final int exploredNodes;
        public final int generatedColumns;

        Result(
                boolean feasible,
                boolean optimal,
                boolean provenInfeasible,
                String status,
                double objective,
                int exploredNodes,
                int generatedColumns
        ) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.provenInfeasible = provenInfeasible;
            this.status = status;
            this.objective = objective;
            this.exploredNodes = exploredNodes;
            this.generatedColumns = generatedColumns;
        }
    }

    private static final class ActiveSet {
        final int[] activeCustomers;
        final double[] qLocal;
        final int[] localIndexByGlobal;
        final double totalLoad;

        ActiveSet(int[] activeCustomers, double[] qLocal, int[] localIndexByGlobal, double totalLoad) {
            this.activeCustomers = activeCustomers;
            this.qLocal = qLocal;
            this.localIndexByGlobal = localIndexByGlobal;
            this.totalLoad = totalLoad;
        }

        int size() {
            return activeCustomers.length;
        }
    }

    private static final class BranchPair {
        final int customerA;
        final int customerB;
        final double sameRouteValue;

        BranchPair(int customerA, int customerB, double sameRouteValue) {
            this.customerA = customerA;
            this.customerB = customerB;
            this.sameRouteValue = sameRouteValue;
        }
    }

    private static final class BranchNode {
        final HashSet<Long> togetherPairs;
        final HashSet<Long> separatePairs;
        final ArrayList<RouteColumn> seedColumns;
        final int depth;

        BranchNode(
                HashSet<Long> togetherPairs,
                HashSet<Long> separatePairs,
                ArrayList<RouteColumn> seedColumns,
                int depth
        ) {
            this.togetherPairs = togetherPairs;
            this.separatePairs = separatePairs;
            this.seedColumns = seedColumns;
            this.depth = depth;
        }

        BranchNode childTogether(BranchPair pair, ArrayList<RouteColumn> filteredSeedColumns) {
            long key = pairKey(pair.customerA, pair.customerB);
            if (separatePairs.contains(key)) {
                return null;
            }
            HashSet<Long> newTogether = new HashSet<Long>(togetherPairs);
            newTogether.add(key);
            return new BranchNode(newTogether, new HashSet<Long>(separatePairs), filteredSeedColumns, depth + 1);
        }

        BranchNode childSeparate(BranchPair pair, ArrayList<RouteColumn> filteredSeedColumns) {
            long key = pairKey(pair.customerA, pair.customerB);
            if (togetherPairs.contains(key)) {
                return null;
            }
            HashSet<Long> newSeparate = new HashSet<Long>(separatePairs);
            newSeparate.add(key);
            return new BranchNode(new HashSet<Long>(togetherPairs), newSeparate, filteredSeedColumns, depth + 1);
        }
    }

    private static final class RouteVarData {
        final RouteColumn route;
        final IloNumVar var;
        final int[] localCustomers;

        RouteVarData(RouteColumn route, IloNumVar var, int[] localCustomers) {
            this.route = route;
            this.var = var;
            this.localCustomers = localCustomers;
        }
    }

    private static final class ActiveRouteValue {
        final RouteColumn route;
        final double value;
        final int[] localCustomers;

        ActiveRouteValue(RouteColumn route, double value, int[] localCustomers) {
            this.route = route;
            this.value = value;
            this.localCustomers = localCustomers;
        }
    }

    private static final class NodeLpResult {
        final boolean feasible;
        final boolean optimal;
        final boolean provenInfeasible;
        final String status;
        final double objective;
        final boolean integral;
        final BranchPair branchPair;
        final ArrayList<RouteColumn> generatedRoutes;
        final int generatedColumns;

        NodeLpResult(
                boolean feasible,
                boolean optimal,
                boolean provenInfeasible,
                String status,
                double objective,
                boolean integral,
                BranchPair branchPair,
                ArrayList<RouteColumn> generatedRoutes,
                int generatedColumns
        ) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.provenInfeasible = provenInfeasible;
            this.status = status;
            this.objective = objective;
            this.integral = integral;
            this.branchPair = branchPair;
            this.generatedRoutes = generatedRoutes;
            this.generatedColumns = generatedColumns;
        }
    }

    private static final class PendingNode {
        final BranchNode node;
        final NodeLpResult lpResult;

        PendingNode(BranchNode node, NodeLpResult lpResult) {
            this.node = node;
            this.lpResult = lpResult;
        }
    }

    private final PricingEspprcSolver pricingSolver = new PricingEspprcSolver();

    public Result solve(Instance ins, int t, double[] qBar, int[] zBar) {
        long startNs = System.nanoTime();
        long deadlineNs = startNs + Math.max(1L, Math.round(CplexConfig.TIME_LIMIT_SEC * NANOS_PER_SECOND));
        ActiveSet active = buildActiveSet(ins, qBar, zBar);
        if (active == null) {
            return new Result(false, false, true, "InfeasibleByPositivePickupWithoutVisit", Double.NaN, 0, 0);
        }
        if (active.size() == 0) {
            return new Result(true, true, false, "TrivialZeroBP", 0.0, 0, 0);
        }
        if (active.totalLoad > ins.K * ins.Q + EPS) {
            return new Result(false, false, true, "InfeasibleByTotalCapacity", Double.NaN, 0, 0);
        }
        for (int local = 0; local < active.qLocal.length; local++) {
            if (active.qLocal[local] > ins.Q + EPS) {
                return new Result(false, false, true, "InfeasibleBySingleCustomerCapacity", Double.NaN, 0, 0);
            }
        }

        ArrayList<RouteColumn> rootSeeds = buildSingletonSeeds(ins, active);
        BranchNode root = new BranchNode(new HashSet<Long>(), new HashSet<Long>(), rootSeeds, 0);

        int exploredNodes = 0;
        int totalGeneratedColumns = 0;
        NodeLpResult rootLp = solveNodeLp(ins, t, active, root, deadlineNs);
        exploredNodes++;
        totalGeneratedColumns += rootLp.generatedColumns;

        if (!rootLp.optimal) {
            return new Result(false, false, rootLp.provenInfeasible, rootLp.status, Double.NaN, exploredNodes, totalGeneratedColumns);
        }
        if (!rootLp.feasible) {
            return new Result(false, false, true, rootLp.status, Double.NaN, exploredNodes, totalGeneratedColumns);
        }
        if (rootLp.integral) {
            return new Result(true, true, false, "OptimalBP(root)", rootLp.objective, exploredNodes, totalGeneratedColumns);
        }

        double incumbent = Double.POSITIVE_INFINITY;
        PriorityQueue<PendingNode> queue = new PriorityQueue<PendingNode>(
                Comparator.comparingDouble((PendingNode x) -> x.lpResult.objective).thenComparingInt(x -> x.node.depth)
        );
        queue.add(new PendingNode(root, rootLp));

        while (!queue.isEmpty()) {
            PendingNode pending = queue.poll();
            if (pending.lpResult.objective >= incumbent - EPS) {
                continue;
            }

            BranchPair pair = pending.lpResult.branchPair;
            if (pair == null) {
                incumbent = Math.min(incumbent, pending.lpResult.objective);
                continue;
            }

            ArrayList<RouteColumn> togetherSeeds = filterAllowedColumns(pending.lpResult.generatedRoutes, pending.node, pair, true);
            BranchNode togetherChild = pending.node.childTogether(pair, togetherSeeds);
            if (togetherChild != null) {
                if (isPastDeadline(deadlineNs)) {
                    return timeoutResult(exploredNodes, totalGeneratedColumns, incumbent);
                }
                NodeLpResult childLp = solveNodeLp(ins, t, active, togetherChild, deadlineNs);
                exploredNodes++;
                totalGeneratedColumns += childLp.generatedColumns;
                if (childLp.optimal) {
                    if (childLp.feasible) {
                        if (childLp.integral) {
                            incumbent = Math.min(incumbent, childLp.objective);
                        } else if (childLp.objective < incumbent - EPS) {
                            queue.add(new PendingNode(togetherChild, childLp));
                        }
                    }
                } else {
                    return new Result(false, false, childLp.provenInfeasible, childLp.status, Double.NaN, exploredNodes, totalGeneratedColumns);
                }
            }

            ArrayList<RouteColumn> separateSeeds = filterAllowedColumns(pending.lpResult.generatedRoutes, pending.node, pair, false);
            BranchNode separateChild = pending.node.childSeparate(pair, separateSeeds);
            if (separateChild != null) {
                if (isPastDeadline(deadlineNs)) {
                    return timeoutResult(exploredNodes, totalGeneratedColumns, incumbent);
                }
                NodeLpResult childLp = solveNodeLp(ins, t, active, separateChild, deadlineNs);
                exploredNodes++;
                totalGeneratedColumns += childLp.generatedColumns;
                if (childLp.optimal) {
                    if (childLp.feasible) {
                        if (childLp.integral) {
                            incumbent = Math.min(incumbent, childLp.objective);
                        } else if (childLp.objective < incumbent - EPS) {
                            queue.add(new PendingNode(separateChild, childLp));
                        }
                    }
                } else {
                    return new Result(false, false, childLp.provenInfeasible, childLp.status, Double.NaN, exploredNodes, totalGeneratedColumns);
                }
            }
        }

        if (Double.isInfinite(incumbent) || Double.isNaN(incumbent)) {
            return new Result(false, false, false, "BranchPriceNoIncumbent", Double.NaN, exploredNodes, totalGeneratedColumns);
        }
        return new Result(true, true, false, "OptimalBP", incumbent, exploredNodes, totalGeneratedColumns);
    }

    private NodeLpResult solveNodeLp(Instance ins, int t, ActiveSet active, BranchNode node, long deadlineNs) {
        if (isPastDeadline(deadlineNs)) {
            return new NodeLpResult(false, false, false, "BP_TimeLimit", Double.NaN, false, null, node.seedColumns, 0);
        }
        final int customerCount = active.size();
        final PricingEspprcSolver.RoutePricingConstraints pricingConstraints = buildPricingConstraints(active, node);
        final double artificialPenalty = computeArtificialPenalty(ins, customerCount);

        try (IloCplex cplex = new IloCplex()) {
            double remainingSec = remainingSeconds(deadlineNs);
            if (remainingSec <= 0.0) {
                return new NodeLpResult(false, false, false, "BP_TimeLimit", Double.NaN, false, null, node.seedColumns, 0);
            }
            configureLp(cplex, remainingSec);

            IloObjective obj = cplex.addMinimize();
            IloRange[] coverEq = new IloRange[customerCount];
            for (int local = 0; local < customerCount; local++) {
                coverEq[local] = cplex.addEq(cplex.linearNumExpr(), 1.0, "BP_Cover_" + t + "_" + local);
            }
            IloRange vehicleLimit = cplex.addLe(cplex.linearNumExpr(), ins.K, "BP_Vehicle_" + t);

            IloNumVar[] artificial = new IloNumVar[customerCount];
            for (int local = 0; local < customerCount; local++) {
                IloColumn col = cplex.column(obj, artificialPenalty).and(cplex.column(coverEq[local], 1.0));
                artificial[local] = cplex.numVar(col, 0.0, Double.MAX_VALUE, "bp_art_" + t + "_" + local);
            }

            HashSet<String> existingRouteKeys = new HashSet<String>();
            ArrayList<RouteVarData> routeVars = new ArrayList<RouteVarData>();
            ArrayList<RouteColumn> generatedRoutes = new ArrayList<RouteColumn>();
            int generatedColumns = 0;

            for (int idx = 0; idx < node.seedColumns.size(); idx++) {
                RouteColumn route = node.seedColumns.get(idx);
                if (!isRouteAllowed(route, node)) {
                    continue;
                }
                if (addRouteColumn(cplex, obj, coverEq, vehicleLimit, active.localIndexByGlobal, route, routeVars, existingRouteKeys,
                        "bp_seed_" + t + "_" + idx)) {
                    generatedRoutes.add(route);
                }
            }

            for (int local = 0; local < customerCount; local++) {
                RouteColumn singleton = new RouteColumn(
                        new int[]{active.activeCustomers[local]},
                        ins.c[0][active.activeCustomers[local]] + ins.c[active.activeCustomers[local]][ins.n + 1],
                        active.qLocal[local]
                );
                if (!isRouteAllowed(singleton, node)) {
                    continue;
                }
                if (addRouteColumn(cplex, obj, coverEq, vehicleLimit, active.localIndexByGlobal, singleton, routeVars, existingRouteKeys,
                        "bp_single_" + t + "_" + local)) {
                    generatedRoutes.add(singleton);
                }
            }

            while (true) {
                if (isPastDeadline(deadlineNs)) {
                    return new NodeLpResult(false, false, false, "BP_TimeLimit", Double.NaN, false, null, generatedRoutes, generatedColumns);
                }
                boolean solved = cplex.solve();
                String status = cplex.getStatus().toString();
                if (!solved || !status.startsWith("Optimal")) {
                    return new NodeLpResult(false, false, false, "BP_LP_" + status, Double.NaN, false, null, generatedRoutes, generatedColumns);
                }

                double[] dualCover = new double[customerCount];
                for (int local = 0; local < customerCount; local++) {
                    dualCover[local] = cplex.getDual(coverEq[local]);
                }
                double dualVehicle = cplex.getDual(vehicleLimit);

                PricingEspprcSolver.Result pricing = pricingSolver.findNegativeRoutes(
                        ins,
                        active.activeCustomers,
                        active.qLocal,
                        dualCover,
                        dualVehicle,
                        existingRouteKeys,
                        maxColumnsPerPricing(customerCount),
                        pricingConstraints
                );

                if (isPastDeadline(deadlineNs)) {
                    return new NodeLpResult(false, false, false, "BP_TimeLimit", Double.NaN, false, null, generatedRoutes, generatedColumns);
                }
                if (!pricing.foundNegativeColumn) {
                    break;
                }

                ArrayList<RouteColumn> newRoutes = pricing.routes;
                if (newRoutes == null || newRoutes.isEmpty()) {
                    if (addRouteColumn(cplex, obj, coverEq, vehicleLimit, active.localIndexByGlobal, pricing.route, routeVars, existingRouteKeys,
                            "bp_cg_" + t + "_" + generatedColumns)) {
                        generatedRoutes.add(pricing.route);
                        generatedColumns++;
                    }
                } else {
                    for (int idx = 0; idx < newRoutes.size(); idx++) {
                        RouteColumn route = newRoutes.get(idx);
                        if (addRouteColumn(cplex, obj, coverEq, vehicleLimit, active.localIndexByGlobal, route, routeVars, existingRouteKeys,
                                "bp_cg_" + t + "_" + generatedColumns)) {
                            generatedRoutes.add(route);
                            generatedColumns++;
                        }
                    }
                }
            }

            for (int local = 0; local < customerCount; local++) {
                if (cplex.getValue(artificial[local]) > ART_EPS) {
                    return new NodeLpResult(false, true, true, "BP_NodeInfeasible", Double.NaN, false, null, generatedRoutes, generatedColumns);
                }
            }

            ArrayList<ActiveRouteValue> positiveRoutes = new ArrayList<ActiveRouteValue>();
            for (int idx = 0; idx < routeVars.size(); idx++) {
                RouteVarData data = routeVars.get(idx);
                double value = cplex.getValue(data.var);
                if (value > EPS) {
                    positiveRoutes.add(new ActiveRouteValue(data.route, value, data.localCustomers));
                }
            }

            BranchPair branchPair = chooseBranchPair(active, positiveRoutes);
            boolean integral = branchPair == null;
            return new NodeLpResult(
                    true,
                    true,
                    false,
                    "BP_NodeOptimal",
                    cplex.getObjValue(),
                    integral,
                    branchPair,
                    generatedRoutes,
                    generatedColumns
            );
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve branch-and-price node for CVRP(t=" + t + ")", e);
        }
    }

    private static ActiveSet buildActiveSet(Instance ins, double[] qBar, int[] zBar) {
        ArrayList<Integer> customers = new ArrayList<Integer>();
        ArrayList<Double> qList = new ArrayList<Double>();
        int[] localIndexByGlobal = new int[ins.n + 1];
        Arrays.fill(localIndexByGlobal, -1);
        double totalLoad = 0.0;

        for (int i = 1; i <= ins.n; i++) {
            if (zBar[i] != 0) {
                localIndexByGlobal[i] = customers.size();
                customers.add(i);
                qList.add(qBar[i]);
                totalLoad += qBar[i];
            } else if (qBar[i] > EPS) {
                return null;
            }
        }

        int[] activeCustomers = new int[customers.size()];
        double[] qLocal = new double[qList.size()];
        for (int local = 0; local < customers.size(); local++) {
            activeCustomers[local] = customers.get(local);
            qLocal[local] = qList.get(local);
        }
        return new ActiveSet(activeCustomers, qLocal, localIndexByGlobal, totalLoad);
    }

    private static ArrayList<RouteColumn> buildSingletonSeeds(Instance ins, ActiveSet active) {
        ArrayList<RouteColumn> seeds = new ArrayList<RouteColumn>();
        for (int local = 0; local < active.size(); local++) {
            int customer = active.activeCustomers[local];
            seeds.add(new RouteColumn(
                    new int[]{customer},
                    ins.c[0][customer] + ins.c[customer][ins.n + 1],
                    active.qLocal[local]
            ));
        }
        return seeds;
    }

    private static ArrayList<RouteColumn> filterAllowedColumns(
            ArrayList<RouteColumn> columns,
            BranchNode parentNode,
            BranchPair pair,
            boolean together
    ) {
        ArrayList<RouteColumn> filtered = new ArrayList<RouteColumn>();
        HashSet<Long> togetherPairs = new HashSet<Long>(parentNode.togetherPairs);
        HashSet<Long> separatePairs = new HashSet<Long>(parentNode.separatePairs);
        long key = pairKey(pair.customerA, pair.customerB);
        if (together) {
            togetherPairs.add(key);
        } else {
            separatePairs.add(key);
        }
        BranchNode childView = new BranchNode(togetherPairs, separatePairs, new ArrayList<RouteColumn>(), parentNode.depth + 1);
        for (int idx = 0; idx < columns.size(); idx++) {
            RouteColumn route = columns.get(idx);
            if (isRouteAllowed(route, childView)) {
                filtered.add(route);
            }
        }
        return filtered;
    }

    private static BranchPair chooseBranchPair(ActiveSet active, ArrayList<ActiveRouteValue> positiveRoutes) {
        int m = active.size();
        if (m <= 1) {
            return null;
        }
        double[][] sameRoute = new double[m][m];
        for (int idx = 0; idx < positiveRoutes.size(); idx++) {
            ActiveRouteValue route = positiveRoutes.get(idx);
            int[] locals = route.localCustomers;
            for (int a = 0; a < locals.length; a++) {
                for (int b = a + 1; b < locals.length; b++) {
                    int i = locals[a];
                    int j = locals[b];
                    sameRoute[i][j] += route.value;
                }
            }
        }

        BranchPair best = null;
        double bestScore = Double.POSITIVE_INFINITY;
        for (int i = 0; i < m; i++) {
            for (int j = i + 1; j < m; j++) {
                double val = sameRoute[i][j];
                if (val > EPS && val < 1.0 - EPS) {
                    double score = Math.abs(0.5 - val);
                    if (score < bestScore - 1e-12) {
                        bestScore = score;
                        best = new BranchPair(active.activeCustomers[i], active.activeCustomers[j], val);
                    }
                }
            }
        }
        return best;
    }

    private static PricingEspprcSolver.RoutePricingConstraints buildPricingConstraints(ActiveSet active, BranchNode node) {
        if (node.togetherPairs.isEmpty() && node.separatePairs.isEmpty()) {
            return null;
        }

        int m = active.size();
        if (m <= 62) {
            long[] required = new long[m];
            long[] forbidden = new long[m];
            populateLongConstraints(required, forbidden, active, node.togetherPairs, true);
            populateLongConstraints(required, forbidden, active, node.separatePairs, false);
            return PricingEspprcSolver.RoutePricingConstraints.fromLongMasks(required, forbidden);
        }

        BitSet[] required = new BitSet[m];
        BitSet[] forbidden = new BitSet[m];
        for (int local = 0; local < m; local++) {
            required[local] = new BitSet(m);
            forbidden[local] = new BitSet(m);
        }
        populateBitSetConstraints(required, forbidden, active, node.togetherPairs, true);
        populateBitSetConstraints(required, forbidden, active, node.separatePairs, false);
        return PricingEspprcSolver.RoutePricingConstraints.fromBitSets(required, forbidden);
    }

    private static void populateLongConstraints(
            long[] required,
            long[] forbidden,
            ActiveSet active,
            HashSet<Long> pairs,
            boolean together
    ) {
        for (Long keyObj : pairs) {
            long key = keyObj.longValue();
            int a = (int) (key >>> 32);
            int b = (int) key;
            int localA = (a >= 0 && a < active.localIndexByGlobal.length) ? active.localIndexByGlobal[a] : -1;
            int localB = (b >= 0 && b < active.localIndexByGlobal.length) ? active.localIndexByGlobal[b] : -1;
            if (localA < 0 || localB < 0) {
                continue;
            }
            if (together) {
                required[localA] |= (1L << localB);
                required[localB] |= (1L << localA);
            } else {
                forbidden[localA] |= (1L << localB);
                forbidden[localB] |= (1L << localA);
            }
        }
    }

    private static void populateBitSetConstraints(
            BitSet[] required,
            BitSet[] forbidden,
            ActiveSet active,
            HashSet<Long> pairs,
            boolean together
    ) {
        for (Long keyObj : pairs) {
            long key = keyObj.longValue();
            int a = (int) (key >>> 32);
            int b = (int) key;
            int localA = (a >= 0 && a < active.localIndexByGlobal.length) ? active.localIndexByGlobal[a] : -1;
            int localB = (b >= 0 && b < active.localIndexByGlobal.length) ? active.localIndexByGlobal[b] : -1;
            if (localA < 0 || localB < 0) {
                continue;
            }
            if (together) {
                required[localA].set(localB);
                required[localB].set(localA);
            } else {
                forbidden[localA].set(localB);
                forbidden[localB].set(localA);
            }
        }
    }

    private static boolean isRouteAllowed(RouteColumn route, BranchNode node) {
        if (node.togetherPairs.isEmpty() && node.separatePairs.isEmpty()) {
            return true;
        }
        BitSet routeCustomers = new BitSet();
        for (int idx = 0; idx < route.globalCustomersInOrder.length; idx++) {
            routeCustomers.set(route.globalCustomersInOrder[idx]);
        }
        for (Long keyObj : node.separatePairs) {
            long key = keyObj.longValue();
            int a = (int) (key >>> 32);
            int b = (int) key;
            if (routeCustomers.get(a) && routeCustomers.get(b)) {
                return false;
            }
        }
        for (Long keyObj : node.togetherPairs) {
            long key = keyObj.longValue();
            int a = (int) (key >>> 32);
            int b = (int) key;
            if (routeCustomers.get(a) != routeCustomers.get(b)) {
                return false;
            }
        }
        return true;
    }

    private static boolean addRouteColumn(
            IloCplex cplex,
            IloObjective obj,
            IloRange[] coverEq,
            IloRange vehicleLimit,
            int[] localIndexByGlobal,
            RouteColumn route,
            ArrayList<RouteVarData> routeVars,
            HashSet<String> existingRouteKeys,
            String varName
    ) throws IloException {
        if (route == null) {
            return false;
        }
        String key = route.key();
        if (existingRouteKeys.contains(key)) {
            return false;
        }
        existingRouteKeys.add(key);

        IloColumn col = cplex.column(obj, route.cost).and(cplex.column(vehicleLimit, 1.0));
        boolean[] covered = new boolean[coverEq.length];
        int[] tmp = new int[route.globalCustomersInOrder.length];
        int count = 0;
        for (int pos = 0; pos < route.globalCustomersInOrder.length; pos++) {
            int global = route.globalCustomersInOrder[pos];
            int local = (global >= 0 && global < localIndexByGlobal.length) ? localIndexByGlobal[global] : -1;
            if (local >= 0 && !covered[local]) {
                col = col.and(cplex.column(coverEq[local], 1.0));
                covered[local] = true;
                tmp[count++] = local;
            }
        }
        int[] locals = new int[count];
        System.arraycopy(tmp, 0, locals, 0, count);
        IloNumVar var = cplex.numVar(col, 0.0, Double.MAX_VALUE, "xi_bp_" + varName);
        routeVars.add(new RouteVarData(route, var, locals));
        return true;
    }

    private static int maxColumnsPerPricing(int customerCount) {
        if (customerCount >= 25) {
            return 24;
        }
        if (customerCount >= 15) {
            return 18;
        }
        if (customerCount >= 8) {
            return 12;
        }
        return 8;
    }

    private static long pairKey(int a, int b) {
        int x = Math.min(a, b);
        int y = Math.max(a, b);
        return (((long) x) << 32) | (y & 0xffffffffL);
    }

    private static double computeArtificialPenalty(Instance ins, int customerCount) {
        double maxArc = 0.0;
        for (int i = 0; i < ins.nodeCount; i++) {
            for (int j = 0; j < ins.nodeCount; j++) {
                if (i != j && ins.c[i][j] > maxArc) {
                    maxArc = ins.c[i][j];
                }
            }
        }
        return 1_000_000.0 + maxArc * (customerCount + 2);
    }

    private static Result timeoutResult(int exploredNodes, int totalGeneratedColumns, double incumbent) {
        if (Double.isFinite(incumbent)) {
            return new Result(true, false, false, "BP_TimeLimitWithIncumbent", incumbent, exploredNodes, totalGeneratedColumns);
        }
        return new Result(false, false, false, "BP_TimeLimit", Double.NaN, exploredNodes, totalGeneratedColumns);
    }

    private static boolean isPastDeadline(long deadlineNs) {
        return System.nanoTime() >= deadlineNs;
    }

    private static double remainingSeconds(long deadlineNs) {
        long remainingNs = deadlineNs - System.nanoTime();
        if (remainingNs <= 0L) {
            return 0.0;
        }
        return remainingNs / (double) NANOS_PER_SECOND;
    }

    private static void configureLp(IloCplex cplex, double timeLimitSec) throws IloException {
        cplex.setParam(IloCplex.Param.TimeLimit, Math.max(1.0, timeLimitSec));
        cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Dual);
        cplex.setParam(IloCplex.Param.Threads, 1);
        if (!LOG_TO_CONSOLE) {
            cplex.setOut(null);
            cplex.setWarning(null);
        }
    }
}
