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
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.PriorityQueue;

/**
 * Exact CVRP subproblem solver by branch-and-price.
 *
 * This implementation follows the LRP BP structure:
 * - set-partitioning/covering style route master solved by CPLEX LP
 * - column generation with exact pricing
 * - vehicle-count branching first
 * - arc branching second
 * - node-local capacity-cut separation
 *
 * The surrounding LBBD interface remains unchanged.
 */
public final class BranchAndPriceCvrpSolver {
    private static final double EPS = 1e-7;
    private static final double ART_EPS = 1e-7;
    private static final double CUT_VIOL_EPS = 5e-2;
    private static final double FRAC_DISTANCE_LIMIT = 0.4;
    private static final double ARC_SUPPORT_EPS = 1e-9;
    private static final int MAX_NEW_CUTS_PER_NODE = 6;
    private static final int MAX_CONNECTED_CUT_SIZE = 10;
    private static final int MAX_STRONG_BRANCH_CANDIDATES = 8;
    private static final double MIN_STRONG_BRANCH_REMAINING_SEC = 5.0;
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

    private enum BranchType {
        LEAF,
        VEHICLE,
        ARC
    }

    private static final class ActiveSet {
        final int[] activeCustomers;
        final double[] qLocal;
        final int[] localIndexByGlobal;
        final double totalLoad;
        final double vehicleCapacity;

        ActiveSet(int[] activeCustomers, double[] qLocal, int[] localIndexByGlobal, double totalLoad, double vehicleCapacity) {
            this.activeCustomers = activeCustomers;
            this.qLocal = qLocal;
            this.localIndexByGlobal = localIndexByGlobal;
            this.totalLoad = totalLoad;
            this.vehicleCapacity = vehicleCapacity;
        }

        int size() {
            return activeCustomers.length;
        }
    }

    private static final class ArcBranch {
        final int fromNode; // 0 = depot, 1..m = local customers
        final int toNode;   // 0 = depot, 1..m = local customers
        final double value;

        ArcBranch(int fromNode, int toNode, double value) {
            this.fromNode = fromNode;
            this.toNode = toNode;
            this.value = value;
        }
    }

    private static final class BranchDecision {
        final BranchType type;
        final double vehicleValue;
        final ArcBranch arc;
        final ArrayList<RouteColumn> leftSeedColumns;
        final ArrayList<RouteColumn> rightSeedColumns;

        private BranchDecision(
                BranchType type,
                double vehicleValue,
                ArcBranch arc,
                ArrayList<RouteColumn> leftSeedColumns,
                ArrayList<RouteColumn> rightSeedColumns
        ) {
            this.type = type;
            this.vehicleValue = vehicleValue;
            this.arc = arc;
            this.leftSeedColumns = leftSeedColumns;
            this.rightSeedColumns = rightSeedColumns;
        }

        static BranchDecision leaf() {
            return new BranchDecision(BranchType.LEAF, Double.NaN, null, null, null);
        }

        static BranchDecision vehicle(double vehicleValue) {
            return new BranchDecision(BranchType.VEHICLE, vehicleValue, null, null, null);
        }

        static BranchDecision arc(ArcBranch arc) {
            return new BranchDecision(BranchType.ARC, Double.NaN, arc, null, null);
        }

        static BranchDecision arc(ArcBranch arc, ArrayList<RouteColumn> leftSeedColumns, ArrayList<RouteColumn> rightSeedColumns) {
            return new BranchDecision(BranchType.ARC, Double.NaN, arc, leftSeedColumns, rightSeedColumns);
        }
    }

    private static final class CapacityCutDef {
        final BitSet members; // customer locals 0..m-1
        final int[] locals;
        final double rhs;
        final String signature;

        CapacityCutDef(BitSet members, int[] locals, double rhs, String signature) {
            this.members = members;
            this.locals = locals;
            this.rhs = rhs;
            this.signature = signature;
        }
    }

    private static final class BranchNode {
        final int depth;
        final int vehicleLowerBound;
        final int vehicleUpperBound;
        final boolean[][] forbiddenArcs; // local node ids 0..m, 0 = depot
        final ArrayList<RouteColumn> seedColumns;
        final ArrayList<CapacityCutDef> capacityCuts;

        BranchNode(
                int depth,
                int vehicleLowerBound,
                int vehicleUpperBound,
                boolean[][] forbiddenArcs,
                ArrayList<RouteColumn> seedColumns,
                ArrayList<CapacityCutDef> capacityCuts
        ) {
            this.depth = depth;
            this.vehicleLowerBound = vehicleLowerBound;
            this.vehicleUpperBound = vehicleUpperBound;
            this.forbiddenArcs = forbiddenArcs;
            this.seedColumns = seedColumns;
            this.capacityCuts = capacityCuts;
        }

        static BranchNode root(int localNodeCount, int vehicleUpperBound, ArrayList<RouteColumn> seedColumns) {
            return new BranchNode(
                    0,
                    0,
                    vehicleUpperBound,
                    new boolean[localNodeCount + 1][localNodeCount + 1],
                    seedColumns,
                    new ArrayList<CapacityCutDef>()
            );
        }

        BranchNode childVehicleUpper(int newUpperBound, ArrayList<RouteColumn> childSeeds, ArrayList<CapacityCutDef> childCuts) {
            return new BranchNode(
                    depth + 1,
                    vehicleLowerBound,
                    newUpperBound,
                    copyForbiddenArcs(forbiddenArcs),
                    childSeeds,
                    childCuts
            );
        }

        BranchNode childVehicleLower(int newLowerBound, ArrayList<RouteColumn> childSeeds, ArrayList<CapacityCutDef> childCuts) {
            return new BranchNode(
                    depth + 1,
                    newLowerBound,
                    vehicleUpperBound,
                    copyForbiddenArcs(forbiddenArcs),
                    childSeeds,
                    childCuts
            );
        }

        BranchNode childForbidArc(ArcBranch arc, ArrayList<RouteColumn> childSeeds, ArrayList<CapacityCutDef> childCuts) {
            boolean[][] childForbid = copyForbiddenArcs(forbiddenArcs);
            forbidArc(childForbid, arc.fromNode, arc.toNode);
            return new BranchNode(depth + 1, vehicleLowerBound, vehicleUpperBound, childForbid, childSeeds, childCuts);
        }

        BranchNode childForceArc(ArcBranch arc, ArrayList<RouteColumn> childSeeds, ArrayList<CapacityCutDef> childCuts) {
            boolean[][] childForbid = copyForbiddenArcs(forbiddenArcs);
            int nodeCount = childForbid.length - 1;
            if (arc.fromNode != 0) {
                for (int to = 0; to <= nodeCount; to++) {
                    if (to != arc.toNode) {
                        forbidArc(childForbid, arc.fromNode, to);
                    }
                }
            }
            if (arc.toNode != 0) {
                for (int from = 0; from <= nodeCount; from++) {
                    if (from != arc.fromNode) {
                        forbidArc(childForbid, from, arc.toNode);
                    }
                }
            }
            return new BranchNode(depth + 1, vehicleLowerBound, vehicleUpperBound, childForbid, childSeeds, childCuts);
        }

        private static boolean[][] copyForbiddenArcs(boolean[][] source) {
            boolean[][] copy = new boolean[source.length][];
            for (int i = 0; i < source.length; i++) {
                copy[i] = Arrays.copyOf(source[i], source[i].length);
            }
            return copy;
        }

        private static void forbidArc(boolean[][] matrix, int fromNode, int toNode) {
            if (fromNode < 0 || toNode < 0 || fromNode >= matrix.length || toNode >= matrix.length) {
                return;
            }
            if (fromNode == toNode) {
                return;
            }
            matrix[fromNode][toNode] = true;
        }
    }

    private static final class RouteVarData {
        final RouteColumn route;
        final IloNumVar var;
        final int[] localCustomers;
        final int[] localNodePath;

        RouteVarData(RouteColumn route, IloNumVar var, int[] localCustomers, int[] localNodePath) {
            this.route = route;
            this.var = var;
            this.localCustomers = localCustomers;
            this.localNodePath = localNodePath;
        }
    }

    private static final class ActiveRouteValue {
        final RouteColumn route;
        final double value;
        final int[] localCustomers;
        final int[] localNodePath;

        ActiveRouteValue(RouteColumn route, double value, int[] localCustomers, int[] localNodePath) {
            this.route = route;
            this.value = value;
            this.localCustomers = localCustomers;
            this.localNodePath = localNodePath;
        }
    }

    private static final class CutRowData {
        final CapacityCutDef def;
        final IloRange range;
        final IloNumVar slack;

        CutRowData(CapacityCutDef def, IloRange range, IloNumVar slack) {
            this.def = def;
            this.range = range;
            this.slack = slack;
        }
    }

    private static final class DualData {
        final double[] coverDuals;
        final double vehicleLowerDual;
        final double vehicleUpperDual;
        final double[] returnArcDual;
        final double[][] internalArcDual;

        DualData(
                double[] coverDuals,
                double vehicleLowerDual,
                double vehicleUpperDual,
                double[] returnArcDual,
                double[][] internalArcDual
        ) {
            this.coverDuals = coverDuals;
            this.vehicleLowerDual = vehicleLowerDual;
            this.vehicleUpperDual = vehicleUpperDual;
            this.returnArcDual = returnArcDual;
            this.internalArcDual = internalArcDual;
        }
    }

    private static final class NodeLpResult {
        final boolean feasible;
        final boolean optimal;
        final boolean provenInfeasible;
        final String status;
        final double objective;
        final boolean integral;
        final BranchDecision branchDecision;
        final ArrayList<RouteColumn> routePool;
        final ArrayList<CapacityCutDef> capacityCuts;
        final int generatedColumns;

        NodeLpResult(
                boolean feasible,
                boolean optimal,
                boolean provenInfeasible,
                String status,
                double objective,
                boolean integral,
                BranchDecision branchDecision,
                ArrayList<RouteColumn> routePool,
                ArrayList<CapacityCutDef> capacityCuts,
                int generatedColumns
        ) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.provenInfeasible = provenInfeasible;
            this.status = status;
            this.objective = objective;
            this.integral = integral;
            this.branchDecision = branchDecision;
            this.routePool = routePool;
            this.capacityCuts = capacityCuts;
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

    private static final class StrongBranchEval {
        final boolean infeasible;
        final double objective;
        final ArrayList<RouteColumn> routePool;

        StrongBranchEval(boolean infeasible, double objective, ArrayList<RouteColumn> routePool) {
            this.infeasible = infeasible;
            this.objective = objective;
            this.routePool = routePool;
        }
    }

    private static final class NodeMaster implements AutoCloseable {
        final IloCplex cplex;
        final IloObjective obj;
        final ActiveSet active;
        final boolean[][] forbiddenArcs;
        final IloRange[] coverEq;
        final IloRange vehicleLower;
        final IloRange vehicleUpper;
        final IloNumVar[] coverArtificial;
        final IloNumVar vehicleLowerSlack;
        final ArrayList<CutRowData> cutRows;
        final ArrayList<RouteVarData> routeVars;
        final ArrayList<RouteColumn> routePool;
        final HashSet<String> existingRouteKeys;
        final int[] localIndexByGlobal;
        final double artificialPenalty;

        NodeMaster(ActiveSet active, BranchNode node, int t, long deadlineNs, double artificialPenalty) throws IloException {
            this.active = active;
            this.forbiddenArcs = node.forbiddenArcs;
            this.localIndexByGlobal = active.localIndexByGlobal;
            this.artificialPenalty = artificialPenalty;
            this.cplex = new IloCplex();
            configureLp(cplex, remainingSeconds(deadlineNs));

            int customerCount = active.size();
            this.obj = cplex.addMinimize();
            this.coverEq = new IloRange[customerCount];
            for (int local = 0; local < customerCount; local++) {
                coverEq[local] = cplex.addEq(cplex.linearNumExpr(), 1.0, "BP_Cover_" + t + "_" + local);
            }
            this.vehicleLower = cplex.addGe(cplex.linearNumExpr(), node.vehicleLowerBound, "BP_VehLB_" + t);
            this.vehicleUpper = cplex.addLe(cplex.linearNumExpr(), node.vehicleUpperBound, "BP_VehUB_" + t);

            this.coverArtificial = new IloNumVar[customerCount];
            for (int local = 0; local < customerCount; local++) {
                IloColumn col = cplex.column(obj, artificialPenalty).and(cplex.column(coverEq[local], 1.0));
                coverArtificial[local] = cplex.numVar(col, 0.0, Double.MAX_VALUE, "bp_art_" + t + "_" + local);
            }
            IloColumn vlbCol = cplex.column(obj, artificialPenalty).and(cplex.column(vehicleLower, 1.0));
            this.vehicleLowerSlack = cplex.numVar(vlbCol, 0.0, Double.MAX_VALUE, "bp_vlb_art_" + t);

            this.cutRows = new ArrayList<CutRowData>();
            this.routeVars = new ArrayList<RouteVarData>();
            this.routePool = new ArrayList<RouteColumn>();
            this.existingRouteKeys = new HashSet<String>();

            for (int idx = 0; idx < node.capacityCuts.size(); idx++) {
                addCapacityCut(node.capacityCuts.get(idx), "bp_cut_" + t + "_" + idx);
            }
        }

        boolean addRoute(RouteColumn route, String varName) throws IloException {
            if (route == null) {
                return false;
            }
            String key = route.key();
            if (existingRouteKeys.contains(key)) {
                return false;
            }
            int[] localNodePath = toLocalNodePath(route, localIndexByGlobal);
            if (localNodePath == null || !isRouteAllowed(localNodePath, forbiddenArcs)) {
                return false;
            }

            existingRouteKeys.add(key);
            int[] localCustomers = extractLocalCustomers(localNodePath);
            routePool.add(route);

            IloColumn col = cplex.column(obj, route.cost)
                    .and(cplex.column(vehicleLower, 1.0))
                    .and(cplex.column(vehicleUpper, 1.0));
            for (int idx = 0; idx < localCustomers.length; idx++) {
                col = col.and(cplex.column(coverEq[localCustomers[idx]], 1.0));
            }
            for (int cutIdx = 0; cutIdx < cutRows.size(); cutIdx++) {
                CutRowData cut = cutRows.get(cutIdx);
                int coeff = countOutgoingCrossings(localNodePath, cut.def.members);
                if (coeff > 0) {
                    col = col.and(cplex.column(cut.range, coeff));
                }
            }

            IloNumVar var = cplex.numVar(col, 0.0, Double.MAX_VALUE, "xi_bp_" + varName);
            routeVars.add(new RouteVarData(route, var, localCustomers, localNodePath));
            return true;
        }

        void addCapacityCut(CapacityCutDef cut, String name) throws IloException {
            IloRange range = cplex.addGe(cplex.linearNumExpr(), cut.rhs, name);
            for (int idx = 0; idx < routeVars.size(); idx++) {
                RouteVarData data = routeVars.get(idx);
                int coeff = countOutgoingCrossings(data.localNodePath, cut.members);
                if (coeff > 0) {
                    cplex.setLinearCoef(range, data.var, coeff);
                }
            }
            IloColumn slackCol = cplex.column(obj, artificialPenalty).and(cplex.column(range, 1.0));
            IloNumVar slack = cplex.numVar(slackCol, 0.0, Double.MAX_VALUE, name + "_art");
            cutRows.add(new CutRowData(cut, range, slack));
        }

        boolean hasPositiveCoverArtificial() throws IloException {
            for (int local = 0; local < coverArtificial.length; local++) {
                if (cplex.getValue(coverArtificial[local]) > ART_EPS) {
                    return true;
                }
            }
            return false;
        }

        boolean hasPositiveVehicleLowerSlack() throws IloException {
            return cplex.getValue(vehicleLowerSlack) > ART_EPS;
        }

        boolean hasPositiveCutSlack() throws IloException {
            for (int idx = 0; idx < cutRows.size(); idx++) {
                if (cplex.getValue(cutRows.get(idx).slack) > ART_EPS) {
                    return true;
                }
            }
            return false;
        }

        DualData extractDuals() throws IloException {
            int customerCount = active.size();
            double[] coverDuals = new double[customerCount];
            for (int local = 0; local < customerCount; local++) {
                coverDuals[local] = cplex.getDual(coverEq[local]);
            }
            double vehicleLowerDual = cplex.getDual(vehicleLower);
            double vehicleUpperDual = cplex.getDual(vehicleUpper);
            double[] returnArcDual = new double[customerCount];
            double[][] internalArcDual = new double[customerCount][customerCount];
            for (int idx = 0; idx < cutRows.size(); idx++) {
                CutRowData cut = cutRows.get(idx);
                double dual = cplex.getDual(cut.range);
                if (Math.abs(dual) <= 1e-12) {
                    continue;
                }
                for (int a = 0; a < cut.def.locals.length; a++) {
                    int from = cut.def.locals[a];
                    returnArcDual[from] += dual;
                    for (int to = 0; to < customerCount; to++) {
                        if (!cut.def.members.get(to)) {
                            internalArcDual[from][to] += dual;
                        }
                    }
                }
            }
            return new DualData(coverDuals, vehicleLowerDual, vehicleUpperDual, returnArcDual, internalArcDual);
        }

        ArrayList<ActiveRouteValue> collectPositiveRoutes() throws IloException {
            ArrayList<ActiveRouteValue> positiveRoutes = new ArrayList<ActiveRouteValue>();
            for (int idx = 0; idx < routeVars.size(); idx++) {
                RouteVarData data = routeVars.get(idx);
                double value = cplex.getValue(data.var);
                if (value > EPS) {
                    positiveRoutes.add(new ActiveRouteValue(data.route, value, data.localCustomers, data.localNodePath));
                }
            }
            return positiveRoutes;
        }

        boolean isIntegralRouteSolution() throws IloException {
            for (int idx = 0; idx < routeVars.size(); idx++) {
                double value = cplex.getValue(routeVars.get(idx).var);
                if (value > EPS && isFractional(value)) {
                    return false;
                }
            }
            return true;
        }

        ArrayList<RouteColumn> copyRoutePool() {
            return new ArrayList<RouteColumn>(routePool);
        }

        ArrayList<CapacityCutDef> copyCapacityCuts() {
            return new ArrayList<CapacityCutDef>(capacityCutsView());
        }

        ArrayList<CapacityCutDef> capacityCutsView() {
            ArrayList<CapacityCutDef> cuts = new ArrayList<CapacityCutDef>(cutRows.size());
            for (int idx = 0; idx < cutRows.size(); idx++) {
                cuts.add(cutRows.get(idx).def);
            }
            return cuts;
        }

        @Override
        public void close() {
            cplex.end();
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
        BranchNode root = BranchNode.root(active.size(), ins.K, rootSeeds);

        int exploredNodes = 0;
        int totalGeneratedColumns = 0;
        double incumbent = Double.POSITIVE_INFINITY;

        NodeLpResult rootLp = solveNodeLp(ins, t, active, root, deadlineNs);
        exploredNodes++;
        totalGeneratedColumns += rootLp.generatedColumns;

        if (!rootLp.optimal) {
            return unresolvedResult(rootLp.status, exploredNodes, totalGeneratedColumns, incumbent);
        }
        if (!rootLp.feasible) {
            return new Result(false, false, rootLp.provenInfeasible, rootLp.status, Double.NaN, exploredNodes, totalGeneratedColumns);
        }
        if (rootLp.integral) {
            return new Result(true, true, false, "OptimalBP(root)", rootLp.objective, exploredNodes, totalGeneratedColumns);
        }

        PriorityQueue<PendingNode> queue = new PriorityQueue<PendingNode>(
                16,
                new Comparator<PendingNode>() {
                    @Override
                    public int compare(PendingNode a, PendingNode b) {
                        int cmp = Double.compare(a.lpResult.objective, b.lpResult.objective);
                        if (cmp != 0) {
                            return cmp;
                        }
                        return Integer.compare(a.node.depth, b.node.depth);
                    }
                }
        );
        queue.add(new PendingNode(root, rootLp));

        while (!queue.isEmpty()) {
            if (isPastDeadline(deadlineNs)) {
                return unresolvedResult("BP_TimeLimit", exploredNodes, totalGeneratedColumns, incumbent);
            }

            PendingNode pending = queue.poll();
            if (pending.lpResult.objective >= incumbent - EPS) {
                continue;
            }

            BranchDecision decision = pending.lpResult.branchDecision;
            if (decision.type == BranchType.LEAF) {
                incumbent = Math.min(incumbent, pending.lpResult.objective);
                continue;
            }

            ArrayList<CapacityCutDef> childCuts = pending.lpResult.capacityCuts;
            ArrayList<RouteColumn> leftSeeds;
            ArrayList<RouteColumn> rightSeeds;
            BranchNode leftNode;
            BranchNode rightNode;

            if (decision.type == BranchType.VEHICLE) {
                int newUpper = (int) Math.floor(decision.vehicleValue);
                int newLower = (int) Math.ceil(decision.vehicleValue);
                leftSeeds = filterAllowedRoutes(pending.lpResult.routePool, pending.node, active.localIndexByGlobal);
                rightSeeds = filterAllowedRoutes(pending.lpResult.routePool, pending.node, active.localIndexByGlobal);
                leftNode = pending.node.childVehicleUpper(newUpper, leftSeeds, childCuts);
                rightNode = pending.node.childVehicleLower(newLower, rightSeeds, childCuts);
            } else {
                leftNode = null;
                rightNode = null;
                leftSeeds = decision.leftSeedColumns != null
                        ? decision.leftSeedColumns
                        : filterAllowedRoutesAfterForbid(
                                pending.lpResult.routePool,
                                pending.node,
                                decision.arc,
                                false,
                                active.localIndexByGlobal
                        );
                rightSeeds = decision.rightSeedColumns != null
                        ? decision.rightSeedColumns
                        : filterAllowedRoutesAfterForbid(
                                pending.lpResult.routePool,
                                pending.node,
                                decision.arc,
                                true,
                                active.localIndexByGlobal
                        );
                leftNode = pending.node.childForbidArc(decision.arc, leftSeeds, childCuts);
                rightNode = pending.node.childForceArc(decision.arc, rightSeeds, childCuts);
            }

            if (leftNode != null) {
                NodeLpResult leftLp = solveNodeLp(ins, t, active, leftNode, deadlineNs);
                exploredNodes++;
                totalGeneratedColumns += leftLp.generatedColumns;
                if (!leftLp.optimal) {
                    return unresolvedResult(leftLp.status, exploredNodes, totalGeneratedColumns, incumbent);
                }
                if (leftLp.feasible) {
                    if (leftLp.integral) {
                        incumbent = Math.min(incumbent, leftLp.objective);
                    } else if (leftLp.objective < incumbent - EPS) {
                        queue.add(new PendingNode(leftNode, leftLp));
                    }
                }
            }

            if (rightNode != null) {
                NodeLpResult rightLp = solveNodeLp(ins, t, active, rightNode, deadlineNs);
                exploredNodes++;
                totalGeneratedColumns += rightLp.generatedColumns;
                if (!rightLp.optimal) {
                    return unresolvedResult(rightLp.status, exploredNodes, totalGeneratedColumns, incumbent);
                }
                if (rightLp.feasible) {
                    if (rightLp.integral) {
                        incumbent = Math.min(incumbent, rightLp.objective);
                    } else if (rightLp.objective < incumbent - EPS) {
                        queue.add(new PendingNode(rightNode, rightLp));
                    }
                }
            }
        }

        if (isFinite(incumbent)) {
            return new Result(true, true, false, "OptimalBP", incumbent, exploredNodes, totalGeneratedColumns);
        }
        return new Result(false, false, false, "BranchPriceNoIncumbent", Double.NaN, exploredNodes, totalGeneratedColumns);
    }

    private NodeLpResult solveNodeLp(Instance ins, int t, ActiveSet active, BranchNode node, long deadlineNs) {
        if (isPastDeadline(deadlineNs)) {
            return new NodeLpResult(false, false, false, "BP_TimeLimit", Double.NaN, false, BranchDecision.leaf(),
                    node.seedColumns, node.capacityCuts, 0);
        }
        if (node.vehicleLowerBound > node.vehicleUpperBound) {
            return new NodeLpResult(false, true, true, "BP_NodeInfeasibleByVehicleBounds", Double.NaN, false, BranchDecision.leaf(),
                    node.seedColumns, node.capacityCuts, 0);
        }
        if (node.vehicleLowerBound > active.size()) {
            return new NodeLpResult(false, true, true, "BP_NodeInfeasibleByVehicleLowerBound", Double.NaN, false, BranchDecision.leaf(),
                    node.seedColumns, node.capacityCuts, 0);
        }
        int minVehiclesByLoad = (int) Math.ceil((active.totalLoad - EPS) / active.vehicleCapacity);
        if (minVehiclesByLoad > node.vehicleUpperBound) {
            return new NodeLpResult(false, true, true, "BP_NodeInfeasibleByVehicleUpperBound", Double.NaN, false, BranchDecision.leaf(),
                    node.seedColumns, node.capacityCuts, 0);
        }

        int generatedColumns = 0;
        double artificialPenalty = computeArtificialPenalty(ins, active.size());
        try (NodeMaster master = new NodeMaster(active, node, t, deadlineNs, artificialPenalty)) {
            for (int idx = 0; idx < node.seedColumns.size(); idx++) {
                master.addRoute(node.seedColumns.get(idx), "seed_" + t + "_" + idx);
            }
            for (int local = 0; local < active.size(); local++) {
                RouteColumn singleton = new RouteColumn(
                        new int[]{active.activeCustomers[local]},
                        ins.c[0][active.activeCustomers[local]] + ins.c[active.activeCustomers[local]][ins.n + 1],
                        active.qLocal[local]
                );
                master.addRoute(singleton, "single_" + t + "_" + local);
            }

            while (true) {
                if (isPastDeadline(deadlineNs)) {
                    return new NodeLpResult(false, false, false, "BP_TimeLimit", Double.NaN, false, BranchDecision.leaf(),
                            master.copyRoutePool(), master.copyCapacityCuts(), generatedColumns);
                }

                boolean solved = master.cplex.solve();
                String status = master.cplex.getStatus().toString();
                if (!solved || !status.startsWith("Optimal")) {
                    if (status.startsWith("Infeasible")) {
                        return new NodeLpResult(false, true, true, "BP_NodeInfeasibleLP", Double.NaN, false, BranchDecision.leaf(),
                                master.copyRoutePool(), master.copyCapacityCuts(), generatedColumns);
                    }
                    return new NodeLpResult(false, false, false, "BP_LP_" + status, Double.NaN, false, BranchDecision.leaf(),
                            master.copyRoutePool(), master.copyCapacityCuts(), generatedColumns);
                }

                DualData duals = master.extractDuals();
                PricingEspprcSolver.RoutePricingConstraints pricingConstraints = buildPricingConstraints(active, node, duals);
                PricingEspprcSolver.Result pricing = pricingSolver.findNegativeRoutesRelaxedThenExact(
                        ins,
                        active.activeCustomers,
                        active.qLocal,
                        duals.coverDuals,
                        duals.vehicleLowerDual + duals.vehicleUpperDual,
                        master.existingRouteKeys,
                        maxColumnsPerPricing(active.size()),
                        pricingConstraints
                );

                boolean addedColumn = false;
                if (pricing.routes != null && !pricing.routes.isEmpty()) {
                    for (int idx = 0; idx < pricing.routes.size(); idx++) {
                        if (master.addRoute(pricing.routes.get(idx), "cg_" + t + "_" + generatedColumns)) {
                            generatedColumns++;
                            addedColumn = true;
                        }
                    }
                } else if (pricing.foundNegativeColumn && pricing.route != null) {
                    if (master.addRoute(pricing.route, "cg_" + t + "_" + generatedColumns)) {
                        generatedColumns++;
                        addedColumn = true;
                    }
                }
                if (addedColumn) {
                    continue;
                }

                ArrayList<ActiveRouteValue> positiveRoutes = master.collectPositiveRoutes();
                ArrayList<CapacityCutDef> newCuts = separateCapacityCuts(active, positiveRoutes, master.capacityCutsView());
                if (!newCuts.isEmpty()) {
                    for (int idx = 0; idx < newCuts.size(); idx++) {
                        master.addCapacityCut(newCuts.get(idx), "bp_cut_" + t + "_" + master.cutRows.size() + "_" + idx);
                    }
                    continue;
                }

                if (master.hasPositiveCoverArtificial() || master.hasPositiveVehicleLowerSlack() || master.hasPositiveCutSlack()) {
                    return new NodeLpResult(false, true, true, "BP_NodeInfeasible", Double.NaN, false, BranchDecision.leaf(),
                            master.copyRoutePool(), master.copyCapacityCuts(), generatedColumns);
                }

                boolean integral = master.isIntegralRouteSolution();
                BranchDecision decision = integral
                        ? BranchDecision.leaf()
                        : chooseBranchDecision(
                                ins,
                                t,
                                active,
                                node,
                                positiveRoutes,
                                master.copyRoutePool(),
                                master.copyCapacityCuts(),
                                deadlineNs
                        );
                return new NodeLpResult(
                        true,
                        true,
                        false,
                        "BP_NodeOptimal",
                        master.cplex.getObjValue(),
                        integral,
                        decision,
                        master.copyRoutePool(),
                        master.copyCapacityCuts(),
                        generatedColumns
                );
            }
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
        return new ActiveSet(activeCustomers, qLocal, localIndexByGlobal, totalLoad, ins.Q);
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

    private static PricingEspprcSolver.RoutePricingConstraints buildPricingConstraints(
            ActiveSet active,
            BranchNode node,
            DualData duals
    ) {
        int m = active.size();
        boolean anyForbidden = false;
        boolean[] forbiddenStart = new boolean[m];
        boolean[] forbiddenReturn = new boolean[m];
        boolean[][] forbiddenInternal = new boolean[m][m];
        for (int local = 0; local < m; local++) {
            forbiddenStart[local] = node.forbiddenArcs[0][local + 1];
            forbiddenReturn[local] = node.forbiddenArcs[local + 1][0];
            anyForbidden |= forbiddenStart[local] || forbiddenReturn[local];
            for (int nxt = 0; nxt < m; nxt++) {
                forbiddenInternal[local][nxt] = node.forbiddenArcs[local + 1][nxt + 1];
                anyForbidden |= forbiddenInternal[local][nxt];
            }
        }

        boolean anyCutDual = false;
        if (duals.returnArcDual != null) {
            for (int local = 0; local < duals.returnArcDual.length; local++) {
                if (Math.abs(duals.returnArcDual[local]) > 1e-12) {
                    anyCutDual = true;
                    break;
                }
            }
        }
        if (!anyCutDual && duals.internalArcDual != null) {
            for (int i = 0; i < duals.internalArcDual.length && !anyCutDual; i++) {
                for (int j = 0; j < duals.internalArcDual[i].length; j++) {
                    if (Math.abs(duals.internalArcDual[i][j]) > 1e-12) {
                        anyCutDual = true;
                        break;
                    }
                }
            }
        }

        if (!anyForbidden && !anyCutDual) {
            return null;
        }
        return PricingEspprcSolver.RoutePricingConstraints.create(
                null,
                null,
                null,
                null,
                forbiddenStart,
                forbiddenReturn,
                forbiddenInternal,
                anyCutDual ? duals.returnArcDual : null,
                anyCutDual ? duals.internalArcDual : null
        );
    }

    private BranchDecision chooseBranchDecision(
            Instance ins,
            int t,
            ActiveSet active,
            BranchNode node,
            ArrayList<ActiveRouteValue> positiveRoutes,
            ArrayList<RouteColumn> routePool,
            ArrayList<CapacityCutDef> capacityCuts,
            long deadlineNs
    ) {
        double vehicleValue = 0.0;
        for (int idx = 0; idx < positiveRoutes.size(); idx++) {
            vehicleValue += positiveRoutes.get(idx).value;
        }
        if (isBranchFraction(vehicleValue)) {
            return BranchDecision.vehicle(vehicleValue);
        }

        BranchDecision strongArc = chooseArcByStrongBranching(
                ins,
                t,
                active,
                node,
                positiveRoutes,
                routePool,
                capacityCuts,
                deadlineNs
        );
        if (strongArc != null) {
            return strongArc;
        }
        return chooseArcHeuristically(active, positiveRoutes);
    }

    private BranchDecision chooseArcByStrongBranching(
            Instance ins,
            int t,
            ActiveSet active,
            BranchNode node,
            ArrayList<ActiveRouteValue> positiveRoutes,
            ArrayList<RouteColumn> routePool,
            ArrayList<CapacityCutDef> capacityCuts,
            long deadlineNs
    ) {
        if (active.size() > 62) {
            return null;
        }
        if (remainingSeconds(deadlineNs) < MIN_STRONG_BRANCH_REMAINING_SEC) {
            return null;
        }

        double[][] arcFlow = buildArcFlow(active.size(), positiveRoutes);
        ArrayList<ArcBranch> candidates = buildStrongBranchCandidates(arcFlow);
        if (candidates.isEmpty()) {
            return null;
        }

        ArcBranch bestArc = null;
        double bestScore = Double.NEGATIVE_INFINITY;
        ArrayList<RouteColumn> bestLeftSeeds = null;
        ArrayList<RouteColumn> bestRightSeeds = null;
        for (int idx = 0; idx < candidates.size(); idx++) {
            if (remainingSeconds(deadlineNs) < 1.0) {
                break;
            }
            ArcBranch candidate = candidates.get(idx);
            StrongBranchEval leftEval = evaluateBranchChild(
                    ins,
                    t,
                    active,
                    node.childForbidArc(
                            candidate,
                            filterAllowedRoutesAfterForbid(routePool, node, candidate, false, active.localIndexByGlobal),
                            capacityCuts
                    ),
                    deadlineNs
            );
            StrongBranchEval rightEval = evaluateBranchChild(
                    ins,
                    t,
                    active,
                    node.childForceArc(
                            candidate,
                            filterAllowedRoutesAfterForbid(routePool, node, candidate, true, active.localIndexByGlobal),
                            capacityCuts
                    ),
                    deadlineNs
            );

            if (leftEval.infeasible || rightEval.infeasible) {
                return BranchDecision.arc(candidate, leftEval.routePool, rightEval.routePool);
            }
            double score = leftEval.objective + rightEval.objective;
            if (score > bestScore + 1e-9) {
                bestScore = score;
                bestArc = candidate;
                bestLeftSeeds = leftEval.routePool;
                bestRightSeeds = rightEval.routePool;
            }
        }
        return bestArc == null ? null : BranchDecision.arc(bestArc, bestLeftSeeds, bestRightSeeds);
    }

    private static BranchDecision chooseArcHeuristically(ActiveSet active, ArrayList<ActiveRouteValue> positiveRoutes) {
        double[][] arcFlow = buildArcFlow(active.size(), positiveRoutes);
        ArcBranch bestArc = null;
        double bestScore = Double.NEGATIVE_INFINITY;
        for (int from = 0; from < arcFlow.length; from++) {
            for (int to = 0; to < arcFlow.length; to++) {
                if (from == to || (from == 0 && to == 0)) {
                    continue;
                }
                double value = arcFlow[from][to];
                if (!isBranchFraction(value)) {
                    continue;
                }
                double score = 0.5 - Math.abs(0.5 - value);
                if (from != 0 && to != 0) {
                    score += 1e-6;
                }
                if (score > bestScore + 1e-12) {
                    bestScore = score;
                    bestArc = new ArcBranch(from, to, value);
                }
            }
        }
        if (bestArc != null) {
            return BranchDecision.arc(bestArc);
        }

        ActiveRouteValue fallbackRoute = null;
        double fallbackScore = Double.POSITIVE_INFINITY;
        for (int idx = 0; idx < positiveRoutes.size(); idx++) {
            ActiveRouteValue route = positiveRoutes.get(idx);
            if (!isFractional(route.value)) {
                continue;
            }
            double score = Math.abs(0.5 - route.value);
            if (score < fallbackScore - 1e-12) {
                fallbackScore = score;
                fallbackRoute = route;
            }
        }
        if (fallbackRoute != null) {
            for (int p = 0; p < fallbackRoute.localNodePath.length - 1; p++) {
                int from = fallbackRoute.localNodePath[p];
                int to = fallbackRoute.localNodePath[p + 1];
                if (from != to) {
                    return BranchDecision.arc(new ArcBranch(from, to, arcFlow[from][to]));
                }
            }
        }
        return BranchDecision.leaf();
    }

    private static ArrayList<ArcBranch> buildStrongBranchCandidates(double[][] arcFlow) {
        ArrayList<ArcBranch> candidates = new ArrayList<ArcBranch>();
        for (int from = 0; from < arcFlow.length; from++) {
            for (int to = 0; to < arcFlow.length; to++) {
                if (from == to || (from == 0 && to == 0)) {
                    continue;
                }
                double value = arcFlow[from][to];
                if (!isBranchFraction(value)) {
                    continue;
                }
                candidates.add(new ArcBranch(from, to, value));
            }
        }
        Collections.sort(candidates, new Comparator<ArcBranch>() {
            @Override
            public int compare(ArcBranch a, ArcBranch b) {
                double sa = Math.abs(0.5 - a.value);
                double sb = Math.abs(0.5 - b.value);
                int cmp = Double.compare(sa, sb);
                if (cmp != 0) {
                    return cmp;
                }
                return Double.compare(b.value, a.value);
            }
        });
        if (candidates.size() > MAX_STRONG_BRANCH_CANDIDATES) {
            return new ArrayList<ArcBranch>(candidates.subList(0, MAX_STRONG_BRANCH_CANDIDATES));
        }
        return candidates;
    }

    private StrongBranchEval evaluateBranchChild(
            Instance ins,
            int t,
            ActiveSet active,
            BranchNode node,
            long deadlineNs
    ) {
        if (node.vehicleLowerBound > node.vehicleUpperBound || node.vehicleLowerBound > active.size()) {
            return new StrongBranchEval(true, Double.POSITIVE_INFINITY, node.seedColumns);
        }
        int minVehiclesByLoad = (int) Math.ceil((active.totalLoad - EPS) / active.vehicleCapacity);
        if (minVehiclesByLoad > node.vehicleUpperBound) {
            return new StrongBranchEval(true, Double.POSITIVE_INFINITY, node.seedColumns);
        }

        double artificialPenalty = computeArtificialPenalty(ins, active.size());
        try (NodeMaster master = new NodeMaster(active, node, t, deadlineNs, artificialPenalty)) {
            for (int idx = 0; idx < node.seedColumns.size(); idx++) {
                master.addRoute(node.seedColumns.get(idx), "sb_seed_" + t + "_" + idx);
            }
            for (int local = 0; local < active.size(); local++) {
                RouteColumn singleton = new RouteColumn(
                        new int[]{active.activeCustomers[local]},
                        ins.c[0][active.activeCustomers[local]] + ins.c[active.activeCustomers[local]][ins.n + 1],
                        active.qLocal[local]
                );
                master.addRoute(singleton, "sb_single_" + t + "_" + local);
            }

            while (remainingSeconds(deadlineNs) > 1.0) {
                boolean solved = master.cplex.solve();
                String status = master.cplex.getStatus().toString();
                if (!solved || !status.startsWith("Optimal")) {
                    return new StrongBranchEval(true, Double.POSITIVE_INFINITY, master.copyRoutePool());
                }

                DualData duals = master.extractDuals();
                PricingEspprcSolver.RoutePricingConstraints pricingConstraints = buildPricingConstraints(active, node, duals);
                PricingEspprcSolver.Result pricing = pricingSolver.findNegativeRoutesRelaxedOnly(
                        ins,
                        active.activeCustomers,
                        active.qLocal,
                        duals.coverDuals,
                        duals.vehicleLowerDual + duals.vehicleUpperDual,
                        master.existingRouteKeys,
                        Math.max(4, Math.min(12, maxColumnsPerPricing(active.size()))),
                        pricingConstraints
                );
                boolean addedColumn = false;
                if (pricing.routes != null && !pricing.routes.isEmpty()) {
                    for (int idx = 0; idx < pricing.routes.size(); idx++) {
                        if (master.addRoute(pricing.routes.get(idx), "sb_cg_" + t + "_" + idx)) {
                            addedColumn = true;
                        }
                    }
                } else if (pricing.foundNegativeColumn && pricing.route != null) {
                    addedColumn = master.addRoute(pricing.route, "sb_cg_" + t);
                }
                if (addedColumn) {
                    continue;
                }

                if (master.hasPositiveCoverArtificial() || master.hasPositiveVehicleLowerSlack() || master.hasPositiveCutSlack()) {
                    return new StrongBranchEval(true, Double.POSITIVE_INFINITY, master.copyRoutePool());
                }
                return new StrongBranchEval(false, master.cplex.getObjValue(), master.copyRoutePool());
            }
        } catch (IloException e) {
            return new StrongBranchEval(true, Double.POSITIVE_INFINITY, node.seedColumns);
        }
        return new StrongBranchEval(true, Double.POSITIVE_INFINITY, node.seedColumns);
    }

    private static ArrayList<RouteColumn> filterAllowedRoutes(
            ArrayList<RouteColumn> routes,
            BranchNode node,
            int[] localIndexByGlobal
    ) {
        ArrayList<RouteColumn> filtered = new ArrayList<RouteColumn>();
        HashSet<String> seen = new HashSet<String>();
        for (int idx = 0; idx < routes.size(); idx++) {
            RouteColumn route = routes.get(idx);
            int[] path = toLocalNodePath(route, localIndexByGlobal);
            if (path != null && isRouteAllowed(path, node.forbiddenArcs) && seen.add(route.key())) {
                filtered.add(route);
            }
        }
        return filtered;
    }

    private static ArrayList<RouteColumn> filterAllowedRoutesAfterForbid(
            ArrayList<RouteColumn> routes,
            BranchNode node,
            ArcBranch arc,
            boolean forceArc,
            int[] localIndexByGlobal
    ) {
        BranchNode view = forceArc
                ? node.childForceArc(arc, new ArrayList<RouteColumn>(), node.capacityCuts)
                : node.childForbidArc(arc, new ArrayList<RouteColumn>(), node.capacityCuts);
        return filterAllowedRoutes(routes, view, localIndexByGlobal);
    }

    private static int[] toLocalNodePath(RouteColumn route, int[] localIndexByGlobal) {
        if (route == null) {
            return null;
        }
        int len = route.globalCustomersInOrder.length;
        int[] path = new int[len + 2];
        path[0] = 0;
        for (int pos = 0; pos < len; pos++) {
            int global = route.globalCustomersInOrder[pos];
            if (global < 0 || global >= localIndexByGlobal.length || localIndexByGlobal[global] < 0) {
                return null;
            }
            path[pos + 1] = localIndexByGlobal[global] + 1;
        }
        path[len + 1] = 0;
        return path;
    }

    private static int[] extractLocalCustomers(int[] localNodePath) {
        int[] locals = new int[localNodePath.length - 2];
        for (int pos = 1; pos < localNodePath.length - 1; pos++) {
            locals[pos - 1] = localNodePath[pos] - 1;
        }
        return locals;
    }

    private static boolean isRouteAllowed(int[] localNodePath, boolean[][] forbiddenArcs) {
        for (int pos = 0; pos < localNodePath.length - 1; pos++) {
            if (forbiddenArcs[localNodePath[pos]][localNodePath[pos + 1]]) {
                return false;
            }
        }
        return true;
    }

    private static int countOutgoingCrossings(int[] localNodePath, BitSet members) {
        int count = 0;
        for (int pos = 0; pos < localNodePath.length - 1; pos++) {
            int fromNode = localNodePath[pos];
            if (fromNode == 0) {
                continue;
            }
            int fromLocal = fromNode - 1;
            if (!members.get(fromLocal)) {
                continue;
            }
            int toNode = localNodePath[pos + 1];
            if (toNode == 0) {
                count++;
            } else if (!members.get(toNode - 1)) {
                count++;
            }
        }
        return count;
    }

    private static double[][] buildArcFlow(int customerCount, ArrayList<ActiveRouteValue> positiveRoutes) {
        double[][] arcFlow = new double[customerCount + 1][customerCount + 1];
        for (int idx = 0; idx < positiveRoutes.size(); idx++) {
            ActiveRouteValue route = positiveRoutes.get(idx);
            for (int pos = 0; pos < route.localNodePath.length - 1; pos++) {
                arcFlow[route.localNodePath[pos]][route.localNodePath[pos + 1]] += route.value;
            }
        }
        return arcFlow;
    }

    private static ArrayList<CapacityCutDef> separateCapacityCuts(
            ActiveSet active,
            ArrayList<ActiveRouteValue> positiveRoutes,
            ArrayList<CapacityCutDef> existingCuts
    ) {
        ArrayList<CapacityCutDef> cuts = new ArrayList<CapacityCutDef>();
        if (active.size() <= 2 || positiveRoutes.isEmpty()) {
            return cuts;
        }

        double[][] arcFlow = buildArcFlow(active.size(), positiveRoutes);
        HashSet<String> seen = new HashSet<String>();
        for (int idx = 0; idx < existingCuts.size(); idx++) {
            seen.add(existingCuts.get(idx).signature);
        }

        shrinkSeparate(active, arcFlow, seen, cuts);
        if (cuts.isEmpty()) {
            connectedSeparate(active, arcFlow, seen, cuts);
        }
        return cuts;
    }

    private static void shrinkSeparate(
            ActiveSet active,
            double[][] arcFlow,
            HashSet<String> seen,
            ArrayList<CapacityCutDef> out
    ) {
        ArrayList<BitSet> clusters = new ArrayList<BitSet>();
        for (int local = 0; local < active.size(); local++) {
            BitSet cluster = new BitSet(active.size());
            cluster.set(local);
            clusters.add(cluster);
        }

        while (clusters.size() > 1 && out.size() < MAX_NEW_CUTS_PER_NODE) {
            double bestWeight = 0.0;
            int bestI = -1;
            int bestJ = -1;
            for (int i = 0; i < clusters.size(); i++) {
                for (int j = i + 1; j < clusters.size(); j++) {
                    double weight = symmetricFlowBetween(clusters.get(i), clusters.get(j), arcFlow);
                    if (weight > bestWeight + 1e-12) {
                        bestWeight = weight;
                        bestI = i;
                        bestJ = j;
                    }
                }
            }
            if (bestI < 0 || bestWeight <= 0.5 + 1e-12) {
                break;
            }

            BitSet merged = (BitSet) clusters.get(bestI).clone();
            merged.or(clusters.get(bestJ));
            if (merged.cardinality() >= 2 && merged.cardinality() < active.size()) {
                maybeAddCapacityCut(active, merged, arcFlow, seen, out);
            }
            clusters.remove(bestJ);
            clusters.set(bestI, merged);
        }
    }

    private static void connectedSeparate(
            ActiveSet active,
            double[][] arcFlow,
            HashSet<String> seen,
            ArrayList<CapacityCutDef> out
    ) {
        BitSet[] neighbors = new BitSet[active.size()];
        for (int i = 0; i < active.size(); i++) {
            neighbors[i] = new BitSet(active.size());
            for (int j = 0; j < active.size(); j++) {
                if (i == j) {
                    continue;
                }
                double support = arcFlow[i + 1][j + 1] + arcFlow[j + 1][i + 1];
                if (support > ARC_SUPPORT_EPS) {
                    neighbors[i].set(j);
                }
            }
        }

        HashSet<String> visitedStates = new HashSet<String>();
        for (int root = 0; root < active.size() && out.size() < MAX_NEW_CUTS_PER_NODE; root++) {
            BitSet current = new BitSet(active.size());
            current.set(root);
            BitSet frontier = (BitSet) neighbors[root].clone();
            dfsConnectedCuts(active, arcFlow, neighbors, current, frontier, seen, visitedStates, out);
        }
    }

    private static void dfsConnectedCuts(
            ActiveSet active,
            double[][] arcFlow,
            BitSet[] neighbors,
            BitSet current,
            BitSet frontier,
            HashSet<String> seenCuts,
            HashSet<String> visitedStates,
            ArrayList<CapacityCutDef> out
    ) {
        if (out.size() >= MAX_NEW_CUTS_PER_NODE || current.cardinality() > MAX_CONNECTED_CUT_SIZE) {
            return;
        }
        String stateKey = bitSetSignature(current);
        if (!visitedStates.add(stateKey)) {
            return;
        }

        if (current.cardinality() >= 2 && current.cardinality() < active.size()) {
            double lhs = outgoingFlow(current, arcFlow, active.size());
            double rhs = cutRhs(current, active.qLocal, active.vehicleCapacity);
            if (lhs + CUT_VIOL_EPS < rhs) {
                maybeAddCapacityCut(active, current, arcFlow, seenCuts, out);
                if (out.size() >= MAX_NEW_CUTS_PER_NODE) {
                    return;
                }
            }
            if (lhs - rhs > 1.2) {
                return;
            }
        }

        int next = frontier.nextSetBit(0);
        while (next >= 0 && out.size() < MAX_NEW_CUTS_PER_NODE) {
            BitSet nextCurrent = (BitSet) current.clone();
            nextCurrent.set(next);
            BitSet nextFrontier = (BitSet) frontier.clone();
            nextFrontier.or(neighbors[next]);
            nextFrontier.andNot(nextCurrent);
            dfsConnectedCuts(active, arcFlow, neighbors, nextCurrent, nextFrontier, seenCuts, visitedStates, out);
            next = frontier.nextSetBit(next + 1);
        }
    }

    private static void maybeAddCapacityCut(
            ActiveSet active,
            BitSet members,
            double[][] arcFlow,
            HashSet<String> seen,
            ArrayList<CapacityCutDef> out
    ) {
        String signature = bitSetSignature(members);
        if (seen.contains(signature)) {
            return;
        }
        double lhs = outgoingFlow(members, arcFlow, active.size());
        double rhs = cutRhs(members, active.qLocal, active.vehicleCapacity);
        if (lhs + CUT_VIOL_EPS >= rhs) {
            return;
        }

        int[] locals = new int[members.cardinality()];
        int idx = 0;
        int bit = members.nextSetBit(0);
        while (bit >= 0) {
            locals[idx++] = bit;
            bit = members.nextSetBit(bit + 1);
        }
        out.add(new CapacityCutDef((BitSet) members.clone(), locals, rhs, signature));
        seen.add(signature);
    }

    private static double symmetricFlowBetween(BitSet a, BitSet b, double[][] arcFlow) {
        double weight = 0.0;
        int i = a.nextSetBit(0);
        while (i >= 0) {
            int j = b.nextSetBit(0);
            while (j >= 0) {
                weight += arcFlow[i + 1][j + 1] + arcFlow[j + 1][i + 1];
                j = b.nextSetBit(j + 1);
            }
            i = a.nextSetBit(i + 1);
        }
        return weight;
    }

    private static double outgoingFlow(BitSet members, double[][] arcFlow, int customerCount) {
        double lhs = 0.0;
        int from = members.nextSetBit(0);
        while (from >= 0) {
            lhs += arcFlow[from + 1][0];
            for (int to = 0; to < customerCount; to++) {
                if (!members.get(to)) {
                    lhs += arcFlow[from + 1][to + 1];
                }
            }
            from = members.nextSetBit(from + 1);
        }
        return lhs;
    }

    private static double cutRhs(BitSet members, double[] qLocal, double capacity) {
        double demand = 0.0;
        int bit = members.nextSetBit(0);
        while (bit >= 0) {
            demand += qLocal[bit];
            bit = members.nextSetBit(bit + 1);
        }
        return Math.ceil((demand - 1e-9) / capacity);
    }

    private static String bitSetSignature(BitSet set) {
        StringBuilder sb = new StringBuilder();
        int bit = set.nextSetBit(0);
        boolean first = true;
        while (bit >= 0) {
            if (!first) {
                sb.append(',');
            }
            sb.append(bit);
            first = false;
            bit = set.nextSetBit(bit + 1);
        }
        return sb.toString();
    }

    private static int maxColumnsPerPricing(int customerCount) {
        return 120;
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

    private static boolean isFractional(double value) {
        double floor = Math.floor(value);
        return value > floor + EPS && value < floor + 1.0 - EPS;
    }

    private static boolean isBranchFraction(double value) {
        if (!isFractional(value)) {
            return false;
        }
        return fracCost(value) <= FRAC_DISTANCE_LIMIT + 1e-12;
    }

    private static double fracCost(double value) {
        double line = Math.floor(value) + 0.5;
        return Math.abs(value - line);
    }

    private static Result unresolvedResult(String status, int exploredNodes, int totalGeneratedColumns, double incumbent) {
        if (isFinite(incumbent)) {
            return new Result(true, false, false, status + "WithIncumbent", incumbent, exploredNodes, totalGeneratedColumns);
        }
        return new Result(false, false, false, status, Double.NaN, exploredNodes, totalGeneratedColumns);
    }

    private static boolean isFinite(double value) {
        return !Double.isNaN(value) && !Double.isInfinite(value);
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
