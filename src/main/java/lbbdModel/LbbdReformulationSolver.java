package lbbdModel;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import instance.Instance;
import lbbdModel.rmp.PeriodRouteMasterLpSolver;
import lbbdModel.rmp.RoutingLowerBoundSolver;
import lbbdModel.rmp.T1InitialCutSolver;
import lbbdModel.rmp.T2InitialCutSolver;
import model.CplexConfig;
import model.SolveResult;

import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public final class LbbdReformulationSolver {

    private static final double CUT_EPS = 1e-6;
    private static final double RESULT_TOL = 1e-4;
    private static final double INCUMBENT_PRUNE_TOL = 1e-4;
    private static final int MAX_ITERATIONS = 5000;
    private static final int MAX_RMP_PERIODS_PER_ITER = 3;
    private static final boolean ENABLE_INCUMBENT_PRUNE_NOGOOD = false;
    private static final boolean ENABLE_GLOBAL_COMBINATION_OPT_CUT = true;
    private static final boolean ENABLE_ITER_TIMING_LOG = true;
    private static final boolean LOG_TO_CONSOLE = false;
    private static final boolean ENABLE_HINT_SHORTCUT = false;
    private static final int MASTER_THREADS = 0;

    // Adaptive parameters based on instance size
    private double masterMipGapEarly;
    private double masterMipGapMid;
    private double masterMipGapLate;
    private double masterMipGapVeryLate;
    private int masterMipGapEarlyIterCutoff;
    private int masterMipGapMidIterCutoff;
    private int masterMipGapVeryLateIterCutoff;

    private int activeRmpDualCutsCapBase;
    private int activeRmpDualCutsCapLate;
    private int activePeriodOptCutsCapBase;
    private int activePeriodOptCutsCapLate;
    private int activeGlobalOptCutsCapBase;
    private int activeGlobalOptCutsCapLate;
    private int activeIncPruneCutsCapBase;
    private int activeIncPruneCutsCapLate;
    private int activeCutCapLateStartIter;

    private int maxPeriodOptCutsPerIter;
    private int latePeriodOptCutsPerIter;
    private int latePeriodOptStartIter;
    private double latePeriodOptMinDelta;
    private int globalOptCutStartIter;
    private double globalOptCutMinDelta;
    private int incumbentPruneStartIter;
    private double incumbentPruneMinExcess;

    private final Instance ins;

    public LbbdReformulationSolver(Instance ins) {
        this.ins = ins;

        // Adaptive parameter configuration based on instance size
        if (ins.n <= 12) {
            // Small scale: preserve original parameters to avoid performance degradation
            initSmallScaleParams();
        } else if (ins.n <= 20) {
            // Medium scale: moderate optimization
            initMediumScaleParams();
        } else {
            // Large scale: aggressive optimization
            initLargeScaleParams();
        }
    }

    private void initSmallScaleParams() {
        // MIP gap: keep relaxed to avoid over-solving
        this.masterMipGapEarly = 3e-3;
        this.masterMipGapMid = 7e-4;
        this.masterMipGapLate = 2.5e-4;
        this.masterMipGapVeryLate = 1.0e-4;
        this.masterMipGapEarlyIterCutoff = 100;
        this.masterMipGapMidIterCutoff = 280;
        this.masterMipGapVeryLateIterCutoff = 360;

        // Cut capacity: preserve original values to avoid oversized master problem
        this.activeRmpDualCutsCapBase = 600;
        this.activeRmpDualCutsCapLate = 1000;
        this.activePeriodOptCutsCapBase = 120;
        this.activePeriodOptCutsCapLate = 80;
        this.activeGlobalOptCutsCapBase = 90;
        this.activeGlobalOptCutsCapLate = 60;
        this.activeIncPruneCutsCapBase = 120;
        this.activeIncPruneCutsCapLate = 90;
        this.activeCutCapLateStartIter = 420;

        // Cut generation: conservative strategy
        this.maxPeriodOptCutsPerIter = 1;
        this.latePeriodOptCutsPerIter = 2;
        this.latePeriodOptStartIter = 220;
        this.latePeriodOptMinDelta = 2500.0;
        this.globalOptCutStartIter = 140;
        this.globalOptCutMinDelta = 2000.0;
        this.incumbentPruneStartIter = 180;
        this.incumbentPruneMinExcess = 200.0;
    }

    private void initMediumScaleParams() {
        // MIP gap: moderately tightened
        this.masterMipGapEarly = 2e-3;
        this.masterMipGapMid = 6e-4;
        this.masterMipGapLate = 2.2e-4;
        this.masterMipGapVeryLate = 1.0e-4;
        this.masterMipGapEarlyIterCutoff = 90;
        this.masterMipGapMidIterCutoff = 240;
        this.masterMipGapVeryLateIterCutoff = 330;

        // Cut capacity: moderately increased
        this.activeRmpDualCutsCapBase = 700;
        this.activeRmpDualCutsCapLate = 1100;
        this.activePeriodOptCutsCapBase = 135;
        this.activePeriodOptCutsCapLate = 90;
        this.activeGlobalOptCutsCapBase = 105;
        this.activeGlobalOptCutsCapLate = 70;
        this.activeIncPruneCutsCapBase = 135;
        this.activeIncPruneCutsCapLate = 100;
        this.activeCutCapLateStartIter = 360;

        // Cut generation: moderately aggressive
        this.maxPeriodOptCutsPerIter = 1;
        this.latePeriodOptCutsPerIter = 2;
        this.latePeriodOptStartIter = 185;
        this.latePeriodOptMinDelta = 2000.0;
        this.globalOptCutStartIter = 120;
        this.globalOptCutMinDelta = 1600.0;
        this.incumbentPruneStartIter = 150;
        this.incumbentPruneMinExcess = 175.0;
    }

    private void initLargeScaleParams() {
        // MIP gap: aggressively tightened
        this.masterMipGapEarly = 1.5e-3;
        this.masterMipGapMid = 5e-4;
        this.masterMipGapLate = 2e-4;
        this.masterMipGapVeryLate = 1.0e-4;
        this.masterMipGapEarlyIterCutoff = 80;
        this.masterMipGapMidIterCutoff = 200;
        this.masterMipGapVeryLateIterCutoff = 300;

        // Cut capacity: significantly increased
        this.activeRmpDualCutsCapBase = 800;
        this.activeRmpDualCutsCapLate = 1200;
        this.activePeriodOptCutsCapBase = 150;
        this.activePeriodOptCutsCapLate = 100;
        this.activeGlobalOptCutsCapBase = 120;
        this.activeGlobalOptCutsCapLate = 80;
        this.activeIncPruneCutsCapBase = 150;
        this.activeIncPruneCutsCapLate = 110;
        this.activeCutCapLateStartIter = 300;

        // Cut generation: aggressive strategy
        this.maxPeriodOptCutsPerIter = 2;
        this.latePeriodOptCutsPerIter = 3;
        this.latePeriodOptStartIter = 150;
        this.latePeriodOptMinDelta = 1500.0;
        this.globalOptCutStartIter = 100;
        this.globalOptCutMinDelta = 1200.0;
        this.incumbentPruneStartIter = 120;
        this.incumbentPruneMinExcess = 150.0;
    }

    private static final class ActiveRmpDualCut {
        final String key;
        final IloRange range;

        ActiveRmpDualCut(String key, IloRange range) {
            this.key = key;
            this.range = range;
        }
    }

    private static final class ActivePeriodOptCut {
        final String key;
        final IloRange range;

        ActivePeriodOptCut(String key, IloRange range) {
            this.key = key;
            this.range = range;
        }
    }

    private static final class ActiveGlobalOptCut {
        final String key;
        final IloRange range;

        ActiveGlobalOptCut(String key, IloRange range) {
            this.key = key;
            this.range = range;
        }
    }

    private static final class ActiveIncumbentPruneCut {
        final String key;
        final IloRange range;

        ActiveIncumbentPruneCut(String key, IloRange range) {
            this.key = key;
            this.range = range;
        }
    }

    public SolveResult solve(double targetObjective, double targetTol) {
        return solve(targetObjective, targetTol, null);
    }

    public SolveResult solve(double targetObjective, double targetTol, int[][] prevVisitHint) {
        long startNs = System.nanoTime();
        if (ENABLE_HINT_SHORTCUT && !Double.isNaN(targetObjective) && prevVisitHint != null) {
            SolveResult fast = tryValidateByHint(ins, targetObjective, targetTol, prevVisitHint, startNs);
            if (fast != null) {
                return fast;
            }
        }

        HashMap<String, PeriodCvrpSubproblemSolver.Result> subproblemCache =
                new HashMap<String, PeriodCvrpSubproblemSolver.Result>();
        HashMap<String, PeriodRouteMasterLpSolver.Result> rmpCutCache =
                new HashMap<String, PeriodRouteMasterLpSolver.Result>();
        HashSet<String> addedRmpCutKeys = new HashSet<String>();
        ArrayDeque<ActiveRmpDualCut> activeRmpDualCuts = new ArrayDeque<ActiveRmpDualCut>();
        HashSet<String> addedPeriodOptCutKeys = new HashSet<String>();
        ArrayDeque<ActivePeriodOptCut> activePeriodOptCuts = new ArrayDeque<ActivePeriodOptCut>();
        HashSet<String> addedGlobalOptCutKeys = new HashSet<String>();
        ArrayDeque<ActiveGlobalOptCut> activeGlobalOptCuts = new ArrayDeque<ActiveGlobalOptCut>();
        HashSet<String> activeIncumbentPruneCutKeys = new HashSet<String>();
        ArrayDeque<ActiveIncumbentPruneCut> activeIncumbentPruneCuts = new ArrayDeque<ActiveIncumbentPruneCut>();

        int subThreads = Math.min(ins.l, Runtime.getRuntime().availableProcessors());
        final PeriodCvrpSubproblemSolver[] subSolvers = new PeriodCvrpSubproblemSolver[ins.l + 1];
        for (int tt = 1; tt <= ins.l; tt++) {
            subSolvers[tt] = new PeriodCvrpSubproblemSolver();
        }
        ExecutorService subExecutor = subThreads > 1 ? Executors.newFixedThreadPool(subThreads) : null;

        try (IloCplex cplex = new IloCplex()) {
            configure(cplex);
            MasterModel master = new MasterModel(cplex, ins);
            PeriodRouteMasterLpSolver rmpLpSolver = new PeriodRouteMasterLpSolver(ins);
            RoutingLowerBoundSolver lbRSolver = new RoutingLowerBoundSolver(ins);
            RoutingLowerBoundSolver.Result lbRResult = lbRSolver.solve();
            double lbR = (lbRResult.feasible && !Double.isNaN(lbRResult.lbR) && !Double.isInfinite(lbRResult.lbR))
                    ? Math.max(0.0, lbRResult.lbR)
                    : 0.0;
            System.out.println("[LBBD] LB_R="
                    + fmt(lbR)
                    + " status=" + lbRResult.status
                    + " iters=" + lbRResult.iterations
                    + " cols=" + lbRResult.generatedColumns);

            // === T1 Initial Cuts (speed.tex Section B.2, eq 4-9), parallel across periods ===
            int t1CutsAdded = 0;
            T1InitialCutSolver.T1Result[] t1Results = new T1InitialCutSolver.T1Result[ins.l + 1];
            int t1Threads = Math.min(ins.l, Math.min(Runtime.getRuntime().availableProcessors(), 4));
            ExecutorService t1Executor = t1Threads > 1 ? Executors.newFixedThreadPool(t1Threads) : null;
            @SuppressWarnings("unchecked")
            Future<T1InitialCutSolver.T1Result>[] t1Futures =
                    (t1Executor != null) ? new Future[ins.l + 1] : null;
            if (t1Executor != null) {
                try {
                    for (int t = 1; t <= ins.l; t++) {
                        final int period = t;
                        t1Futures[t] = t1Executor.submit(new Callable<T1InitialCutSolver.T1Result>() {
                            @Override
                            public T1InitialCutSolver.T1Result call() {
                                return new T1InitialCutSolver(ins).solve(period);
                            }
                        });
                    }
                    for (int t = 1; t <= ins.l; t++) {
                        try {
                            t1Results[t] = t1Futures[t].get();
                        } catch (Exception e) {
                            System.err.println("[LBBD] T1 parallel solve failed for t=" + t + ": " + e.getMessage());
                        }
                    }
                } finally {
                    t1Executor.shutdownNow();
                }
                T1InitialCutSolver fallbackSolver = new T1InitialCutSolver(ins);
                for (int t = 1; t <= ins.l; t++) {
                    if (t1Results[t] == null) {
                        t1Results[t] = fallbackSolver.solve(t);
                    }
                }
            } else {
                T1InitialCutSolver t1Solver = new T1InitialCutSolver(ins);
                for (int t = 1; t <= ins.l; t++) {
                    t1Results[t] = t1Solver.solve(t);
                }
            }
            for (int t = 1; t <= ins.l; t++) {
                T1InitialCutSolver.T1Result t1 = t1Results[t];
                if (t1 != null
                        && t1.feasible
                        && t1.optimal
                        && t1.pricingProvedOptimal
                        && t1.artificialClean
                        && t1.dualU0 <= 1e-7
                        && hasMeaningfulDualCoeffs(t, t1.dualW, t1.dualU0)) {
                    master.addDualInitialCut(t, t1.dualW, t1.dualU0, "T1_InitCut_t" + t);
                    t1CutsAdded++;
                }
            }
            System.out.println("[LBBD] T1 initial cuts added: " + t1CutsAdded + "/" + ins.l);

            // === T2 Initial Cuts (speed.tex Section B.3, eq 4-12 in per-period form) ===
            T2InitialCutSolver t2Solver = new T2InitialCutSolver(ins);
            T2InitialCutSolver.T2Result t2 = t2Solver.solve();
            int t2CutsAdded = 0;
            if (t2.feasible && t2.optimal && t2.pricingProvedOptimal && t2.artificialClean
                    && t2.dualWByPeriod != null && t2.dualU0ByPeriod != null) {
                for (int t = 1; t <= ins.l; t++) {
                    if (t2.dualU0ByPeriod[t] <= 1e-7
                            && hasMeaningfulDualCoeffs(t, t2.dualWByPeriod[t], t2.dualU0ByPeriod[t])) {
                        master.addDualInitialCut(t, t2.dualWByPeriod[t], t2.dualU0ByPeriod[t], "T2_InitCut_t" + t);
                        t2CutsAdded++;
                    }
                }
            }
            System.out.println("[LBBD] T2 initial cuts added: " + t2CutsAdded + "/" + ins.l);

            int iteration = 0;
            int feasibilityCuts = 0;
            int optimalityCuts = 0;
            int lpDualCuts = 0;
            int lpDualCutSkips = 0;
            int targetPruneCuts = 0;
            int incumbentPruneCuts = 0;
            int globalOptCuts = 0;
            long totalMasterSolveNs = 0L;
            long totalSubSolveNs = 0L;
            long totalRmpSolveNs = 0L;
            long totalSubCalls = 0L;
            long totalSubCacheMisses = 0L;
            long totalRmpCalls = 0L;
            long totalRmpCacheMisses = 0L;

            double bestUpperBound = Double.POSITIVE_INFINITY;
            double bestUpperBoundGapRef = Double.NaN;
            double bestLowerBound = Double.NEGATIVE_INFINITY;

            boolean allMasterOptimal = true;
            boolean allSubproblemsOptimal = true;
            boolean forceStrictMasterGap = false;
            MasterPoint prevPoint = null;

            while (iteration < MAX_ITERATIONS) {
                iteration++;
                trimActiveRmpDualCuts(master, activeRmpDualCuts, addedRmpCutKeys, activeRmpDualCutCap(iteration));
                trimActivePeriodOptCuts(master, activePeriodOptCuts, addedPeriodOptCutKeys, activePeriodOptCutCap(iteration));
                trimActiveGlobalOptCuts(master, activeGlobalOptCuts, addedGlobalOptCutKeys, activeGlobalOptCutCap(iteration));
                trimActiveIncPruneCuts(master, activeIncumbentPruneCuts, activeIncumbentPruneCutKeys, activeIncPruneCutCap(iteration));

                long iterMasterSolveNs = 0L;
                long iterSubSolveNs = 0L;
                long iterRmpSolveNs = 0L;
                int iterSubCalls = 0;
                int iterSubCacheMisses = 0;
                int iterRmpCalls = 0;
                int iterRmpCacheMisses = 0;

                if (isFinite(bestUpperBound)) {
                    setMasterUpperCutoff(cplex, bestUpperBound - INCUMBENT_PRUNE_TOL);
                }
                double masterMipGapThisIter = recommendedMasterMipGap(iteration, forceStrictMasterGap);
                cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, masterMipGapThisIter);

                if (prevPoint != null) {
                    master.setMIPStart(prevPoint);
                }

                long masterStartNs = System.nanoTime();
                boolean solved = cplex.solve();
                iterMasterSolveNs = System.nanoTime() - masterStartNs;
                totalMasterSolveNs += iterMasterSolveNs;
                String masterStatus = safeStatus(cplex);
                if (!solved) {
                    if (isFinite(bestUpperBound) && isInfeasibleLikeStatus(masterStatus)) {
                        return buildCutoffExhaustedResult(
                                startNs,
                                iteration,
                                feasibilityCuts,
                                optimalityCuts,
                                lpDualCuts,
                                allSubproblemsOptimal,
                                bestUpperBound,
                                bestLowerBound
                        );
                    }
                    return buildFailedResult(startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts, masterStatus);
                }
                if (!masterStatus.startsWith("Optimal")) {
                    allMasterOptimal = false;
                }

                MasterPoint point = master.extractPoint(ins);
                prevPoint = point;
                bestLowerBound = Math.max(bestLowerBound, point.bmpBestBound);

                double phi = 0.0;
                double[] phiByPeriod = new double[ins.l + 1];
                String[] subKeysByPeriod = new String[ins.l + 1];
                int infeasiblePeriod = -1;

                // Phase 1: check cache and submit parallel tasks for misses
                PeriodCvrpSubproblemSolver.Result[] subResults = new PeriodCvrpSubproblemSolver.Result[ins.l + 1];
                @SuppressWarnings("unchecked")
                Future<PeriodCvrpSubproblemSolver.Result>[] subFutures = (subExecutor != null)
                        ? new Future[ins.l + 1] : null;
                int cacheMissCount = 0;

                for (int t = 1; t <= ins.l; t++) {
                    iterSubCalls++;
                    totalSubCalls++;
                    String key = buildSubproblemCacheKey(point, t, ins.n);
                    subKeysByPeriod[t] = key;
                    PeriodCvrpSubproblemSolver.Result cached = subproblemCache.get(key);
                    if (cached != null) {
                        subResults[t] = cached;
                    } else {
                        cacheMissCount++;
                        iterSubCacheMisses++;
                        totalSubCacheMisses++;
                        if (subExecutor != null) {
                            final int period = t;
                            final double[] qBarT = point.qBar[t];
                            final int[] zBarT = point.zBar[t];
                            subFutures[t] = subExecutor.submit(new Callable<PeriodCvrpSubproblemSolver.Result>() {
                                @Override
                                public PeriodCvrpSubproblemSolver.Result call() {
                                    return subSolvers[period].solve(ins, period, qBarT, zBarT);
                                }
                            });
                        }
                    }
                }

                // Phase 2: collect results (parallel or sequential fallback)
                if (cacheMissCount > 0) {
                    long subStartNs = System.nanoTime();
                    if (subExecutor != null) {
                        for (int t = 1; t <= ins.l; t++) {
                            if (subFutures[t] != null) {
                                try {
                                    subResults[t] = subFutures[t].get();
                                } catch (Exception e) {
                                    throw new RuntimeException("Subproblem solve failed for t=" + t, e);
                                }
                            }
                        }
                    } else {
                        for (int t = 1; t <= ins.l; t++) {
                            if (subResults[t] == null) {
                                subResults[t] = subSolvers[t].solve(ins, t, point.qBar[t], point.zBar[t]);
                            }
                        }
                    }
                    long subElapsed = System.nanoTime() - subStartNs;
                    iterSubSolveNs += subElapsed;
                    totalSubSolveNs += subElapsed;

                    for (int t = 1; t <= ins.l; t++) {
                        if (subKeysByPeriod[t] != null && subResults[t] != null
                                && !subproblemCache.containsKey(subKeysByPeriod[t])) {
                            PeriodCvrpSubproblemSolver.Result sub = subResults[t];
                            if (sub.optimal || sub.provenInfeasible) {
                                subproblemCache.put(subKeysByPeriod[t], sub);
                            }
                        }
                    }
                }

                // Phase 3: process results
                int unresolvedPeriod = -1;
                String unresolvedSubproblemStatus = null;
                for (int t = 1; t <= ins.l; t++) {
                    PeriodCvrpSubproblemSolver.Result sub = subResults[t];
                    if (sub.provenInfeasible) {
                        infeasiblePeriod = t;
                        break;
                    }
                    if (!sub.feasible || !sub.optimal || Double.isNaN(sub.objective) || Double.isInfinite(sub.objective)) {
                        unresolvedPeriod = t;
                        unresolvedSubproblemStatus = sub.status;
                        allSubproblemsOptimal = false;
                        break;
                    }
                    phiByPeriod[t] = sub.objective;
                    phi += sub.objective;
                }

                if (infeasiblePeriod >= 0) {
                    System.out.println("[LBBD] iter=" + iteration
                            + " infeasiblePeriod=" + infeasiblePeriod
                            + " masterObj=" + fmt(point.bmpObjective)
                            + " cacheSize=" + subproblemCache.size()
                            + " lpCacheSize=" + rmpCutCache.size()
                            + " lpDualCuts=" + lpDualCuts
                            + " lpDualSkips=" + lpDualCutSkips
                            + " activeRmpCuts=" + activeRmpDualCuts.size()
                            + " activePeriodOptCuts=" + activePeriodOptCuts.size()
                            + " activeGlobalOptCuts=" + activeGlobalOptCuts.size()
                            + " activeIncPruneCuts=" + activeIncumbentPruneCuts.size()
                            + " incPruneCuts=" + incumbentPruneCuts
                            + " globalOptCuts=" + globalOptCuts
                            + " targetPruneCuts=" + targetPruneCuts);
                    boolean addedStrong = master.addStrongFeasibilityCut(point, infeasiblePeriod);
                    if (!addedStrong) {
                        master.addNoGoodCut(point);
                    }
                    feasibilityCuts++;
                    forceStrictMasterGap = false;
                    continue;
                }
                if (unresolvedPeriod >= 0) {
                    return buildSubproblemUnresolvedResult(
                            startNs,
                            iteration,
                            unresolvedPeriod,
                            unresolvedSubproblemStatus,
                            feasibilityCuts,
                            optimalityCuts,
                            lpDualCuts
                    );
                }

                double exactObjective = point.nonRouteCost + phi;
                if (exactObjective < bestUpperBound) {
                    bestUpperBound = exactObjective;
                    bestUpperBoundGapRef = point.bmpObjective;
                    setMasterUpperCutoff(cplex, bestUpperBound - INCUMBENT_PRUNE_TOL);
                }

                double delta = phi - point.omegaSum;
                boolean violatedOptimality = phi > point.omegaSum + CUT_EPS;
                boolean[] periodCoveredByFreshRmpCut = null;
                int freshRmpCutsAddedThisIter = 0;

                if (violatedOptimality) {
                    periodCoveredByFreshRmpCut = new boolean[ins.l + 1];
                    int periodsToProcess = Math.min(MAX_RMP_PERIODS_PER_ITER, ins.l);
                    boolean[] picked = new boolean[ins.l + 1];
                    for (int pick = 0; pick < periodsToProcess; pick++) {
                        int bestT = -1;
                        double bestGap = CUT_EPS;
                        for (int t = 1; t <= ins.l; t++) {
                            if (picked[t]) {
                                continue;
                            }
                            String key = subKeysByPeriod[t];
                            if (key == null || addedRmpCutKeys.contains(key)) {
                                continue;
                            }
                            double periodGap = phiByPeriod[t] - point.omega[t];
                            if (periodGap > bestGap) {
                                bestGap = periodGap;
                                bestT = t;
                            }
                        }
                        if (bestT < 0) {
                            break;
                        }
                        picked[bestT] = true;
                        String key = subKeysByPeriod[bestT];
                        PeriodRouteMasterLpSolver.Result rmpCut = rmpCutCache.get(key);
                        if (rmpCut == null) {
                            iterRmpCalls++;
                            totalRmpCalls++;
                            iterRmpCacheMisses++;
                            totalRmpCacheMisses++;
                            long rmpStartNs = System.nanoTime();
                            rmpCut = rmpLpSolver.solve(
                                    ins,
                                    bestT,
                                    point.qBar[bestT],
                                    point.zBar[bestT],
                                    buildPrevVisitForPeriod(point, bestT, ins.n)
                            );
                            long rmpElapsed = System.nanoTime() - rmpStartNs;
                            iterRmpSolveNs += rmpElapsed;
                            totalRmpSolveNs += rmpElapsed;
                            if (rmpCut.feasible && rmpCut.optimal && rmpCut.pricingProvedOptimal) {
                                rmpCutCache.put(key, rmpCut);
                            }
                        }
                        if (rmpCut == null
                                || !rmpCut.feasible
                                || !rmpCut.optimal
                                || !rmpCut.pricingProvedOptimal
                                || !rmpCut.artificialClean
                                || rmpCut.lpObjective <= point.omega[bestT] + CUT_EPS) {
                            continue;
                        }

                        double cutAtPoint = ins.K * rmpCut.dualU0;
                        for (int i = 1; i <= ins.n; i++) {
                            for (int v = ins.pi[i][bestT]; v <= bestT - 1; v++) {
                                if (!point.lambdaOne[i][v][bestT]) {
                                    continue;
                                }
                                if (rmpCut.dualUiByPrev != null
                                        && i < rmpCut.dualUiByPrev.length
                                        && v < rmpCut.dualUiByPrev[i].length) {
                                    cutAtPoint += rmpCut.dualUiByPrev[i][v];
                                }
                            }
                        }
                        double cutPointTol = 1e-5 * Math.max(1.0, Math.abs(rmpCut.lpObjective));
                        boolean validAtPoint = rmpCut.dualU0 <= 1e-7
                                && Math.abs(cutAtPoint - rmpCut.lpObjective) <= cutPointTol
                                && cutAtPoint <= rmpCut.lpObjective + cutPointTol;
                        if (!validAtPoint) {
                            lpDualCutSkips++;
                            continue;
                        }
                        IloRange addedRange = master.addPeriodRmpDualCut(bestT, point, rmpCut);
                        addedRmpCutKeys.add(key);
                        activeRmpDualCuts.addLast(new ActiveRmpDualCut(key, addedRange));
                        trimActiveRmpDualCuts(master, activeRmpDualCuts, addedRmpCutKeys, activeRmpDualCutCap(iteration));
                        periodCoveredByFreshRmpCut[bestT] = true;
                        freshRmpCutsAddedThisIter++;
                        lpDualCuts++;
                    }
                }

                if (shouldLogIteration(iteration)) {
                    System.out.println("[LBBD] iter=" + iteration
                            + " masterObj=" + fmt(point.bmpObjective)
                            + " exactObj=" + fmt(exactObjective)
                            + " omegaSum=" + fmt(point.omegaSum)
                            + " phi=" + fmt(phi)
                            + " delta(phi-omega)=" + fmt(delta)
                            + " cacheSize=" + subproblemCache.size()
                            + " lpCacheSize=" + rmpCutCache.size()
                            + " lpDualCuts=" + lpDualCuts
                            + " lpDualSkips=" + lpDualCutSkips
                            + " activeRmpCuts=" + activeRmpDualCuts.size()
                            + " activePeriodOptCuts=" + activePeriodOptCuts.size()
                            + " activeGlobalOptCuts=" + activeGlobalOptCuts.size()
                            + " activeIncPruneCuts=" + activeIncumbentPruneCuts.size()
                            + " incPruneCuts=" + incumbentPruneCuts
                            + " globalOptCuts=" + globalOptCuts
                            + " targetPruneCuts=" + targetPruneCuts
                            + formatTimingSuffix(
                            iteration,
                            iterMasterSolveNs,
                            iterSubSolveNs,
                            iterRmpSolveNs,
                            iterSubCalls,
                            iterSubCacheMisses,
                            iterRmpCalls,
                            iterRmpCacheMisses,
                            totalMasterSolveNs,
                            totalSubSolveNs,
                            totalRmpSolveNs,
                            totalSubCalls,
                            totalSubCacheMisses,
                            totalRmpCalls,
                            totalRmpCacheMisses
                    ));
                }

                if (!Double.isNaN(targetObjective) && Math.abs(exactObjective - targetObjective) <= targetTol) {
                    String status = "LBBD_TargetMatched(iter=" + iteration
                            + ",target=" + fmt(targetObjective)
                            + ",feasCuts=" + feasibilityCuts
                            + ",optCuts=" + optimalityCuts
                            + ",lpDualCuts=" + lpDualCuts
                            + ",lpDualSkips=" + lpDualCutSkips
                            + ",incPruneCuts=" + incumbentPruneCuts
                            + ",targetPruneCuts=" + targetPruneCuts + ")";
                    double sec = elapsedSec(startNs);
                    return new SolveResult(
                            "LBBDReformulation",
                            true,
                            false,
                            status,
                            exactObjective,
                            point.bmpObjective,
                            relativeGap(exactObjective, point.bmpObjective),
                            sec
                    );
                }
                if (!Double.isNaN(targetObjective) && exactObjective < targetObjective - targetTol) {
                    throw new IllegalStateException(
                            "LBBD found objective better than reform target: exact="
                                    + exactObjective + ", target=" + targetObjective
                    );
                }

                if (violatedOptimality) {
                    if (maxPeriodOptCutsPerIter > 0) {
                        int periodOptBudget = (freshRmpCutsAddedThisIter == 0)
                                ? Math.max(maxPeriodOptCutsPerIter, ins.l)
                                : maxPeriodOptCutsPerIter;
                        if (iteration >= latePeriodOptStartIter && delta >= latePeriodOptMinDelta) {
                            periodOptBudget = Math.min(ins.l, Math.max(periodOptBudget, latePeriodOptCutsPerIter));
                        }
                        int addedPeriodOptThisIter = 0;
                        while (addedPeriodOptThisIter < periodOptBudget) {
                            int bestT = -1;
                            double bestGap = CUT_EPS;
                            for (int t = 1; t <= ins.l; t++) {
                                if (phiByPeriod[t] <= point.omega[t] + CUT_EPS) {
                                    continue;
                                }
                                String periodKey = buildPeriodOptCutKey(subKeysByPeriod[t], t);
                                if (periodKey == null || addedPeriodOptCutKeys.contains(periodKey)) {
                                    continue;
                                }
                                double gapT = phiByPeriod[t] - point.omega[t];
                                if (gapT > bestGap) {
                                    bestGap = gapT;
                                    bestT = t;
                                }
                            }
                            if (bestT < 0) {
                                break;
                            }
                            String periodKey = buildPeriodOptCutKey(subKeysByPeriod[bestT], bestT);
                            IloRange addedRange = master.addSinglePeriodOptimalityCut(point, phiByPeriod[bestT], bestT);
                            if (addedRange == null) {
                                break;
                            }
                            addedPeriodOptCutKeys.add(periodKey);
                            activePeriodOptCuts.addLast(new ActivePeriodOptCut(periodKey, addedRange));
                            trimActivePeriodOptCuts(master, activePeriodOptCuts, addedPeriodOptCutKeys, activePeriodOptCutCap(iteration));
                            addedPeriodOptThisIter++;
                        }
                    }
                    if (ENABLE_GLOBAL_COMBINATION_OPT_CUT) {
                        if (iteration >= globalOptCutStartIter && delta >= globalOptCutMinDelta) {
                            String globalOptKey = buildGlobalOptCutKey(point, ins.n, ins.l);
                            if (!addedGlobalOptCutKeys.contains(globalOptKey)) {
                                IloRange range = master.addOptimalityCutRange(point, phi, lbR);
                                if (range != null) {
                                    addedGlobalOptCutKeys.add(globalOptKey);
                                    activeGlobalOptCuts.addLast(new ActiveGlobalOptCut(globalOptKey, range));
                                    trimActiveGlobalOptCuts(master, activeGlobalOptCuts, addedGlobalOptCutKeys, activeGlobalOptCutCap(iteration));
                                    globalOptCuts++;
                                }
                            }
                        }
                    }
                    optimalityCuts++;

                    if (isFinite(bestUpperBound)
                            && exactObjective >= bestUpperBound - INCUMBENT_PRUNE_TOL) {
                    if (ENABLE_INCUMBENT_PRUNE_NOGOOD
                            && iteration >= incumbentPruneStartIter
                            && exactObjective >= bestUpperBound + incumbentPruneMinExcess) {
                            String incPruneKey = buildGlobalRoutingPatternKey(point, ins.n, ins.l);
                            if (!activeIncumbentPruneCutKeys.contains(incPruneKey)) {
                                IloRange ng = master.addNoGoodCutRange(point);
                                activeIncumbentPruneCutKeys.add(incPruneKey);
                                activeIncumbentPruneCuts.addLast(new ActiveIncumbentPruneCut(incPruneKey, ng));
                                trimActiveIncPruneCuts(master, activeIncumbentPruneCuts, activeIncumbentPruneCutKeys, activeIncPruneCutCap(iteration));
                                incumbentPruneCuts++;
                            }
                        }
                        continue;
                    }
                }

                if (!Double.isNaN(targetObjective) && exactObjective > targetObjective + targetTol) {
                    master.addNoGoodCut(point);
                    targetPruneCuts++;
                    forceStrictMasterGap = false;
                    continue;
                }

                if (violatedOptimality) {
                    forceStrictMasterGap = false;
                    continue;
                }

                if (masterMipGapThisIter > CplexConfig.MIP_GAP + 1e-12) {
                    forceStrictMasterGap = true;
                    continue;
                }

                double finalGap = relativeGap(exactObjective, point.bmpObjective);
                boolean optimal = allSubproblemsOptimal && finalGap <= RESULT_TOL;
                String status = "LBBD_Converged(iter=" + iteration
                        + ",feasCuts=" + feasibilityCuts
                        + ",optCuts=" + optimalityCuts
                        + ",lpDualCuts=" + lpDualCuts
                        + ",lpDualSkips=" + lpDualCutSkips
                        + ",incPruneCuts=" + incumbentPruneCuts
                        + ",targetPruneCuts=" + targetPruneCuts
                        + ",masterStatus=" + masterStatus + ")";

                double sec = elapsedSec(startNs);
                return new SolveResult(
                        "LBBDReformulation",
                        true,
                        optimal,
                        status,
                        exactObjective,
                        point.bmpObjective,
                        finalGap,
                        sec
                );
            }

            double sec = elapsedSec(startNs);
            boolean feasible = !Double.isInfinite(bestUpperBound) && !Double.isNaN(bestUpperBound);
            String status = "LBBD_MaxIterations(iter=" + iteration
                    + ",feasCuts=" + feasibilityCuts
                    + ",optCuts=" + optimalityCuts
                    + ",lpDualCuts=" + lpDualCuts
                    + ",lpDualSkips=" + lpDualCutSkips
                    + ",incPruneCuts=" + incumbentPruneCuts
                    + ",targetPruneCuts=" + targetPruneCuts + ")";
            double objective = feasible ? bestUpperBound : Double.NaN;
            double bestBound = Double.isInfinite(bestLowerBound) ? Double.NaN : bestLowerBound;
            double gap = (feasible && !Double.isNaN(bestBound))
                    ? relativeGap(objective, bestBound)
                    : Double.NaN;
            if (feasible && !Double.isNaN(bestUpperBoundGapRef)) {
                gap = relativeGap(bestUpperBound, bestUpperBoundGapRef);
            }

            return new SolveResult("LBBDReformulation", feasible, false, status, objective, bestBound, gap, sec);
        } catch (IloException e) {
            throw new RuntimeException("Failed to solve LBBDReformulation", e);
        } finally {
            if (subExecutor != null) {
                subExecutor.shutdownNow();
            }
            for (int tt = 1; tt < subSolvers.length; tt++) {
                if (subSolvers[tt] != null) {
                    subSolvers[tt].close();
                }
            }
        }
    }

    /**
     * Skip degenerate initial cuts: all-zero coefficients imply omega[t] >= 0.
     */
    private boolean hasMeaningfulDualCoeffs(int t, double[][] dualW, double dualU0) {
        if (Math.abs(dualU0) > 1e-9) {
            return true;
        }
        if (dualW == null) {
            return false;
        }
        for (int i = 1; i <= ins.n; i++) {
            for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                if (i < dualW.length && dualW[i] != null && v < dualW[i].length && Math.abs(dualW[i][v]) > 1e-9) {
                    return true;
                }
            }
        }
        return false;
    }

    private SolveResult tryValidateByHint(
            Instance ins,
            double targetObjective,
            double targetTol,
            int[][] prevVisitHint,
            long startNs
    ) {
        if (prevVisitHint.length < ins.n + 1) {
            return null;
        }
        try (PeriodCvrpSubproblemSolver subproblemSolver = new PeriodCvrpSubproblemSolver()) {
            double phi = 0.0;
            for (int t = 1; t <= ins.l; t++) {
                double[] qBar = new double[ins.n + 1];
                int[] zBar = new int[ins.n + 1];
                for (int i = 1; i <= ins.n; i++) {
                    if (prevVisitHint[i] == null || prevVisitHint[i].length <= t) {
                        return null;
                    }
                    int v = prevVisitHint[i][t];
                    if (v < 0) {
                        continue;
                    }
                    if (v < ins.pi[i][t] || v > t - 1) {
                        return null;
                    }
                    zBar[i] = 1;
                    qBar[i] = ins.g(i, v, t);
                }
                PeriodCvrpSubproblemSolver.Result sub = subproblemSolver.solve(ins, t, qBar, zBar);
                if (!sub.feasible || !sub.optimal) {
                    return null;
                }
                phi += sub.objective;
            }

            double sec = elapsedSec(startNs);
            String status = "LBBD_TargetMatchedByHint(target=" + fmt(targetObjective)
                    + ",phi=" + fmt(phi) + ")";
            return new SolveResult(
                    "LBBDReformulation",
                    true,
                    false,
                    status,
                    targetObjective,
                    Double.NaN,
                    Double.NaN,
                    sec
            );
        }
    }

    private static boolean shouldLogIteration(int iteration) {
        return iteration <= 10 || iteration % 10 == 0;
    }

    private static String buildSubproblemCacheKey(MasterPoint point, int t, int n) {
        StringBuilder sb = new StringBuilder(n * 3 + 4);
        sb.append(t).append(':');
        for (int i = 1; i <= n; i++) {
            if (i > 1) {
                sb.append(',');
            }
            sb.append(point.prevVisit[i][t]);
        }
        return sb.toString();
    }

    private static int[] buildPrevVisitForPeriod(MasterPoint point, int t, int n) {
        int[] prev = new int[n + 1];
        for (int i = 1; i <= n; i++) {
            prev[i] = point.prevVisit[i][t];
        }
        return prev;
    }

    private static String buildPeriodOptCutKey(String subKey, int t) {
        if (subKey == null) {
            return null;
        }
        return "P|" + t + "|" + subKey;
    }

    private static String buildGlobalRoutingPatternKey(MasterPoint point, int n, int l) {
        StringBuilder sb = new StringBuilder(l * (n * 3 + 4));
        for (int t = 1; t <= l; t++) {
            if (t > 1) {
                sb.append('|');
            }
            sb.append(t).append(':');
            for (int i = 1; i <= n; i++) {
                if (i > 1) {
                    sb.append(',');
                }
                sb.append(point.prevVisit[i][t]);
            }
        }
        return sb.toString();
    }

    private static String buildGlobalOptCutKey(MasterPoint point, int n, int l) {
        return "GOPT|" + buildGlobalRoutingPatternKey(point, n, l);
    }

    private static SolveResult buildFailedResult(
            long startNs,
            int iteration,
            int feasibilityCuts,
            int optimalityCuts,
            int lpDualCuts,
            String masterStatus
    ) {
        double sec = elapsedSec(startNs);
        String status = "LBBD_MasterSolveFailed(iter=" + iteration
                + ",feasCuts=" + feasibilityCuts
                + ",optCuts=" + optimalityCuts
                + ",lpDualCuts=" + lpDualCuts
                + ",masterStatus=" + masterStatus + ")";
        return new SolveResult("LBBDReformulation", false, false, status, Double.NaN, Double.NaN, Double.NaN, sec);
    }

    private static SolveResult buildSubproblemUnresolvedResult(
            long startNs,
            int iteration,
            int period,
            String subStatus,
            int feasibilityCuts,
            int optimalityCuts,
            int lpDualCuts
    ) {
        double sec = elapsedSec(startNs);
        String status = "LBBD_SubproblemUnresolved(iter=" + iteration
                + ",period=" + period
                + ",feasCuts=" + feasibilityCuts
                + ",optCuts=" + optimalityCuts
                + ",lpDualCuts=" + lpDualCuts
                + ",subStatus=" + subStatus + ")";
        return new SolveResult("LBBDReformulation", false, false, status, Double.NaN, Double.NaN, Double.NaN, sec);
    }

    private static SolveResult buildCutoffExhaustedResult(
            long startNs,
            int iteration,
            int feasibilityCuts,
            int optimalityCuts,
            int lpDualCuts,
            boolean allSubproblemsOptimal,
            double bestUpperBound,
            double bestLowerBound
    ) {
        double sec = elapsedSec(startNs);
        double impliedLowerBound = bestUpperBound - INCUMBENT_PRUNE_TOL;
        double finalBestBound = isFinite(bestLowerBound) ? Math.max(bestLowerBound, impliedLowerBound) : impliedLowerBound;
        double finalGap = relativeGap(bestUpperBound, finalBestBound);
        boolean optimal = allSubproblemsOptimal && finalGap <= RESULT_TOL;
        String status = "LBBD_CutoffExhausted(iter=" + iteration
                + ",feasCuts=" + feasibilityCuts
                + ",optCuts=" + optimalityCuts
                + ",lpDualCuts=" + lpDualCuts
                + ",cutoff=" + fmt(impliedLowerBound) + ")";
        return new SolveResult(
                "LBBDReformulation",
                true,
                optimal,
                status,
                bestUpperBound,
                finalBestBound,
                finalGap,
                sec
        );
    }

    private static double elapsedSec(long startNs) {
        return (System.nanoTime() - startNs) / 1_000_000_000.0;
    }

    private static String fmt(double v) {
        return String.format(Locale.US, "%.6f", v);
    }

    private static String formatTimingSuffix(
            int iteration,
            long iterMasterSolveNs,
            long iterSubSolveNs,
            long iterRmpSolveNs,
            int iterSubCalls,
            int iterSubCacheMisses,
            int iterRmpCalls,
            int iterRmpCacheMisses,
            long totalMasterSolveNs,
            long totalSubSolveNs,
            long totalRmpSolveNs,
            long totalSubCalls,
            long totalSubCacheMisses,
            long totalRmpCalls,
            long totalRmpCacheMisses
    ) {
        if (!ENABLE_ITER_TIMING_LOG) {
            return "";
        }
        double iterMasterMs = iterMasterSolveNs / 1_000_000.0;
        double iterSubMs = iterSubSolveNs / 1_000_000.0;
        double iterRmpMs = iterRmpSolveNs / 1_000_000.0;
        double avgMasterMs = iteration > 0 ? (totalMasterSolveNs / 1_000_000.0) / iteration : 0.0;
        double avgSubMs = iteration > 0 ? (totalSubSolveNs / 1_000_000.0) / iteration : 0.0;
        double avgRmpMs = iteration > 0 ? (totalRmpSolveNs / 1_000_000.0) / iteration : 0.0;
        return String.format(
                Locale.US,
                " | timeMs(iter: master=%.1f,sub=%.1f,rmp=%.1f; avg: master=%.1f,sub=%.1f,rmp=%.1f) | cacheMiss(iter: sub=%d/%d,rmp=%d) | cacheMissTotal(sub=%d/%d,rmp=%d)",
                iterMasterMs,
                iterSubMs,
                iterRmpMs,
                avgMasterMs,
                avgSubMs,
                avgRmpMs,
                iterSubCacheMisses,
                iterSubCalls,
                iterRmpCacheMisses,
                totalSubCacheMisses,
                totalSubCalls,
                totalRmpCacheMisses
        );
    }

    private static double relativeGap(double upper, double lower) {
        if (Double.isNaN(upper) || Double.isNaN(lower)) {
            return Double.NaN;
        }
        double denom = Math.max(1.0, Math.abs(upper));
        return Math.abs(upper - lower) / denom;
    }

    private static boolean isInfeasibleLikeStatus(String status) {
        return status != null && status.startsWith("Infeasible");
    }

    private static boolean isFinite(double value) {
        return !Double.isNaN(value) && !Double.isInfinite(value);
    }

    private static void setMasterUpperCutoff(IloCplex cplex, double cutoff) throws IloException {
        if (!isFinite(cutoff)) {
            return;
        }
        cplex.setParam(IloCplex.DoubleParam.CutUp, cutoff);
    }

    private static String safeStatus(IloCplex cplex) {
        try {
            return cplex.getStatus().toString();
        } catch (IloException e) {
            return "Unknown";
        }
    }

    private double recommendedMasterMipGap(int iteration, boolean forceStrict) {
        if (forceStrict) {
            return CplexConfig.MIP_GAP;
        }
        if (iteration <= masterMipGapEarlyIterCutoff) {
            return masterMipGapEarly;
        }
        if (iteration <= masterMipGapMidIterCutoff) {
            return masterMipGapMid;
        }
        if (iteration <= masterMipGapVeryLateIterCutoff) {
            return masterMipGapLate;
        }
        return masterMipGapVeryLate;
    }

    private int activeRmpDualCutCap(int iteration) {
        return iteration >= activeCutCapLateStartIter
                ? activeRmpDualCutsCapLate
                : activeRmpDualCutsCapBase;
    }

    private int activePeriodOptCutCap(int iteration) {
        return iteration >= activeCutCapLateStartIter
                ? activePeriodOptCutsCapLate
                : activePeriodOptCutsCapBase;
    }

    private int activeGlobalOptCutCap(int iteration) {
        return iteration >= activeCutCapLateStartIter
                ? activeGlobalOptCutsCapLate
                : activeGlobalOptCutsCapBase;
    }

    private int activeIncPruneCutCap(int iteration) {
        return iteration >= activeCutCapLateStartIter
                ? activeIncPruneCutsCapLate
                : activeIncPruneCutsCapBase;
    }

    private static void trimActiveRmpDualCuts(
            MasterModel master,
            ArrayDeque<ActiveRmpDualCut> activeCuts,
            HashSet<String> keys,
            int cap
    ) throws IloException {
        while (activeCuts.size() > cap) {
            ActiveRmpDualCut evicted = activeCuts.removeFirst();
            master.deleteCut(evicted.range);
            keys.remove(evicted.key);
        }
    }

    private static void trimActivePeriodOptCuts(
            MasterModel master,
            ArrayDeque<ActivePeriodOptCut> activeCuts,
            HashSet<String> keys,
            int cap
    ) throws IloException {
        while (activeCuts.size() > cap) {
            ActivePeriodOptCut evicted = activeCuts.removeFirst();
            master.deleteCut(evicted.range);
            keys.remove(evicted.key);
        }
    }

    private static void trimActiveGlobalOptCuts(
            MasterModel master,
            ArrayDeque<ActiveGlobalOptCut> activeCuts,
            HashSet<String> keys,
            int cap
    ) throws IloException {
        while (activeCuts.size() > cap) {
            ActiveGlobalOptCut evicted = activeCuts.removeFirst();
            master.deleteCut(evicted.range);
            keys.remove(evicted.key);
        }
    }

    private static void trimActiveIncPruneCuts(
            MasterModel master,
            ArrayDeque<ActiveIncumbentPruneCut> activeCuts,
            HashSet<String> keys,
            int cap
    ) throws IloException {
        while (activeCuts.size() > cap) {
            ActiveIncumbentPruneCut evicted = activeCuts.removeFirst();
            master.deleteCut(evicted.range);
            keys.remove(evicted.key);
        }
    }

    private static boolean isPickupNode(int node, int n) {
        return node >= 1 && node <= n;
    }

    private static double holdingCostOnArc(Instance ins, int i, int v, int t) {
        return ins.e(i, v, t);
    }

    private static void configure(IloCplex cplex) throws IloException {
        cplex.setParam(IloCplex.Param.TimeLimit, CplexConfig.TIME_LIMIT_SEC);
        cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, CplexConfig.MIP_GAP);
        if (MASTER_THREADS > 0) {
            cplex.setParam(IloCplex.Param.Threads, MASTER_THREADS);
        }
        if (!LOG_TO_CONSOLE) {
            cplex.setOut(null);
            cplex.setWarning(null);
        }
    }

    private static final class LambdaRef {
        final int i;
        final int v;
        final int t;
        final IloNumVar var;

        LambdaRef(int i, int v, int t, IloNumVar var) {
            this.i = i;
            this.v = v;
            this.t = t;
            this.var = var;
        }
    }

    private static final class MasterPoint {
        final double bmpObjective;
        final double bmpBestBound;
        final double nonRouteCost;
        final double omegaSum;
        final double[] omega;
        final double[][] qBar; // [t][i]
        final int[][] zBar;    // [t][i]
        final int[][] prevVisit; // [i][t]
        final boolean[][][] lambdaOne; // [i][v][t], t in 1..l only
        final int routeLambdaOneCount;

        MasterPoint(
                double bmpObjective,
                double bmpBestBound,
                double nonRouteCost,
                double omegaSum,
                double[] omega,
                double[][] qBar,
                int[][] zBar,
                int[][] prevVisit,
                boolean[][][] lambdaOne,
                int routeLambdaOneCount
        ) {
            this.bmpObjective = bmpObjective;
            this.bmpBestBound = bmpBestBound;
            this.nonRouteCost = nonRouteCost;
            this.omegaSum = omegaSum;
            this.omega = omega;
            this.qBar = qBar;
            this.zBar = zBar;
            this.prevVisit = prevVisit;
            this.lambdaOne = lambdaOne;
            this.routeLambdaOneCount = routeLambdaOneCount;
        }
    }

    private static final class MasterModel {
        private final IloCplex cplex;
        private final Instance ins;

        private final IloNumVar[] omega;
        private final IloNumVar[] m;
        private final IloNumVar[][][] lambda;
        private final ArrayList<LambdaRef> routingLambdaRefs;

        private int feasCutCount = 0;
        private int noGoodCutCount = 0;
        private int optCutCount = 0;

        MasterModel(IloCplex cplex, Instance ins) throws IloException {
            this.cplex = cplex;
            this.ins = ins;
            this.omega = new IloNumVar[ins.l + 1];
            this.m = new IloNumVar[ins.l + 1];
            this.lambda = new IloNumVar[ins.n + 1][ins.l + 2][ins.l + 2];
            this.routingLambdaRefs = new ArrayList<LambdaRef>();
            build();
        }

        private void build() throws IloException {
            int n = ins.n;
            int l = ins.l;
            double[] minIncomingToSupplier = new double[n + 1];
            double[] minOutgoingFromSupplier = new double[n + 1];
            for (int j = 1; j <= n; j++) {
                double minIn = Double.POSITIVE_INFINITY;
                double minOut = Double.POSITIVE_INFINITY;
                for (int i = 0; i < ins.nodeCount; i++) {
                    if (i == j) {
                        continue;
                    }
                    if (ins.c[i][j] < minIn) {
                        minIn = ins.c[i][j];
                    }
                    if (ins.c[j][i] < minOut) {
                        minOut = ins.c[j][i];
                    }
                }
                minIncomingToSupplier[j] = minIn;
                minOutgoingFromSupplier[j] = minOut;
            }
            double minRoundTrip = Double.POSITIVE_INFINITY;
            for (int j = 1; j <= n; j++) {
                double rt = ins.c[0][j] + ins.c[j][n + 1];
                if (rt < minRoundTrip) {
                    minRoundTrip = rt;
                }
            }
            double loadBasedOmegaCoeff = (minRoundTrip > 0.0 && ins.Q > 0.0) ? (minRoundTrip / ins.Q) : 0.0;

            IloNumVar[] y = new IloNumVar[l + 1];
            IloNumVar[] p = new IloNumVar[l + 1];
            IloNumVar[] p0 = new IloNumVar[l + 1];
            IloNumVar[] i0 = new IloNumVar[l + 1];

            for (int t = 1; t <= l; t++) {
                y[t] = cplex.boolVar("y_" + t);
                p[t] = cplex.numVar(0.0, Double.MAX_VALUE, "p_" + t);
                p0[t] = cplex.numVar(0.0, Double.MAX_VALUE, "P0_" + t);
                i0[t] = cplex.numVar(0.0, Double.MAX_VALUE, "I0_" + t);
                m[t] = cplex.intVar(0, ins.K, "m_" + t);
                omega[t] = cplex.numVar(0.0, Double.MAX_VALUE, "omega_" + t);
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l + 1; t++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (t <= l && ins.g(i, v, t) > ins.Q + 1e-9) {
                            continue;
                        }
                        IloNumVar var = cplex.boolVar("lambda_" + i + "_" + v + "_" + t);
                        lambda[i][v][t] = var;
                        if (t <= l) {
                            routingLambdaRefs.add(new LambdaRef(i, v, t, var));
                        }
                    }
                }
            }

            IloLinearNumExpr obj = cplex.linearNumExpr();
            for (int t = 1; t <= l; t++) {
                obj.addTerm(1.0, omega[t]);
                obj.addTerm(ins.u, p[t]);
                obj.addTerm(ins.f, y[t]);
                obj.addTerm(ins.h0, i0[t]);
                obj.addTerm(ins.hp, p0[t]);
            }
            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l + 1; t++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (lambda[i][v][t] != null) {
                            obj.addTerm(holdingCostOnArc(ins, i, v, t), lambda[i][v][t]);
                        }
                    }
                }
            }
            cplex.addMinimize(obj);

            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr finishedBalance = cplex.linearNumExpr();
                double rhsFinished = ins.dt[t];
                if (t == 1) {
                    rhsFinished -= ins.P00;
                } else {
                    finishedBalance.addTerm(1.0, p0[t - 1]);
                }
                finishedBalance.addTerm(1.0, p[t]);
                finishedBalance.addTerm(-1.0, p0[t]);
                cplex.addEq(finishedBalance, rhsFinished, "FinishedBalance_" + t);

                IloLinearNumExpr prodCap = cplex.linearNumExpr();
                prodCap.addTerm(1.0, p[t]);
                prodCap.addTerm(-ins.C, y[t]);
                cplex.addLe(prodCap, 0.0, "ProdCap_" + t);

                cplex.addLe(p0[t], ins.Lp, "FinishedCap_" + t);
            }

            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr factoryBalance = cplex.linearNumExpr();
                double rhsFactory = 0.0;
                if (t == 1) {
                    rhsFactory = -ins.I00;
                } else {
                    factoryBalance.addTerm(1.0, i0[t - 1]);
                }
                for (int i = 1; i <= n; i++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (lambda[i][v][t] != null) {
                            factoryBalance.addTerm(ins.g(i, v, t), lambda[i][v][t]);
                        }
                    }
                }
                factoryBalance.addTerm(-ins.k, p[t]);
                factoryBalance.addTerm(-1.0, i0[t]);
                cplex.addEq(factoryBalance, rhsFactory, "FactoryBalance_" + t);

                cplex.addLe(i0[t], ins.L0, "FactoryCap_" + t);
            }

            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr vehCap = cplex.linearNumExpr();
                for (int i = 1; i <= n; i++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (lambda[i][v][t] != null) {
                            vehCap.addTerm(ins.g(i, v, t), lambda[i][v][t]);
                        }
                    }
                }
                vehCap.addTerm(-ins.Q, m[t]);
                cplex.addLe(vehCap, 0.0, "VehicleCap_" + t);
                cplex.addLe(m[t], ins.K, "VehicleCount_" + t);

                IloLinearNumExpr lbIn = cplex.linearNumExpr();
                IloLinearNumExpr lbOut = cplex.linearNumExpr();
                IloLinearNumExpr lbLoad = cplex.linearNumExpr();
                IloLinearNumExpr lbM = cplex.linearNumExpr();
                lbIn.addTerm(1.0, omega[t]);
                lbOut.addTerm(1.0, omega[t]);
                lbLoad.addTerm(1.0, omega[t]);
                lbM.addTerm(1.0, omega[t]);
                lbM.addTerm(-minRoundTrip, m[t]);

                for (int i = 1; i <= n; i++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (lambda[i][v][t] != null) {
                            lbIn.addTerm(-minIncomingToSupplier[i], lambda[i][v][t]);
                            lbOut.addTerm(-minOutgoingFromSupplier[i], lambda[i][v][t]);
                            lbLoad.addTerm(-loadBasedOmegaCoeff * ins.g(i, v, t), lambda[i][v][t]);
                        }
                    }
                }
                cplex.addGe(lbIn, 0.0, "OmegaLB_In_" + t);
                cplex.addGe(lbOut, 0.0, "OmegaLB_Out_" + t);
                cplex.addGe(lbLoad, 0.0, "OmegaLB_Load_" + t);
                cplex.addGe(lbM, 0.0, "OmegaLB_M_" + t);
            }

            for (int i = 1; i <= n; i++) {
                IloLinearNumExpr startFlow = cplex.linearNumExpr();
                for (int t = 1; t <= ins.mu[i][0]; t++) {
                    if (lambda[i][0][t] != null) {
                        startFlow.addTerm(1.0, lambda[i][0][t]);
                    }
                }
                cplex.addEq(startFlow, 1.0, "StartFlow_" + i);
            }

            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l; t++) {
                    IloLinearNumExpr balance = cplex.linearNumExpr();
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (lambda[i][v][t] != null) {
                            balance.addTerm(1.0, lambda[i][v][t]);
                        }
                    }
                    for (int tau = t + 1; tau <= ins.mu[i][t]; tau++) {
                        if (lambda[i][t][tau] != null) {
                            balance.addTerm(-1.0, lambda[i][t][tau]);
                        }
                    }
                    cplex.addEq(balance, 0.0, "LambdaFlow_" + i + "_" + t);
                }
            }

            for (int i = 1; i <= n; i++) {
                IloLinearNumExpr terminalFlow = cplex.linearNumExpr();
                for (int t = ins.pi[i][l + 1]; t <= l; t++) {
                    if (lambda[i][t][l + 1] != null) {
                        terminalFlow.addTerm(1.0, lambda[i][t][l + 1]);
                    }
                }
                cplex.addEq(terminalFlow, 1.0, "TerminalFlow_" + i);
            }

            // === Valid Inequality A.1: Minimum collection start time (eq 4-1/4-2) ===
            {
                double initialCap = ins.P00;
                if (Math.abs(ins.k) <= 1e-12) {
                    initialCap = Double.POSITIVE_INFINITY;
                } else {
                    initialCap += ins.I00 / ins.k;
                }
                double cumDemand = 0.0;
                int tPrime = -1;
                for (int t = 1; t <= l; t++) {
                    cumDemand += ins.dt[t];
                    if (initialCap < cumDemand - 1e-9) {
                        tPrime = t;
                        break;
                    }
                }
                if (tPrime > 0) {
                    IloLinearNumExpr expr = cplex.linearNumExpr();
                    for (int tau = 1; tau <= tPrime; tau++) {
                        for (int i = 1; i <= n; i++) {
                            for (int v = ins.pi[i][tau]; v <= tau - 1; v++) {
                                if (lambda[i][v][tau] != null) {
                                    expr.addTerm(1.0, lambda[i][v][tau]);
                                }
                            }
                        }
                    }
                    cplex.addGe(expr, 1.0, "VI_MinCollectStart");
                }
            }

            // === Valid Inequality A.2: Per-period total capacity (eq 4-3) ===
            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr expr = cplex.linearNumExpr();
                boolean hasTerms = false;
                for (int i = 1; i <= n; i++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (lambda[i][v][t] != null) {
                            expr.addTerm(ins.g(i, v, t), lambda[i][v][t]);
                            hasTerms = true;
                        }
                    }
                }
                if (hasTerms) {
                    cplex.addLe(expr, (double) ins.K * ins.Q, "VI_TotalCap_" + t);
                }
            }

            // === Valid Inequality A.3: Large-item constraint (eq 4-4) ===
            for (int t = 1; t <= l; t++) {
                IloLinearNumExpr expr = cplex.linearNumExpr();
                boolean hasTerms = false;
                for (int i = 1; i <= n; i++) {
                    for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                        if (lambda[i][v][t] != null) {
                            double gVal = ins.g(i, v, t);
                            if (gVal > ins.Q / 2.0 + 1e-9) {
                                expr.addTerm(1.0, lambda[i][v][t]);
                                hasTerms = true;
                            }
                        }
                    }
                }
                if (hasTerms) {
                    cplex.addLe(expr, ins.K, "VI_BigItem_" + t);
                }
            }
        }

        MasterPoint extractPoint(Instance ins) throws IloException {
            int n = ins.n;
            int l = ins.l;

            double[] omegaVal = new double[l + 1];
            double omegaSum = 0.0;
            for (int t = 1; t <= l; t++) {
                omegaVal[t] = cplex.getValue(omega[t]);
                omegaSum += omegaVal[t];
            }

            boolean[][][] lambdaOne = new boolean[n + 1][l + 1][l + 1];
            int[][] prevVisit = new int[n + 1][l + 1];
            for (int i = 1; i <= n; i++) {
                for (int t = 1; t <= l; t++) {
                    prevVisit[i][t] = -1;
                }
            }
            double[][] qBar = new double[l + 1][n + 1];
            int[][] zBar = new int[l + 1][n + 1];

            int routeLambdaOneCount = 0;
            for (int idx = 0; idx < routingLambdaRefs.size(); idx++) {
                LambdaRef ref = routingLambdaRefs.get(idx);
                double value = cplex.getValue(ref.var);
                if (value > 0.5) {
                    if (!lambdaOne[ref.i][ref.v][ref.t]) {
                        lambdaOne[ref.i][ref.v][ref.t] = true;
                        routeLambdaOneCount++;
                    }
                    if (prevVisit[ref.i][ref.t] >= 0 && prevVisit[ref.i][ref.t] != ref.v) {
                        throw new IllegalStateException(
                                "Multiple active lambda predecessors for i=" + ref.i + ", t=" + ref.t
                        );
                    }
                    prevVisit[ref.i][ref.t] = ref.v;
                    zBar[ref.t][ref.i] = 1;
                    qBar[ref.t][ref.i] += ins.g(ref.i, ref.v, ref.t);
                }
            }

            double bmpObjective = cplex.getObjValue();
            double bmpBestBound;
            try {
                bmpBestBound = cplex.getBestObjValue();
            } catch (IloException ignored) {
                bmpBestBound = bmpObjective;
            }
            double nonRouteCost = bmpObjective - omegaSum;

            return new MasterPoint(
                    bmpObjective,
                    bmpBestBound,
                    nonRouteCost,
                    omegaSum,
                    omegaVal,
                    qBar,
                    zBar,
                    prevVisit,
                    lambdaOne,
                    routeLambdaOneCount
            );
        }

        boolean addStrongFeasibilityCut(MasterPoint point, int t) throws IloException {
            IloLinearNumExpr expr = cplex.linearNumExpr();
            int visited = 0;

            for (int i = 1; i <= ins.n; i++) {
                int vBar = point.prevVisit[i][t];
                if (vBar >= 0) {
                    for (int v = ins.pi[i][t]; v <= vBar; v++) {
                        if (lambda[i][v][t] != null) {
                            expr.addTerm(1.0, lambda[i][v][t]);
                        }
                    }
                    visited++;
                }
            }
            if (visited == 0) {
                return false;
            }

            feasCutCount++;
            cplex.addLe(expr, visited - 1.0, "LBBD_StrongFeasCut_" + feasCutCount + "_t" + t);
            return true;
        }

        void addNoGoodCut(MasterPoint point) throws IloException {
            addNoGoodCutRange(point);
        }

        IloRange addNoGoodCutRange(MasterPoint point) throws IloException {
            IloLinearNumExpr expr = cplex.linearNumExpr();
            for (int idx = 0; idx < routingLambdaRefs.size(); idx++) {
                LambdaRef ref = routingLambdaRefs.get(idx);
                if (point.lambdaOne[ref.i][ref.v][ref.t]) {
                    expr.addTerm(-1.0, ref.var);
                } else {
                    expr.addTerm(1.0, ref.var);
                }
            }
            double rhs = 1.0 - point.routeLambdaOneCount;
            noGoodCutCount++;
            return cplex.addGe(expr, rhs, "LBBD_NoGoodCut_" + noGoodCutCount);
        }

        IloRange addOptimalityCutRange(MasterPoint point, double phi, double lbR) throws IloException {
            double slope = phi - lbR;
            if (slope <= CUT_EPS) {
                return null;
            }
            IloLinearNumExpr expr = cplex.linearNumExpr();
            for (int t = 1; t <= ins.l; t++) {
                expr.addTerm(1.0, omega[t]);
            }
            for (int idx = 0; idx < routingLambdaRefs.size(); idx++) {
                LambdaRef ref = routingLambdaRefs.get(idx);
                if (point.lambdaOne[ref.i][ref.v][ref.t]) {
                    expr.addTerm(-slope, ref.var);
                } else {
                    expr.addTerm(slope, ref.var);
                }
            }
            double rhs = phi - slope * point.routeLambdaOneCount;
            optCutCount++;
            return cplex.addGe(expr, rhs, "LBBD_OptCut_" + optCutCount);
        }

        void addPeriodOptimalityCuts(MasterPoint point, double[] phiByPeriod, boolean[] skipIfTrue) throws IloException {
            for (int t = 1; t <= ins.l; t++) {
                double phiT = phiByPeriod[t];
                if (phiT <= CUT_EPS) {
                    continue;
                }
                if (phiT <= point.omega[t] + CUT_EPS) {
                    continue;
                }
                if (skipIfTrue != null && t < skipIfTrue.length && skipIfTrue[t]) {
                    continue;
                }
                addSinglePeriodOptimalityCut(point, phiT, t);
            }
        }

        IloRange addSinglePeriodOptimalityCut(MasterPoint point, double phiT, int t) throws IloException {
            if (phiT <= CUT_EPS) {
                return null;
            }
            if (phiT <= point.omega[t] + CUT_EPS) {
                return null;
            }

            IloLinearNumExpr expr = cplex.linearNumExpr();
            expr.addTerm(1.0, omega[t]);

            int onesCount = 0;
            for (int idx = 0; idx < routingLambdaRefs.size(); idx++) {
                LambdaRef ref = routingLambdaRefs.get(idx);
                if (ref.t != t) {
                    continue;
                }
                if (point.lambdaOne[ref.i][ref.v][ref.t]) {
                    expr.addTerm(-phiT, ref.var);
                    onesCount++;
                } else {
                    expr.addTerm(phiT, ref.var);
                }
            }

            double rhs = phiT - phiT * onesCount;
            optCutCount++;
            return cplex.addGe(expr, rhs, "LBBD_PeriodOptCut_" + optCutCount + "_t" + t);
        }

        IloRange addDualInitialCut(int t, double[][] dualW, double dualU0, String tag) throws IloException {
            IloLinearNumExpr expr = cplex.linearNumExpr();
            expr.addTerm(1.0, omega[t]);
            double rhs = ins.K * dualU0;

            for (int i = 1; i <= ins.n; i++) {
                for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                    IloNumVar var = lambda[i][v][t];
                    if (var == null) {
                        continue;
                    }
                    double w = 0.0;
                    if (dualW != null && i < dualW.length && v < dualW[i].length) {
                        w = dualW[i][v];
                    }
                    if (Math.abs(w) > 1e-9) {
                        expr.addTerm(-w, var);
                    }
                }
            }
            return cplex.addGe(expr, rhs, tag);
        }

        IloRange addPeriodRmpDualCut(int t, MasterPoint point, PeriodRouteMasterLpSolver.Result cut) throws IloException {
            IloLinearNumExpr expr = cplex.linearNumExpr();
            expr.addTerm(1.0, omega[t]);
            double rhs = ins.K * cut.dualU0;

            for (int i = 1; i <= ins.n; i++) {
                for (int v = ins.pi[i][t]; v <= t - 1; v++) {
                    IloNumVar var = lambda[i][v][t];
                    if (var == null) {
                        continue;
                    }
                    double ui = (cut.dualUiByPrev != null
                            && i < cut.dualUiByPrev.length
                            && v < cut.dualUiByPrev[i].length)
                            ? cut.dualUiByPrev[i][v]
                            : 0.0;
                    if (Math.abs(ui) <= 1e-9) {
                        continue;
                    }
                    expr.addTerm(-ui, var);
                }
            }

            optCutCount++;
            return cplex.addGe(expr, rhs, "LBBD_RmpDualCut_" + optCutCount + "_t" + t);
        }

        void deleteCut(IloRange cutRange) throws IloException {
            if (cutRange != null) {
                cplex.delete(cutRange);
            }
        }

        void setMIPStart(MasterPoint point) throws IloException {
            ArrayList<IloNumVar> vars = new ArrayList<IloNumVar>();
            ArrayList<Double> vals = new ArrayList<Double>();
            for (int idx = 0; idx < routingLambdaRefs.size(); idx++) {
                LambdaRef ref = routingLambdaRefs.get(idx);
                vars.add(ref.var);
                vals.add(point.lambdaOne[ref.i][ref.v][ref.t] ? 1.0 : 0.0);
            }
            for (int t = 1; t <= ins.l; t++) {
                vars.add(omega[t]);
                vals.add(point.omega[t]);
            }
            IloNumVar[] varArr = vars.toArray(new IloNumVar[0]);
            double[] valArr = new double[vals.size()];
            for (int idx = 0; idx < vals.size(); idx++) {
                valArr[idx] = vals.get(idx);
            }
            cplex.addMIPStart(varArr, valArr);
        }

    }
}
