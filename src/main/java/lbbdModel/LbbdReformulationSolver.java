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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public final class LbbdReformulationSolver {

    private static final double CUT_EPS = 1e-6;
    private static final double RESULT_TOL = CplexConfig.MIP_GAP;
    private static final double INCUMBENT_PRUNE_TOL = CplexConfig.MIP_GAP;
    private static final int MAX_ITERATIONS = 5000;
    private static final int MAX_RMP_PERIODS_PER_ITER = 3;
    private static final boolean ENABLE_INCUMBENT_PRUNE_NOGOOD = false;
    private static final boolean ENABLE_GLOBAL_COMBINATION_OPT_CUT = true;
    private static final boolean ENABLE_ITER_TIMING_LOG = true;
    private static final boolean LOG_TO_CONSOLE = false;
    private static final boolean ENABLE_HINT_SHORTCUT = false;
    private static final int MASTER_THREADS = 0;
    private static final int CUT_POOL_GRACE_ITERS = 3;
    private static final double CUT_ACTIVITY_TOL = 1e-5;
    private static final long NANOS_PER_SECOND = 1_000_000_000L;

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
    private final boolean enableInitialCuts;

    public LbbdReformulationSolver(Instance ins) {
        this(ins, true);
    }

    public LbbdReformulationSolver(Instance ins, boolean enableInitialCuts) {
        this.ins = ins;
        this.enableInitialCuts = enableInitialCuts;

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

    private abstract static class ActiveCutBase {
        final String key;
        final IloRange range;
        final int bornIteration;
        int lastTouchedIteration;
        int hitCount;
        double lastViolation;
        double maxViolation;

        ActiveCutBase(String key, IloRange range, int bornIteration, double initialViolation) {
            this.key = key;
            this.range = range;
            this.bornIteration = bornIteration;
            this.lastTouchedIteration = bornIteration;
            this.hitCount = 0;
            this.lastViolation = 0.0;
            this.maxViolation = 0.0;
            noteObservation(bornIteration, initialViolation);
        }

        void noteObservation(int iteration, double violation) {
            double positiveViolation = Math.max(0.0, violation);
            this.lastViolation = positiveViolation;
            if (positiveViolation > CUT_ACTIVITY_TOL) {
                this.lastTouchedIteration = iteration;
                this.hitCount++;
                if (positiveViolation > this.maxViolation) {
                    this.maxViolation = positiveViolation;
                }
            }
        }
    }

    private static final class ActiveRmpDualCut extends ActiveCutBase {
        final int period;
        final PeriodRouteMasterLpSolver.Result cut;

        ActiveRmpDualCut(String key, IloRange range, int bornIteration, double initialViolation, int period, PeriodRouteMasterLpSolver.Result cut) {
            super(key, range, bornIteration, initialViolation);
            this.period = period;
            this.cut = cut;
        }
    }

    private static final class ActivePeriodOptCut extends ActiveCutBase {
        final int period;
        final double phiT;
        final int visitedCount;
        final int[] prevVisitBySupplier;

        ActivePeriodOptCut(
                String key,
                IloRange range,
                int bornIteration,
                double initialViolation,
                int period,
                double phiT,
                int visitedCount,
                int[] prevVisitBySupplier
        ) {
            super(key, range, bornIteration, initialViolation);
            this.period = period;
            this.phiT = phiT;
            this.visitedCount = visitedCount;
            this.prevVisitBySupplier = prevVisitBySupplier;
        }
    }

    private static final class ActiveGlobalOptCut extends ActiveCutBase {
        final double slope;
        final double phi;
        final int routeLambdaOneCount;
        final int[][] prevVisit;

        ActiveGlobalOptCut(
                String key,
                IloRange range,
                int bornIteration,
                double initialViolation,
                double slope,
                double phi,
                int routeLambdaOneCount,
                int[][] prevVisit
        ) {
            super(key, range, bornIteration, initialViolation);
            this.slope = slope;
            this.phi = phi;
            this.routeLambdaOneCount = routeLambdaOneCount;
            this.prevVisit = prevVisit;
        }
    }

    private static final class ActiveIncumbentPruneCut extends ActiveCutBase {
        final int routeLambdaOneCount;
        final int[][] prevVisit;

        ActiveIncumbentPruneCut(
                String key,
                IloRange range,
                int bornIteration,
                double initialViolation,
                int routeLambdaOneCount,
                int[][] prevVisit
        ) {
            super(key, range, bornIteration, initialViolation);
            this.routeLambdaOneCount = routeLambdaOneCount;
            this.prevVisit = prevVisit;
        }
    }

    public SolveResult solve(double targetObjective, double targetTol) {
        return solve(targetObjective, targetTol, null);
    }

    public SolveResult solve(double targetObjective, double targetTol, int[][] prevVisitHint) {
        long startNs = System.nanoTime();
        long deadlineNs = startNs + Math.max(1L, Math.round(CplexConfig.TIME_LIMIT_SEC * NANOS_PER_SECOND));
        int iteration = 0;
        int feasibilityCuts = 0;
        int optimalityCuts = 0;
        int lpDualCuts = 0;
        int lpDualCutSkips = 0;
        int targetPruneCuts = 0;
        int incumbentPruneCuts = 0;
        int globalOptCuts = 0;
        double bestUpperBound = Double.POSITIVE_INFINITY;
        double bestLowerBound = Double.NEGATIVE_INFINITY;
        boolean allSubproblemsOptimal = true;
        if (ENABLE_HINT_SHORTCUT && !Double.isNaN(targetObjective) && prevVisitHint != null) {
            SolveResult fast = tryValidateByHint(ins, targetObjective, targetTol, prevVisitHint, startNs, deadlineNs);
            if (fast != null) {
                return fast;
            }
        }

        HashMap<String, PeriodCvrpSubproblemSolver.Result> subproblemCache =
                new HashMap<String, PeriodCvrpSubproblemSolver.Result>();
        HashMap<String, PeriodRouteMasterLpSolver.Result> rmpCutCache =
                new HashMap<String, PeriodRouteMasterLpSolver.Result>();
        HashSet<String> addedRmpCutKeys = new HashSet<String>();
        ArrayList<ActiveRmpDualCut> activeRmpDualCuts = new ArrayList<ActiveRmpDualCut>();
        HashSet<String> addedPeriodOptCutKeys = new HashSet<String>();
        ArrayList<ActivePeriodOptCut> activePeriodOptCuts = new ArrayList<ActivePeriodOptCut>();
        HashSet<String> addedGlobalOptCutKeys = new HashSet<String>();
        ArrayList<ActiveGlobalOptCut> activeGlobalOptCuts = new ArrayList<ActiveGlobalOptCut>();
        HashSet<String> activeIncumbentPruneCutKeys = new HashSet<String>();
        ArrayList<ActiveIncumbentPruneCut> activeIncumbentPruneCuts = new ArrayList<ActiveIncumbentPruneCut>();

        int subThreads = Math.min(ins.l, Runtime.getRuntime().availableProcessors());
        final PeriodCvrpSubproblemSolver[] subSolvers = new PeriodCvrpSubproblemSolver[ins.l + 1];
        for (int tt = 1; tt <= ins.l; tt++) {
            subSolvers[tt] = new PeriodCvrpSubproblemSolver();
        }
        ExecutorService subExecutor = subThreads > 1 ? Executors.newFixedThreadPool(subThreads) : null;
        int rmpThreads = Math.min(ins.l, Math.max(1, Runtime.getRuntime().availableProcessors() / 2));
        ExecutorService rmpExecutor = rmpThreads > 1 ? Executors.newFixedThreadPool(rmpThreads) : null;
        PeriodRouteMasterLpSolver rmpLpSolver = new PeriodRouteMasterLpSolver(ins);

        try (IloCplex cplex = new IloCplex()) {
            configure(cplex);
            MasterModel master = new MasterModel(cplex, ins);
            if (isPastDeadline(deadlineNs)) {
                return buildGlobalTimeLimitResult(
                        startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                        incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                );
            }
            RoutingLowerBoundSolver lbRSolver = new RoutingLowerBoundSolver(ins);
            RoutingLowerBoundSolver.Result lbRResult = lbRSolver.solve(remainingSeconds(deadlineNs));
            if (isPastDeadline(deadlineNs)) {
                return buildGlobalTimeLimitResult(
                        startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                        incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                );
            }
            double lbR = (lbRResult.feasible && !Double.isNaN(lbRResult.lbR) && !Double.isInfinite(lbRResult.lbR))
                    ? Math.max(0.0, lbRResult.lbR)
                    : 0.0;
            System.out.println("[LBBD] LB_R="
                    + fmt(lbR)
                    + " status=" + lbRResult.status
                    + " iters=" + lbRResult.iterations
                    + " cols=" + lbRResult.generatedColumns);

            if (enableInitialCuts) {
                // === T1 Initial Cuts (speed.tex Section B.2, eq 4-9), parallel across periods ===
                int t1CutsAdded = 0;
                T1InitialCutSolver.T1Result[] t1Results = new T1InitialCutSolver.T1Result[ins.l + 1];
                int t1Threads = Math.min(ins.l, Math.min(Runtime.getRuntime().availableProcessors(), 4));
                ExecutorService t1Executor = t1Threads > 1 ? Executors.newFixedThreadPool(t1Threads) : null;
                @SuppressWarnings("unchecked")
                Future<T1InitialCutSolver.T1Result>[] t1Futures =
                        (t1Executor != null) ? new Future[ins.l + 1] : null;
                if (isPastDeadline(deadlineNs)) {
                    return buildGlobalTimeLimitResult(
                            startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                            incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                    );
                }
                if (t1Executor != null) {
                    try {
                        final double t1TimeLimitSec = remainingSeconds(deadlineNs);
                        for (int t = 1; t <= ins.l; t++) {
                            final int period = t;
                            t1Futures[t] = t1Executor.submit(new Callable<T1InitialCutSolver.T1Result>() {
                                @Override
                                public T1InitialCutSolver.T1Result call() {
                                    return new T1InitialCutSolver(ins).solve(period, t1TimeLimitSec);
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
                    if (isPastDeadline(deadlineNs)) {
                        return buildGlobalTimeLimitResult(
                                startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                        );
                    }
                    T1InitialCutSolver fallbackSolver = new T1InitialCutSolver(ins);
                    for (int t = 1; t <= ins.l; t++) {
                        if (t1Results[t] == null) {
                            t1Results[t] = fallbackSolver.solve(t, remainingSeconds(deadlineNs));
                        }
                    }
                } else {
                    T1InitialCutSolver t1Solver = new T1InitialCutSolver(ins);
                    for (int t = 1; t <= ins.l; t++) {
                        if (isPastDeadline(deadlineNs)) {
                            return buildGlobalTimeLimitResult(
                                    startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                    incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                            );
                        }
                        t1Results[t] = t1Solver.solve(t, remainingSeconds(deadlineNs));
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
                if (isPastDeadline(deadlineNs)) {
                    return buildGlobalTimeLimitResult(
                            startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                            incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                    );
                }
                T2InitialCutSolver t2Solver = new T2InitialCutSolver(ins);
                T2InitialCutSolver.T2Result t2 = t2Solver.solve(remainingSeconds(deadlineNs));
                if (isPastDeadline(deadlineNs)) {
                    return buildGlobalTimeLimitResult(
                            startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                            incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                    );
                }
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
            } else {
                System.out.println("[LBBD] initial cuts disabled: skip T1/T2 initial cuts");
            }

            long totalMasterSolveNs = 0L;
            long totalSubSolveNs = 0L;
            long totalRmpSolveNs = 0L;
            long totalSubCalls = 0L;
            long totalSubCacheMisses = 0L;
            long totalRmpCalls = 0L;
            long totalRmpCacheMisses = 0L;

            double bestUpperBoundGapRef = Double.NaN;

            boolean allMasterOptimal = true;
            boolean forceStrictMasterGap = false;
            MasterPoint prevPoint = null;

            while (iteration < MAX_ITERATIONS) {
                if (isPastDeadline(deadlineNs)) {
                    return buildGlobalTimeLimitResult(
                            startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                            incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                    );
                }
                iteration++;

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
                cplex.setParam(IloCplex.Param.TimeLimit, normalizedTimeLimitSec(remainingSeconds(deadlineNs)));
                boolean solved = cplex.solve();
                iterMasterSolveNs = System.nanoTime() - masterStartNs;
                totalMasterSolveNs += iterMasterSolveNs;
                if (isPastDeadline(deadlineNs)) {
                    return buildGlobalTimeLimitResult(
                            startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                            incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                    );
                }
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
                updateActiveCutActivity(point, iteration, activeRmpDualCuts, activePeriodOptCuts, activeGlobalOptCuts, activeIncumbentPruneCuts);
                trimActiveRmpDualCuts(master, activeRmpDualCuts, addedRmpCutKeys, activeRmpDualCutCap(iteration), iteration);
                trimActivePeriodOptCuts(master, activePeriodOptCuts, addedPeriodOptCutKeys, activePeriodOptCutCap(iteration), iteration);
                trimActiveGlobalOptCuts(master, activeGlobalOptCuts, addedGlobalOptCutKeys, activeGlobalOptCutCap(iteration), iteration);
                trimActiveIncPruneCuts(master, activeIncumbentPruneCuts, activeIncumbentPruneCutKeys, activeIncPruneCutCap(iteration), iteration);

                double phi = 0.0;
                double[] phiByPeriod = new double[ins.l + 1];
                String[] patternKeysByPeriod = new String[ins.l + 1];
                String[] subCacheKeysByPeriod = new String[ins.l + 1];
                String[] rmpCacheKeysByPeriod = new String[ins.l + 1];
                int infeasiblePeriod = -1;

                // Phase 1: check cache and submit parallel tasks for misses
                PeriodCvrpSubproblemSolver.Result[] subResults = new PeriodCvrpSubproblemSolver.Result[ins.l + 1];
                @SuppressWarnings("unchecked")
                Future<PeriodCvrpSubproblemSolver.Result>[] subFutures = (subExecutor != null)
                        ? new Future[ins.l + 1] : null;
                int cacheMissCount = 0;

                if (isPastDeadline(deadlineNs)) {
                    return buildGlobalTimeLimitResult(
                            startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                            incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                    );
                }
                for (int t = 1; t <= ins.l; t++) {
                    iterSubCalls++;
                    totalSubCalls++;
                    String patternKey = buildPeriodRoutingPatternKey(point, t, ins.n);
                    String key = buildSubproblemCacheKey(patternKey);
                    patternKeysByPeriod[t] = patternKey;
                    subCacheKeysByPeriod[t] = key;
                    rmpCacheKeysByPeriod[t] = buildRmpCutCacheKey(patternKey);
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
                            final double subTimeLimitSec = remainingSeconds(deadlineNs);
                            subFutures[t] = subExecutor.submit(new Callable<PeriodCvrpSubproblemSolver.Result>() {
                                @Override
                                public PeriodCvrpSubproblemSolver.Result call() {
                                    return subSolvers[period].solve(ins, period, qBarT, zBarT, subTimeLimitSec);
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
                                subResults[t] = subSolvers[t].solve(
                                        ins, t, point.qBar[t], point.zBar[t], remainingSeconds(deadlineNs)
                                );
                            }
                        }
                    }
                    long subElapsed = System.nanoTime() - subStartNs;
                    iterSubSolveNs += subElapsed;
                    totalSubSolveNs += subElapsed;
                    if (isPastDeadline(deadlineNs)) {
                        return buildGlobalTimeLimitResult(
                                startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                        );
                    }

                    for (int t = 1; t <= ins.l; t++) {
                        if (subCacheKeysByPeriod[t] != null && subResults[t] != null
                                && !subproblemCache.containsKey(subCacheKeysByPeriod[t])) {
                            PeriodCvrpSubproblemSolver.Result sub = subResults[t];
                            if (sub.optimal || sub.provenInfeasible) {
                                subproblemCache.put(subCacheKeysByPeriod[t], sub);
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
                    if (isPastDeadline(deadlineNs)) {
                        return buildGlobalTimeLimitResult(
                                startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                        );
                    }
                    continue;
                }
                if (unresolvedPeriod >= 0) {
                    if (isPastDeadline(deadlineNs)) {
                        return buildGlobalTimeLimitResult(
                                startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                        );
                    }
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
                    if (isPastDeadline(deadlineNs)) {
                        return buildGlobalTimeLimitResult(
                                startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                        );
                    }
                    periodCoveredByFreshRmpCut = new boolean[ins.l + 1];
                    int periodsToProcess = Math.min(MAX_RMP_PERIODS_PER_ITER, ins.l);
                    int[] selectedPeriods = selectBestRmpPeriods(point, phiByPeriod, rmpCacheKeysByPeriod, addedRmpCutKeys, periodsToProcess);
                    @SuppressWarnings("unchecked")
                    Future<PeriodRouteMasterLpSolver.Result>[] rmpFutures = (rmpExecutor != null)
                            ? new Future[ins.l + 1] : null;
                    int[][] prevVisitArgPeriods = new int[ins.l + 1][];
                    final PeriodRouteMasterLpSolver finalRmpLpSolver = rmpLpSolver;
                    final MasterPoint finalPoint = point;

                    long rmpStartNs = System.nanoTime();
                    for (int idx = 0; idx < selectedPeriods.length; idx++) {
                        int bestT = selectedPeriods[idx];
                        if (bestT < 0) {
                            continue;
                        }
                        String key = rmpCacheKeysByPeriod[bestT];
                        if (rmpCutCache.containsKey(key)) {
                            continue;
                        }
                        iterRmpCalls++;
                        totalRmpCalls++;
                        iterRmpCacheMisses++;
                        totalRmpCacheMisses++;
                        prevVisitArgPeriods[bestT] = buildPrevVisitForPeriod(point, bestT, ins.n);
                        if (rmpExecutor != null) {
                            final int period = bestT;
                            final int[] prevVisitForPeriod = prevVisitArgPeriods[bestT];
                            final double rmpTimeLimitSec = remainingSeconds(deadlineNs);
                            rmpFutures[bestT] = rmpExecutor.submit(new Callable<PeriodRouteMasterLpSolver.Result>() {
                                @Override
                                public PeriodRouteMasterLpSolver.Result call() {
                                    return finalRmpLpSolver.solve(
                                            ins,
                                            period,
                                            finalPoint.qBar[period],
                                            finalPoint.zBar[period],
                                            prevVisitForPeriod,
                                            rmpTimeLimitSec
                                    );
                                }
                            });
                        } else {
                            PeriodRouteMasterLpSolver.Result solvedCut = rmpLpSolver.solve(
                                    ins,
                                    bestT,
                                    point.qBar[bestT],
                                    point.zBar[bestT],
                                    prevVisitArgPeriods[bestT],
                                    remainingSeconds(deadlineNs)
                            );
                            if (solvedCut.feasible && solvedCut.optimal && solvedCut.pricingProvedOptimal) {
                                rmpCutCache.put(key, solvedCut);
                            }
                        }
                    }
                    if (rmpExecutor != null) {
                        for (int idx = 0; idx < selectedPeriods.length; idx++) {
                            int bestT = selectedPeriods[idx];
                            if (bestT < 0 || rmpFutures[bestT] == null) {
                                continue;
                            }
                            try {
                                PeriodRouteMasterLpSolver.Result solvedCut = rmpFutures[bestT].get();
                                if (solvedCut.feasible && solvedCut.optimal && solvedCut.pricingProvedOptimal) {
                                    rmpCutCache.put(rmpCacheKeysByPeriod[bestT], solvedCut);
                                }
                            } catch (Exception e) {
                                throw new RuntimeException("RMP dual-cut solve failed for t=" + bestT, e);
                            }
                        }
                    }
                    long rmpElapsed = System.nanoTime() - rmpStartNs;
                    iterRmpSolveNs += rmpElapsed;
                    totalRmpSolveNs += rmpElapsed;
                    if (isPastDeadline(deadlineNs)) {
                        return buildGlobalTimeLimitResult(
                                startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                        );
                    }

                    for (int idx = 0; idx < selectedPeriods.length; idx++) {
                        int bestT = selectedPeriods[idx];
                        if (bestT < 0) {
                            continue;
                        }
                        String key = rmpCacheKeysByPeriod[bestT];
                        PeriodRouteMasterLpSolver.Result rmpCut = rmpCutCache.get(key);
                        if (rmpCut == null
                                || !rmpCut.feasible
                                || !rmpCut.optimal
                                || !rmpCut.pricingProvedOptimal
                                || !rmpCut.artificialClean
                                || rmpCut.lpObjective <= point.omega[bestT] + CUT_EPS) {
                            continue;
                        }

                        double cutAtPoint = evaluateRmpDualCutAtPoint(point, bestT, rmpCut, ins);
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
                        activeRmpDualCuts.add(new ActiveRmpDualCut(
                                key,
                                addedRange,
                                iteration,
                                Math.max(0.0, cutAtPoint - point.omega[bestT]),
                                bestT,
                                rmpCut
                        ));
                        trimActiveRmpDualCuts(master, activeRmpDualCuts, addedRmpCutKeys, activeRmpDualCutCap(iteration), iteration);
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
                                String periodKey = buildPeriodOptCutKey(patternKeysByPeriod[t]);
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
                            String periodKey = buildPeriodOptCutKey(patternKeysByPeriod[bestT]);
                            IloRange addedRange = master.addSinglePeriodOptimalityCut(point, phiByPeriod[bestT], bestT);
                            if (addedRange == null) {
                                break;
                            }
                            addedPeriodOptCutKeys.add(periodKey);
                            activePeriodOptCuts.add(new ActivePeriodOptCut(
                                    periodKey,
                                    addedRange,
                                    iteration,
                                    Math.max(0.0, phiByPeriod[bestT] - point.omega[bestT]),
                                    bestT,
                                    phiByPeriod[bestT],
                                    countVisitedForPeriod(point, bestT, ins.n),
                                    copyPrevVisitForPeriod(point, bestT, ins.n)
                            ));
                            trimActivePeriodOptCuts(master, activePeriodOptCuts, addedPeriodOptCutKeys, activePeriodOptCutCap(iteration), iteration);
                            addedPeriodOptThisIter++;
                        }
                    }
                    if (ENABLE_GLOBAL_COMBINATION_OPT_CUT) {
                        if (iteration >= globalOptCutStartIter && delta >= globalOptCutMinDelta) {
                            String globalOptKey = buildGlobalOptCutKey(buildGlobalRoutingPatternKey(point, ins.n, ins.l));
                            if (!addedGlobalOptCutKeys.contains(globalOptKey)) {
                                IloRange range = master.addOptimalityCutRange(point, phi, lbR);
                                if (range != null) {
                                    addedGlobalOptCutKeys.add(globalOptKey);
                                    double slope = phi - lbR;
                                    activeGlobalOptCuts.add(new ActiveGlobalOptCut(
                                            globalOptKey,
                                            range,
                                            iteration,
                                            Math.max(0.0, delta),
                                            slope,
                                            phi,
                                            point.routeLambdaOneCount,
                                            copyPrevVisitMatrix(point, ins.n, ins.l)
                                    ));
                                    trimActiveGlobalOptCuts(master, activeGlobalOptCuts, addedGlobalOptCutKeys, activeGlobalOptCutCap(iteration), iteration);
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
                            String incPruneKey = buildIncumbentPruneCutKey(buildGlobalRoutingPatternKey(point, ins.n, ins.l));
                            if (!activeIncumbentPruneCutKeys.contains(incPruneKey)) {
                                IloRange ng = master.addNoGoodCutRange(point);
                                activeIncumbentPruneCutKeys.add(incPruneKey);
                                activeIncumbentPruneCuts.add(new ActiveIncumbentPruneCut(
                                        incPruneKey,
                                        ng,
                                        iteration,
                                        1.0,
                                        point.routeLambdaOneCount,
                                        copyPrevVisitMatrix(point, ins.n, ins.l)
                                ));
                                trimActiveIncPruneCuts(master, activeIncumbentPruneCuts, activeIncumbentPruneCutKeys, activeIncPruneCutCap(iteration), iteration);
                                incumbentPruneCuts++;
                            }
                        }
                        if (isPastDeadline(deadlineNs)) {
                            return buildGlobalTimeLimitResult(
                                    startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                    incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                            );
                        }
                        continue;
                    }
                }

                if (!Double.isNaN(targetObjective) && exactObjective > targetObjective + targetTol) {
                    master.addNoGoodCut(point);
                    targetPruneCuts++;
                    forceStrictMasterGap = false;
                    if (isPastDeadline(deadlineNs)) {
                        return buildGlobalTimeLimitResult(
                                startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                        );
                    }
                    continue;
                }

                if (violatedOptimality) {
                    forceStrictMasterGap = false;
                    if (isPastDeadline(deadlineNs)) {
                        return buildGlobalTimeLimitResult(
                                startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                        );
                    }
                    continue;
                }

                if (masterMipGapThisIter > CplexConfig.MIP_GAP + 1e-12) {
                    forceStrictMasterGap = true;
                    if (isPastDeadline(deadlineNs)) {
                        return buildGlobalTimeLimitResult(
                                startNs, iteration, feasibilityCuts, optimalityCuts, lpDualCuts,
                                incumbentPruneCuts, targetPruneCuts, bestUpperBound, bestLowerBound
                        );
                    }
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
            if (rmpExecutor != null) {
                rmpExecutor.shutdownNow();
            }
            rmpLpSolver.close();
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
            long startNs,
            long deadlineNs
    ) {
        if (prevVisitHint.length < ins.n + 1) {
            return null;
        }
        try (PeriodCvrpSubproblemSolver subproblemSolver = new PeriodCvrpSubproblemSolver()) {
            double phi = 0.0;
            for (int t = 1; t <= ins.l; t++) {
                if (isPastDeadline(deadlineNs)) {
                    return buildGlobalTimeLimitResult(
                            startNs, 0, 0, 0, 0, 0, 0, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY
                    );
                }
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
                PeriodCvrpSubproblemSolver.Result sub = subproblemSolver.solve(
                        ins, t, qBar, zBar, remainingSeconds(deadlineNs)
                );
                if (!sub.feasible || !sub.optimal) {
                    if (isPastDeadline(deadlineNs)) {
                        return buildGlobalTimeLimitResult(
                                startNs, 0, 0, 0, 0, 0, 0, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY
                        );
                    }
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

    private static String buildPeriodRoutingPatternKey(MasterPoint point, int t, int n) {
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

    private static String buildSubproblemCacheKey(String periodPatternKey) {
        return periodPatternKey == null ? null : ("SUB|" + periodPatternKey);
    }

    private static String buildRmpCutCacheKey(String periodPatternKey) {
        return periodPatternKey == null ? null : ("RMP|" + periodPatternKey);
    }

    private static int[] buildPrevVisitForPeriod(MasterPoint point, int t, int n) {
        int[] prev = new int[n + 1];
        for (int i = 1; i <= n; i++) {
            prev[i] = point.prevVisit[i][t];
        }
        return prev;
    }

    private static int[] copyPrevVisitForPeriod(MasterPoint point, int t, int n) {
        return buildPrevVisitForPeriod(point, t, n);
    }

    private static int[][] copyPrevVisitMatrix(MasterPoint point, int n, int l) {
        int[][] copy = new int[n + 1][l + 1];
        for (int i = 1; i <= n; i++) {
            for (int t = 1; t <= l; t++) {
                copy[i][t] = point.prevVisit[i][t];
            }
        }
        return copy;
    }

    private static int countVisitedForPeriod(MasterPoint point, int t, int n) {
        int count = 0;
        for (int i = 1; i <= n; i++) {
            if (point.prevVisit[i][t] >= 0) {
                count++;
            }
        }
        return count;
    }

    private static String buildPeriodOptCutKey(String periodPatternKey) {
        if (periodPatternKey == null) {
            return null;
        }
        return "POPT|" + periodPatternKey;
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

    private static String buildGlobalOptCutKey(String globalPatternKey) {
        return "GOPT|" + globalPatternKey;
    }

    private static String buildIncumbentPruneCutKey(String globalPatternKey) {
        return "NOGOOD|" + globalPatternKey;
    }

    private static int[] selectBestRmpPeriods(
            MasterPoint point,
            double[] phiByPeriod,
            String[] rmpCacheKeysByPeriod,
            HashSet<String> addedRmpCutKeys,
            int limit
    ) {
        int[] selected = new int[limit];
        for (int idx = 0; idx < selected.length; idx++) {
            selected[idx] = -1;
        }
        boolean[] picked = new boolean[phiByPeriod.length];
        for (int pick = 0; pick < limit; pick++) {
            int bestT = -1;
            double bestGap = CUT_EPS;
            for (int t = 1; t < phiByPeriod.length; t++) {
                if (picked[t]) {
                    continue;
                }
                String key = rmpCacheKeysByPeriod[t];
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
            selected[pick] = bestT;
        }
        return selected;
    }

    private static double evaluateRmpDualCutAtPoint(
            MasterPoint point,
            int t,
            PeriodRouteMasterLpSolver.Result cut,
            Instance ins
    ) {
        double cutAtPoint = ins.K * cut.dualU0;
        for (int i = 1; i <= ins.n; i++) {
            int v = point.prevVisit[i][t];
            if (v < 0) {
                continue;
            }
            if (cut.dualUiByPrev != null
                    && i < cut.dualUiByPrev.length
                    && v < cut.dualUiByPrev[i].length) {
                cutAtPoint += cut.dualUiByPrev[i][v];
            }
        }
        return cutAtPoint;
    }

    private void updateActiveCutActivity(
            MasterPoint point,
            int iteration,
            ArrayList<ActiveRmpDualCut> activeRmpDualCuts,
            ArrayList<ActivePeriodOptCut> activePeriodOptCuts,
            ArrayList<ActiveGlobalOptCut> activeGlobalOptCuts,
            ArrayList<ActiveIncumbentPruneCut> activeIncumbentPruneCuts
    ) {
        for (int idx = 0; idx < activeRmpDualCuts.size(); idx++) {
            ActiveRmpDualCut cut = activeRmpDualCuts.get(idx);
            double violation = evaluateRmpDualCutAtPoint(point, cut.period, cut.cut, ins) - point.omega[cut.period];
            cut.noteObservation(iteration, violation);
        }
        for (int idx = 0; idx < activePeriodOptCuts.size(); idx++) {
            ActivePeriodOptCut cut = activePeriodOptCuts.get(idx);
            double lhs = point.omega[cut.period];
            for (int i = 1; i <= ins.n; i++) {
                int currentPrev = point.prevVisit[i][cut.period];
                if (currentPrev < 0) {
                    continue;
                }
                lhs += (currentPrev == cut.prevVisitBySupplier[i]) ? -cut.phiT : cut.phiT;
            }
            double rhs = cut.phiT - cut.phiT * cut.visitedCount;
            cut.noteObservation(iteration, rhs - lhs);
        }
        for (int idx = 0; idx < activeGlobalOptCuts.size(); idx++) {
            ActiveGlobalOptCut cut = activeGlobalOptCuts.get(idx);
            double lhs = point.omegaSum;
            for (int t = 1; t <= ins.l; t++) {
                for (int i = 1; i <= ins.n; i++) {
                    int currentPrev = point.prevVisit[i][t];
                    if (currentPrev < 0) {
                        continue;
                    }
                    lhs += (currentPrev == cut.prevVisit[i][t]) ? -cut.slope : cut.slope;
                }
            }
            double rhs = cut.phi - cut.slope * cut.routeLambdaOneCount;
            cut.noteObservation(iteration, rhs - lhs);
        }
        for (int idx = 0; idx < activeIncumbentPruneCuts.size(); idx++) {
            ActiveIncumbentPruneCut cut = activeIncumbentPruneCuts.get(idx);
            int matches = 0;
            int differingActives = 0;
            for (int t = 1; t <= ins.l; t++) {
                for (int i = 1; i <= ins.n; i++) {
                    int currentPrev = point.prevVisit[i][t];
                    if (currentPrev < 0) {
                        continue;
                    }
                    if (currentPrev == cut.prevVisit[i][t]) {
                        matches++;
                    } else {
                        differingActives++;
                    }
                }
            }
            double lhs = -matches + differingActives;
            double rhs = 1.0 - cut.routeLambdaOneCount;
            cut.noteObservation(iteration, rhs - lhs);
        }
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

    private static SolveResult buildGlobalTimeLimitResult(
            long startNs,
            int iteration,
            int feasibilityCuts,
            int optimalityCuts,
            int lpDualCuts,
            int incumbentPruneCuts,
            int targetPruneCuts,
            double bestUpperBound,
            double bestLowerBound
    ) {
        double sec = elapsedSec(startNs);
        double finalBestBound = isFinite(bestLowerBound) ? bestLowerBound : Double.NaN;
        if (isFinite(bestUpperBound)) {
            double gap = isFinite(finalBestBound) ? relativeGap(bestUpperBound, finalBestBound) : Double.NaN;
            String status = "LBBD_GlobalTimeLimit(iter=" + iteration
                    + ",feasCuts=" + feasibilityCuts
                    + ",optCuts=" + optimalityCuts
                    + ",lpDualCuts=" + lpDualCuts
                    + ",incPruneCuts=" + incumbentPruneCuts
                    + ",targetPruneCuts=" + targetPruneCuts + ")";
            return new SolveResult(
                    "LBBDReformulation",
                    true,
                    false,
                    status,
                    bestUpperBound,
                    finalBestBound,
                    gap,
                    sec
            );
        }
        String status = "LBBD_GlobalTimeLimit_NoIncumbent(iter=" + iteration
                + ",feasCuts=" + feasibilityCuts
                + ",optCuts=" + optimalityCuts
                + ",lpDualCuts=" + lpDualCuts
                + ",incPruneCuts=" + incumbentPruneCuts
                + ",targetPruneCuts=" + targetPruneCuts + ")";
        return new SolveResult(
                "LBBDReformulation",
                false,
                false,
                status,
                Double.NaN,
                finalBestBound,
                Double.NaN,
                sec
        );
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

    private static double normalizedTimeLimitSec(double timeLimitSec) {
        return Math.max(1e-3, timeLimitSec);
    }

    private static double elapsedSec(long startNs) {
        return (System.nanoTime() - startNs) / (double) NANOS_PER_SECOND;
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

    private static <T extends ActiveCutBase> void trimActiveCuts(
            MasterModel master,
            ArrayList<T> activeCuts,
            HashSet<String> keys,
            int cap,
            int iteration
    ) throws IloException {
        while (activeCuts.size() > cap) {
            int evictIndex = selectEvictionIndex(activeCuts, iteration);
            T evicted = activeCuts.remove(evictIndex);
            master.deleteCut(evicted.range);
            keys.remove(evicted.key);
        }
    }

    private static <T extends ActiveCutBase> int selectEvictionIndex(ArrayList<T> activeCuts, int iteration) {
        int bestIndex = 0;
        for (int idx = 1; idx < activeCuts.size(); idx++) {
            if (compareEvictionPriority(activeCuts.get(idx), activeCuts.get(bestIndex), iteration) < 0) {
                bestIndex = idx;
            }
        }
        return bestIndex;
    }

    private static int compareEvictionPriority(ActiveCutBase a, ActiveCutBase b, int iteration) {
        boolean aYoung = (iteration - a.bornIteration) < CUT_POOL_GRACE_ITERS;
        boolean bYoung = (iteration - b.bornIteration) < CUT_POOL_GRACE_ITERS;
        if (aYoung != bYoung) {
            return aYoung ? 1 : -1;
        }
        boolean aActive = a.lastViolation > CUT_ACTIVITY_TOL;
        boolean bActive = b.lastViolation > CUT_ACTIVITY_TOL;
        if (aActive != bActive) {
            return aActive ? 1 : -1;
        }
        if (a.hitCount != b.hitCount) {
            return Integer.compare(a.hitCount, b.hitCount);
        }
        if (a.lastTouchedIteration != b.lastTouchedIteration) {
            return Integer.compare(a.lastTouchedIteration, b.lastTouchedIteration);
        }
        if (Math.abs(a.maxViolation - b.maxViolation) > 1e-9) {
            return Double.compare(a.maxViolation, b.maxViolation);
        }
        if (a.bornIteration != b.bornIteration) {
            return Integer.compare(a.bornIteration, b.bornIteration);
        }
        return a.key.compareTo(b.key);
    }

    private static void trimActiveRmpDualCuts(
            MasterModel master,
            ArrayList<ActiveRmpDualCut> activeCuts,
            HashSet<String> keys,
            int cap,
            int iteration
    ) throws IloException {
        trimActiveCuts(master, activeCuts, keys, cap, iteration);
    }

    private static void trimActivePeriodOptCuts(
            MasterModel master,
            ArrayList<ActivePeriodOptCut> activeCuts,
            HashSet<String> keys,
            int cap,
            int iteration
    ) throws IloException {
        trimActiveCuts(master, activeCuts, keys, cap, iteration);
    }

    private static void trimActiveGlobalOptCuts(
            MasterModel master,
            ArrayList<ActiveGlobalOptCut> activeCuts,
            HashSet<String> keys,
            int cap,
            int iteration
    ) throws IloException {
        trimActiveCuts(master, activeCuts, keys, cap, iteration);
    }

    private static void trimActiveIncPruneCuts(
            MasterModel master,
            ArrayList<ActiveIncumbentPruneCut> activeCuts,
            HashSet<String> keys,
            int cap,
            int iteration
    ) throws IloException {
        trimActiveCuts(master, activeCuts, keys, cap, iteration);
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
            int mipStartCount = cplex.getNMIPStarts();
            if (mipStartCount > 0) {
                cplex.deleteMIPStarts(0, mipStartCount);
            }
            cplex.addMIPStart(varArr, valArr);
        }

    }
}
