import instance.Instance;
import lbbdModel.LbbdDetailedResult;
import lbbdModel.LbbdReformulationSolver;
import lbbdModel.LbbdSolutionAnalytics;
import model.CplexConfig;
import model.SolveResult;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class ExperimentRunner {
    private static final double DEFAULT_TIME_LIMIT_SEC = 7200.0;
    private static final String[] ABLATION_INSTANCES = new String[]{
            "MVPRP1_10_6_3", "MVPRP2_10_6_3", "MVPRP3_10_6_3", "MVPRP4_10_6_3",
            "MVPRP1_10_9_2", "MVPRP2_10_9_2", "MVPRP3_10_9_2", "MVPRP4_10_9_2",
            "MVPRP1_10_9_3", "MVPRP2_10_9_3", "MVPRP3_10_9_3", "MVPRP4_10_9_3",
            "MVPRP1_15_3_3", "MVPRP2_15_3_3", "MVPRP3_15_3_3", "MVPRP4_15_3_3",
    };
    private static final String[] SENSITIVITY_BASE_DEFAULT = new String[]{
            "MVPRP1_10_6_2", "MVPRP3_10_6_2",
    };
    private static final String[] SENSITIVITY_BASE_FALLBACK = new String[]{
            "MVPRP1_10_3_2", "MVPRP3_10_3_2",
    };
    private static final double[] Q_RATIOS = new double[]{0.80, 1.00, 1.20, 1.50};
    private static final int[] K_VALUES = new int[]{1, 2, 3, 4};
    private static final double[] F_RATIOS = new double[]{0.50, 1.00, 2.00, 5.00};
    private static final double[] H_RATIOS = new double[]{0.50, 1.00, 2.00, 5.00};
    private static final Path ABLATION_OUTPUT = Paths.get("ablation_results.csv");
    private static final Path SENSITIVITY_Q_OUTPUT = Paths.get("sensitivity_Q.csv");
    private static final Path SENSITIVITY_K_OUTPUT = Paths.get("sensitivity_K.csv");
    private static final Path SENSITIVITY_F_OUTPUT = Paths.get("sensitivity_f.csv");
    private static final Path SENSITIVITY_H_OUTPUT = Paths.get("sensitivity_h.csv");
    private static final String ABLATION_HEADER =
            "instance,variant,feasible,optimal,objective,best_bound,gap,time_sec,iterations,feasCuts,optCuts";
    private static final String SENSITIVITY_Q_HEADER =
            "instance,Q_ratio,Q_value,feasible,optimal,objective,time_sec,route_cost,inventory_cost,production_cost,avg_visit_frequency,avg_pickup_per_visit";
    private static final String SENSITIVITY_K_HEADER =
            "instance,K_value,feasible,optimal,objective,time_sec,route_cost,inventory_cost,production_cost,avg_visit_frequency,avg_pickup_per_visit";
    private static final String SENSITIVITY_F_HEADER =
            "instance,f_ratio,f_value,feasible,optimal,objective,time_sec,route_cost,inventory_cost,production_cost,num_production_periods,avg_production_batch,raw_inv_avg,finished_inv_avg";
    private static final String SENSITIVITY_H_HEADER =
            "instance,h_ratio,feasible,optimal,objective,time_sec,route_cost,inventory_cost,production_cost,supplier_inv_avg,raw_inv_avg,finished_inv_avg,avg_visit_frequency";
    private static final Path[] FULL_RESULT_SOURCES = new Path[]{
            Paths.get("10_results_24cases_600s.csv"),
            Paths.get("15_results_24cases_600s.csv"),
    };
    private static final OutputStream DEV_NULL = new OutputStream() {
        @Override
        public void write(int b) {
            // discard
        }
    };
    private static final Pattern ITER_PATTERN = Pattern.compile("iter=(\\d+)");
    private static final Pattern FEAS_CUT_PATTERN = Pattern.compile("feasCuts=(\\d+)");
    private static final Pattern OPT_CUT_PATTERN = Pattern.compile("optCuts=(\\d+)");

    private ExperimentRunner() {
    }

    public static void main(String[] args) throws Exception {
        CliOptions options = CliOptions.parse(args);
        if ("ablation".equals(options.mode)) {
            runAblation(options);
            return;
        }
        if ("sensitivity_q".equals(options.mode)) {
            runSensitivityQ(options);
            return;
        }
        if ("sensitivity_k".equals(options.mode)) {
            runSensitivityK(options);
            return;
        }
        if ("sensitivity_f".equals(options.mode)) {
            runSensitivityF(options);
            return;
        }
        if ("sensitivity_h".equals(options.mode)) {
            runSensitivityH(options);
            return;
        }
        if ("all".equals(options.mode)) {
            runAblation(options);
            runSensitivityQ(options);
            runSensitivityK(options);
            runSensitivityF(options);
            runSensitivityH(options);
            return;
        }
        throw new IllegalArgumentException("Unsupported mode: " + options.mode);
    }

    private static void runAblation(CliOptions options) throws Exception {
        String[] selectedInstances = options.hasExplicitInstances()
                ? options.explicitInstances
                : toInstancePaths(ABLATION_INSTANCES);
        VariantSpec[] variants = new VariantSpec[]{
                new VariantSpec("base", new LbbdReformulationSolver.Config(false, false, false, false, options.timeLimitSec)),
                new VariantSpec("+cuts", new LbbdReformulationSolver.Config(false, true, false, false, options.timeLimitSec)),
                new VariantSpec("+vi", new LbbdReformulationSolver.Config(false, true, true, false, options.timeLimitSec)),
        };
        HashMap<String, FullResultRow> fullResults = loadFullResults();
        ExistingOutput existing = loadExistingOutput(ABLATION_OUTPUT, ABLATION_HEADER, new int[]{0, 1}, 3);
        ensureOutputFile(ABLATION_OUTPUT, ABLATION_HEADER);

        int total = selectedInstances.length * (variants.length + 1);
        int progress = countCompletedTargets(existing, buildAblationTargetKeys(selectedInstances, variants));
        try (BufferedWriter writer = appendWriter(ABLATION_OUTPUT)) {
            for (int idx = 0; idx < selectedInstances.length; idx++) {
                String instancePath = normalizeInstancePath(selectedInstances[idx]);
                String instanceName = instanceName(instancePath);
                Instance baseInstance = null;
                for (int v = 0; v < variants.length; v++) {
                    VariantSpec variant = variants[v];
                    String key = buildKey(instanceName, variant.name);
                    if (existing.rowsByKey.containsKey(key)) {
                        System.out.println("[" + progress + "/" + total + "] skip: ablation " + instanceName + " " + variant.name);
                        continue;
                    }
                    if (baseInstance == null) {
                        baseInstance = loadInstance(instancePath);
                    }
                    SolveWithDetails run = solveQuietly(baseInstance, variant.config);
                    String[] row = new String[]{
                            instanceName,
                            variant.name,
                            fmtBoolean(run.summary.feasible),
                            fmtBoolean(run.summary.optimal),
                            fmt(run.summary.objective),
                            fmt(run.summary.bestBound),
                            fmt(run.summary.mipGap),
                            fmt(run.summary.solveTimeSec),
                            Integer.toString(run.iterations()),
                            Integer.toString(run.feasibilityCuts()),
                            Integer.toString(run.optimalityCuts()),
                    };
                    writeCsvRow(writer, row);
                    existing.rowsByKey.put(key, new ExistingRow(row, 3));
                    progress++;
                    System.out.println("[" + progress + "/" + total + "] ablation " + instanceName + " " + variant.name);
                }

                String fullKey = buildKey(instanceName, "full");
                if (existing.rowsByKey.containsKey(fullKey)) {
                    System.out.println("[" + progress + "/" + total + "] skip: ablation " + instanceName + " full(imported)");
                    continue;
                }
                FullResultRow full = fullResults.get(instanceName);
                if (full == null) {
                    throw new IllegalStateException("Missing full=lbbd_no_init row for instance " + instanceName);
                }
                String[] row = new String[]{
                        instanceName,
                        "full",
                        fmtBoolean(full.feasible),
                        fmtBoolean(full.optimal),
                        fmt(full.objective),
                        fmt(full.bestBound),
                        fmt(full.gap),
                        fmt(full.timeSec),
                        Integer.toString(full.iterations),
                        Integer.toString(full.feasibilityCuts),
                        Integer.toString(full.optimalityCuts),
                };
                writeCsvRow(writer, row);
                existing.rowsByKey.put(fullKey, new ExistingRow(row, 3));
                progress++;
                System.out.println("[" + progress + "/" + total + "] ablation " + instanceName + " full(imported)");
            }
        }
    }

    private static void runSensitivityQ(CliOptions options) throws Exception {
        ExistingOutput existing = loadExistingOutput(SENSITIVITY_Q_OUTPUT, SENSITIVITY_Q_HEADER, new int[]{0, 1}, 4);
        if (options.hasExplicitInstances()) {
            evaluateSensitivityQ(primarySensitivityInstances(options), options.timeLimitSec, existing);
            return;
        }

        ExistingFamily family = detectSensitivityFamily(existing, SENSITIVITY_BASE_DEFAULT, SENSITIVITY_BASE_FALLBACK);
        if (family == ExistingFamily.FALLBACK) {
            validateExistingKeys(existing, buildSensitivityQTargetKeys(toInstancePaths(SENSITIVITY_BASE_FALLBACK)), SENSITIVITY_Q_OUTPUT);
            evaluateSensitivityQ(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, existing);
            return;
        }

        ArrayList<String> primaryKeys = buildSensitivityQTargetKeys(toInstancePaths(SENSITIVITY_BASE_DEFAULT));
        if (family == ExistingFamily.PRIMARY) {
            validateExistingKeys(existing, primaryKeys, SENSITIVITY_Q_OUTPUT);
            if (containsNonOptimal(existing, primaryKeys)) {
                resetOutputFile(SENSITIVITY_Q_OUTPUT, SENSITIVITY_Q_HEADER);
                evaluateSensitivityQ(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, emptyExistingOutput());
                return;
            }
        }

        FamilyRun primary = evaluateSensitivityQ(toInstancePaths(SENSITIVITY_BASE_DEFAULT), options.timeLimitSec, existing);
        if (!primary.allOptimal) {
            resetOutputFile(SENSITIVITY_Q_OUTPUT, SENSITIVITY_Q_HEADER);
            evaluateSensitivityQ(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, emptyExistingOutput());
        }
    }

    private static void runSensitivityK(CliOptions options) throws Exception {
        ExistingOutput existing = loadExistingOutput(SENSITIVITY_K_OUTPUT, SENSITIVITY_K_HEADER, new int[]{0, 1}, 3);
        if (options.hasExplicitInstances()) {
            evaluateSensitivityK(primarySensitivityInstances(options), options.timeLimitSec, existing);
            return;
        }

        ExistingFamily family = detectSensitivityFamily(existing, SENSITIVITY_BASE_DEFAULT, SENSITIVITY_BASE_FALLBACK);
        if (family == ExistingFamily.FALLBACK) {
            validateExistingKeys(existing, buildSensitivityKTargetKeys(toInstancePaths(SENSITIVITY_BASE_FALLBACK)), SENSITIVITY_K_OUTPUT);
            evaluateSensitivityK(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, existing);
            return;
        }

        ArrayList<String> primaryKeys = buildSensitivityKTargetKeys(toInstancePaths(SENSITIVITY_BASE_DEFAULT));
        if (family == ExistingFamily.PRIMARY) {
            validateExistingKeys(existing, primaryKeys, SENSITIVITY_K_OUTPUT);
            if (containsNonOptimal(existing, primaryKeys)) {
                resetOutputFile(SENSITIVITY_K_OUTPUT, SENSITIVITY_K_HEADER);
                evaluateSensitivityK(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, emptyExistingOutput());
                return;
            }
        }

        FamilyRun primary = evaluateSensitivityK(toInstancePaths(SENSITIVITY_BASE_DEFAULT), options.timeLimitSec, existing);
        if (!primary.allOptimal) {
            resetOutputFile(SENSITIVITY_K_OUTPUT, SENSITIVITY_K_HEADER);
            evaluateSensitivityK(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, emptyExistingOutput());
        }
    }

    private static void runSensitivityF(CliOptions options) throws Exception {
        ExistingOutput existing = loadExistingOutput(SENSITIVITY_F_OUTPUT, SENSITIVITY_F_HEADER, new int[]{0, 1}, 4);
        if (options.hasExplicitInstances()) {
            evaluateSensitivityF(primarySensitivityInstances(options), options.timeLimitSec, existing);
            return;
        }

        ExistingFamily family = detectSensitivityFamily(existing, SENSITIVITY_BASE_DEFAULT, SENSITIVITY_BASE_FALLBACK);
        if (family == ExistingFamily.FALLBACK) {
            validateExistingKeys(existing, buildSensitivityFTargetKeys(toInstancePaths(SENSITIVITY_BASE_FALLBACK)), SENSITIVITY_F_OUTPUT);
            evaluateSensitivityF(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, existing);
            return;
        }

        ArrayList<String> primaryKeys = buildSensitivityFTargetKeys(toInstancePaths(SENSITIVITY_BASE_DEFAULT));
        if (family == ExistingFamily.PRIMARY) {
            validateExistingKeys(existing, primaryKeys, SENSITIVITY_F_OUTPUT);
            if (containsNonOptimal(existing, primaryKeys)) {
                resetOutputFile(SENSITIVITY_F_OUTPUT, SENSITIVITY_F_HEADER);
                evaluateSensitivityF(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, emptyExistingOutput());
                return;
            }
        }

        FamilyRun primary = evaluateSensitivityF(toInstancePaths(SENSITIVITY_BASE_DEFAULT), options.timeLimitSec, existing);
        if (!primary.allOptimal) {
            resetOutputFile(SENSITIVITY_F_OUTPUT, SENSITIVITY_F_HEADER);
            evaluateSensitivityF(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, emptyExistingOutput());
        }
    }

    private static void runSensitivityH(CliOptions options) throws Exception {
        ExistingOutput existing = loadExistingOutput(SENSITIVITY_H_OUTPUT, SENSITIVITY_H_HEADER, new int[]{0, 1}, 3);
        if (options.hasExplicitInstances()) {
            evaluateSensitivityH(primarySensitivityInstances(options), options.timeLimitSec, existing);
            return;
        }

        ExistingFamily family = detectSensitivityFamily(existing, SENSITIVITY_BASE_DEFAULT, SENSITIVITY_BASE_FALLBACK);
        if (family == ExistingFamily.FALLBACK) {
            validateExistingKeys(existing, buildSensitivityHTargetKeys(toInstancePaths(SENSITIVITY_BASE_FALLBACK)), SENSITIVITY_H_OUTPUT);
            evaluateSensitivityH(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, existing);
            return;
        }

        ArrayList<String> primaryKeys = buildSensitivityHTargetKeys(toInstancePaths(SENSITIVITY_BASE_DEFAULT));
        if (family == ExistingFamily.PRIMARY) {
            validateExistingKeys(existing, primaryKeys, SENSITIVITY_H_OUTPUT);
            if (containsNonOptimal(existing, primaryKeys)) {
                resetOutputFile(SENSITIVITY_H_OUTPUT, SENSITIVITY_H_HEADER);
                evaluateSensitivityH(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, emptyExistingOutput());
                return;
            }
        }

        FamilyRun primary = evaluateSensitivityH(toInstancePaths(SENSITIVITY_BASE_DEFAULT), options.timeLimitSec, existing);
        if (!primary.allOptimal) {
            resetOutputFile(SENSITIVITY_H_OUTPUT, SENSITIVITY_H_HEADER);
            evaluateSensitivityH(toInstancePaths(SENSITIVITY_BASE_FALLBACK), options.timeLimitSec, emptyExistingOutput());
        }
    }

    private static FamilyRun evaluateSensitivityQ(String[] instancePaths, double timeLimitSec, ExistingOutput existing) throws Exception {
        ArrayList<String[]> rows = new ArrayList<String[]>();
        boolean allOptimal = true;
        LbbdReformulationSolver.Config config = fullRunConfig(timeLimitSec);
        int total = instancePaths.length * Q_RATIOS.length;
        int progress = countCompletedTargets(existing, buildSensitivityQTargetKeys(instancePaths));

        ensureOutputFile(SENSITIVITY_Q_OUTPUT, SENSITIVITY_Q_HEADER);
        try (BufferedWriter writer = appendWriter(SENSITIVITY_Q_OUTPUT)) {
            for (int idx = 0; idx < instancePaths.length; idx++) {
                String instancePath = normalizeInstancePath(instancePaths[idx]);
                Instance base = null;
                String name = instanceName(instancePath);
                for (int qIdx = 0; qIdx < Q_RATIOS.length; qIdx++) {
                    double ratio = Q_RATIOS[qIdx];
                    String key = buildKey(name, fmtRatio(ratio));
                    ExistingRow existingRow = existing.rowsByKey.get(key);
                    if (existingRow != null) {
                        rows.add(existingRow.values);
                        allOptimal &= existingRow.optimal;
                        System.out.println("[" + progress + "/" + total + "] skip: sensitivity_Q " + name + " ratio=" + fmtRatio(ratio));
                        continue;
                    }
                    if (base == null) {
                        base = loadInstance(instancePath);
                    }
                    int qValue = (int) Math.floor(base.Q * ratio);
                    Instance.Overrides overrides = Instance.Overrides.defaults();
                    overrides.Q = (double) qValue;
                    Instance scenario = base.copyWithOverrides(overrides);
                    SolveWithDetails run = solveQuietly(scenario, config);
                    Metrics metrics = Metrics.from(scenario, run.details);
                    String[] row = new String[]{
                            name,
                            fmtRatio(ratio),
                            Integer.toString(qValue),
                            fmtBoolean(run.summary.feasible),
                            fmtBoolean(run.summary.optimal),
                            fmt(run.summary.objective),
                            fmt(run.summary.solveTimeSec),
                            fmt(metrics.routeCost),
                            fmt(metrics.inventoryCost),
                            fmt(metrics.productionCost),
                            fmt(metrics.avgVisitFrequency),
                            fmt(metrics.avgPickupPerVisit),
                    };
                    writeCsvRow(writer, row);
                    existing.rowsByKey.put(key, new ExistingRow(row, 4));
                    rows.add(row);
                    allOptimal &= run.summary.optimal;
                    progress++;
                    System.out.println("[" + progress + "/" + total + "] sensitivity_Q " + name + " Q=" + qValue);
                }
            }
        }
        return new FamilyRun(rows, allOptimal);
    }

    private static FamilyRun evaluateSensitivityK(String[] instancePaths, double timeLimitSec, ExistingOutput existing) throws Exception {
        ArrayList<String[]> rows = new ArrayList<String[]>();
        boolean allOptimal = true;
        LbbdReformulationSolver.Config config = fullRunConfig(timeLimitSec);
        int total = instancePaths.length * K_VALUES.length;
        int progress = countCompletedTargets(existing, buildSensitivityKTargetKeys(instancePaths));

        ensureOutputFile(SENSITIVITY_K_OUTPUT, SENSITIVITY_K_HEADER);
        try (BufferedWriter writer = appendWriter(SENSITIVITY_K_OUTPUT)) {
            for (int idx = 0; idx < instancePaths.length; idx++) {
                String instancePath = normalizeInstancePath(instancePaths[idx]);
                Instance base = null;
                String name = instanceName(instancePath);
                for (int kIdx = 0; kIdx < K_VALUES.length; kIdx++) {
                    int kValue = K_VALUES[kIdx];
                    String key = buildKey(name, Integer.toString(kValue));
                    ExistingRow existingRow = existing.rowsByKey.get(key);
                    if (existingRow != null) {
                        rows.add(existingRow.values);
                        allOptimal &= existingRow.optimal;
                        System.out.println("[" + progress + "/" + total + "] skip: sensitivity_K " + name + " K=" + kValue);
                        continue;
                    }
                    if (base == null) {
                        base = loadInstance(instancePath);
                    }
                    Instance.Overrides overrides = Instance.Overrides.defaults();
                    overrides.K = Integer.valueOf(kValue);
                    Instance scenario = base.copyWithOverrides(overrides);
                    SolveWithDetails run = solveQuietly(scenario, config);
                    Metrics metrics = Metrics.from(scenario, run.details);
                    String[] row = new String[]{
                            name,
                            Integer.toString(kValue),
                            fmtBoolean(run.summary.feasible),
                            fmtBoolean(run.summary.optimal),
                            fmt(run.summary.objective),
                            fmt(run.summary.solveTimeSec),
                            fmt(metrics.routeCost),
                            fmt(metrics.inventoryCost),
                            fmt(metrics.productionCost),
                            fmt(metrics.avgVisitFrequency),
                            fmt(metrics.avgPickupPerVisit),
                    };
                    writeCsvRow(writer, row);
                    existing.rowsByKey.put(key, new ExistingRow(row, 3));
                    rows.add(row);
                    allOptimal &= run.summary.optimal;
                    progress++;
                    System.out.println("[" + progress + "/" + total + "] sensitivity_K " + name + " K=" + kValue);
                }
            }
        }
        return new FamilyRun(rows, allOptimal);
    }

    private static FamilyRun evaluateSensitivityF(String[] instancePaths, double timeLimitSec, ExistingOutput existing) throws Exception {
        ArrayList<String[]> rows = new ArrayList<String[]>();
        boolean allOptimal = true;
        LbbdReformulationSolver.Config config = fullRunConfig(timeLimitSec);
        int total = instancePaths.length * F_RATIOS.length;
        int progress = countCompletedTargets(existing, buildSensitivityFTargetKeys(instancePaths));

        ensureOutputFile(SENSITIVITY_F_OUTPUT, SENSITIVITY_F_HEADER);
        try (BufferedWriter writer = appendWriter(SENSITIVITY_F_OUTPUT)) {
            for (int idx = 0; idx < instancePaths.length; idx++) {
                String instancePath = normalizeInstancePath(instancePaths[idx]);
                Instance base = null;
                String name = instanceName(instancePath);
                for (int fIdx = 0; fIdx < F_RATIOS.length; fIdx++) {
                    double ratio = F_RATIOS[fIdx];
                    String key = buildKey(name, fmtRatio(ratio));
                    ExistingRow existingRow = existing.rowsByKey.get(key);
                    if (existingRow != null) {
                        rows.add(existingRow.values);
                        allOptimal &= existingRow.optimal;
                        System.out.println("[" + progress + "/" + total + "] skip: sensitivity_f " + name + " ratio=" + fmtRatio(ratio));
                        continue;
                    }
                    if (base == null) {
                        base = loadInstance(instancePath);
                    }
                    double fValue = base.f * ratio;
                    Instance.Overrides overrides = Instance.Overrides.defaults();
                    overrides.f = Double.valueOf(fValue);
                    Instance scenario = base.copyWithOverrides(overrides);
                    SolveWithDetails run = solveQuietly(scenario, config);
                    Metrics metrics = Metrics.from(scenario, run.details);
                    String[] row = new String[]{
                            name,
                            fmtRatio(ratio),
                            fmt(fValue),
                            fmtBoolean(run.summary.feasible),
                            fmtBoolean(run.summary.optimal),
                            fmt(run.summary.objective),
                            fmt(run.summary.solveTimeSec),
                            fmt(metrics.routeCost),
                            fmt(metrics.inventoryCost),
                            fmt(metrics.productionCost),
                            Integer.toString(metrics.numProductionPeriods),
                            fmt(metrics.avgProductionBatch),
                            fmt(metrics.rawInventoryAvg),
                            fmt(metrics.finishedInventoryAvg),
                    };
                    writeCsvRow(writer, row);
                    existing.rowsByKey.put(key, new ExistingRow(row, 4));
                    rows.add(row);
                    allOptimal &= run.summary.optimal;
                    progress++;
                    System.out.println("[" + progress + "/" + total + "] sensitivity_f " + name + " f=" + fmt(fValue));
                }
            }
        }
        return new FamilyRun(rows, allOptimal);
    }

    private static FamilyRun evaluateSensitivityH(String[] instancePaths, double timeLimitSec, ExistingOutput existing) throws Exception {
        ArrayList<String[]> rows = new ArrayList<String[]>();
        boolean allOptimal = true;
        LbbdReformulationSolver.Config config = fullRunConfig(timeLimitSec);
        int total = instancePaths.length * H_RATIOS.length;
        int progress = countCompletedTargets(existing, buildSensitivityHTargetKeys(instancePaths));

        ensureOutputFile(SENSITIVITY_H_OUTPUT, SENSITIVITY_H_HEADER);
        try (BufferedWriter writer = appendWriter(SENSITIVITY_H_OUTPUT)) {
            for (int idx = 0; idx < instancePaths.length; idx++) {
                String instancePath = normalizeInstancePath(instancePaths[idx]);
                Instance base = null;
                String name = instanceName(instancePath);
                for (int hIdx = 0; hIdx < H_RATIOS.length; hIdx++) {
                    double ratio = H_RATIOS[hIdx];
                    String key = buildKey(name, fmtRatio(ratio));
                    ExistingRow existingRow = existing.rowsByKey.get(key);
                    if (existingRow != null) {
                        rows.add(existingRow.values);
                        allOptimal &= existingRow.optimal;
                        System.out.println("[" + progress + "/" + total + "] skip: sensitivity_h " + name + " ratio=" + fmtRatio(ratio));
                        continue;
                    }
                    if (base == null) {
                        base = loadInstance(instancePath);
                    }
                    Instance.Overrides overrides = Instance.Overrides.defaults();
                    overrides.h0 = Double.valueOf(base.h0 * ratio);
                    overrides.hp = Double.valueOf(base.hp * ratio);
                    overrides.hi = scaleSupplierHoldCosts(base, ratio);
                    Instance scenario = base.copyWithOverrides(overrides);
                    SolveWithDetails run = solveQuietly(scenario, config);
                    Metrics metrics = Metrics.from(scenario, run.details);
                    String[] row = new String[]{
                            name,
                            fmtRatio(ratio),
                            fmtBoolean(run.summary.feasible),
                            fmtBoolean(run.summary.optimal),
                            fmt(run.summary.objective),
                            fmt(run.summary.solveTimeSec),
                            fmt(metrics.routeCost),
                            fmt(metrics.inventoryCost),
                            fmt(metrics.productionCost),
                            fmt(metrics.supplierInventoryAvg),
                            fmt(metrics.rawInventoryAvg),
                            fmt(metrics.finishedInventoryAvg),
                            fmt(metrics.avgVisitFrequency),
                    };
                    writeCsvRow(writer, row);
                    existing.rowsByKey.put(key, new ExistingRow(row, 3));
                    rows.add(row);
                    allOptimal &= run.summary.optimal;
                    progress++;
                    System.out.println("[" + progress + "/" + total + "] sensitivity_h " + name + " h-ratio=" + fmtRatio(ratio));
                }
            }
        }
        return new FamilyRun(rows, allOptimal);
    }

    private static SolveWithDetails solveQuietly(final Instance ins, final LbbdReformulationSolver.Config config) throws Exception {
        return runQuietly(new CheckedSupplier<SolveWithDetails>() {
            @Override
            public SolveWithDetails get() {
                LbbdReformulationSolver solver = new LbbdReformulationSolver(ins, config);
                SolveResult summary = solver.solve(Double.NaN, CplexConfig.MIP_GAP);
                return new SolveWithDetails(summary, solver.getLastDetailedResult());
            }
        });
    }

    private static HashMap<String, FullResultRow> loadFullResults() throws IOException {
        HashMap<String, FullResultRow> rows = new HashMap<String, FullResultRow>();
        for (int idx = 0; idx < FULL_RESULT_SOURCES.length; idx++) {
            Path source = FULL_RESULT_SOURCES[idx];
            try (BufferedReader reader = Files.newBufferedReader(source, StandardCharsets.UTF_8)) {
                String headerLine = reader.readLine();
                if (headerLine == null) {
                    continue;
                }
                List<String> header = parseCsvLine(headerLine);
                HashMap<String, Integer> indexByName = new HashMap<String, Integer>();
                for (int i = 0; i < header.size(); i++) {
                    indexByName.put(header.get(i), Integer.valueOf(i));
                }
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.trim().isEmpty()) {
                        continue;
                    }
                    List<String> fields = parseCsvLine(line);
                    String method = field(fields, indexByName, "method");
                    if (!"lbbd_no_init".equals(method)) {
                        continue;
                    }
                    String instance = field(fields, indexByName, "instance");
                    String status = field(fields, indexByName, "status");
                    rows.put(instance, new FullResultRow(
                            parseBoolean(field(fields, indexByName, "feasible")),
                            parseBoolean(field(fields, indexByName, "optimal")),
                            parseDouble(field(fields, indexByName, "objective")),
                            parseDouble(field(fields, indexByName, "best_bound")),
                            parseDouble(field(fields, indexByName, "gap")),
                            parseDouble(field(fields, indexByName, "time_sec")),
                            parseStatusMetric(status, ITER_PATTERN),
                            parseStatusMetric(status, FEAS_CUT_PATTERN),
                            parseStatusMetric(status, OPT_CUT_PATTERN)
                    ));
                }
            }
        }
        return rows;
    }

    private static String field(List<String> fields, HashMap<String, Integer> indexByName, String name) {
        Integer index = indexByName.get(name);
        if (index == null || index.intValue() >= fields.size()) {
            throw new IllegalStateException("Missing CSV column: " + name);
        }
        return fields.get(index.intValue());
    }

    private static int parseStatusMetric(String status, Pattern pattern) {
        if (status == null) {
            return 0;
        }
        Matcher matcher = pattern.matcher(status);
        if (!matcher.find()) {
            return 0;
        }
        return Integer.parseInt(matcher.group(1));
    }

    private static String[] primarySensitivityInstances(CliOptions options) {
        return options.hasExplicitInstances()
                ? options.explicitInstances
                : toInstancePaths(SENSITIVITY_BASE_DEFAULT);
    }

    private static double[] scaleSupplierHoldCosts(Instance ins, double ratio) {
        double[] scaled = new double[ins.n + 1];
        scaled[0] = Double.NaN;
        for (int i = 1; i <= ins.n; i++) {
            scaled[i] = ins.hi[i] * ratio;
        }
        return scaled;
    }

    private static LbbdReformulationSolver.Config fullRunConfig(double timeLimitSec) {
        return new LbbdReformulationSolver.Config(false, true, true, true, timeLimitSec);
    }

    private static ExistingOutput emptyExistingOutput() {
        return new ExistingOutput();
    }

    private static ExistingOutput loadExistingOutput(Path path, String expectedHeader, int[] keyColumns, int optimalColumn) throws IOException {
        ExistingOutput output = new ExistingOutput();
        if (!Files.exists(path)) {
            return output;
        }
        try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String header = reader.readLine();
            if (header == null || header.trim().isEmpty()) {
                return output;
            }
            if (!expectedHeader.equals(header)) {
                throw new IllegalStateException("Unexpected CSV header in " + path + ": " + header);
            }
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    continue;
                }
                List<String> fields = parseCsvLine(line);
                String[] values = fields.toArray(new String[0]);
                output.rowsByKey.put(buildKey(values, keyColumns), new ExistingRow(values, optimalColumn));
            }
        }
        return output;
    }

    private static void ensureOutputFile(Path path, String header) throws IOException {
        if (!Files.exists(path) || Files.size(path) == 0L) {
            resetOutputFile(path, header);
            return;
        }
        try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String actualHeader = reader.readLine();
            if (actualHeader == null || actualHeader.trim().isEmpty()) {
                resetOutputFile(path, header);
                return;
            }
            if (!header.equals(actualHeader)) {
                throw new IllegalStateException("Unexpected CSV header in " + path + ": " + actualHeader);
            }
        }
    }

    private static void resetOutputFile(Path path, String header) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(
                path,
                StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING,
                StandardOpenOption.WRITE
        )) {
            writer.write(header);
            writer.write(System.lineSeparator());
            writer.flush();
        }
    }

    private static BufferedWriter appendWriter(Path path) throws IOException {
        return Files.newBufferedWriter(
                path,
                StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,
                StandardOpenOption.APPEND,
                StandardOpenOption.WRITE
        );
    }

    private static Instance loadInstance(String instancePath) throws IOException {
        Instance.Options options = Instance.Options.defaults();
        options.distanceMode = Instance.Options.DistanceMode.EUCLIDEAN_FLOAT;
        options.autoSetDt = true;
        return Instance.fromFile(instancePath, options);
    }

    private static String[] toInstancePaths(String[] instanceNames) {
        String[] paths = new String[instanceNames.length];
        for (int i = 0; i < instanceNames.length; i++) {
            paths[i] = normalizeInstancePath(instanceNames[i]);
        }
        return paths;
    }

    private static String normalizeInstancePath(String raw) {
        if (raw == null) {
            throw new IllegalArgumentException("instance path must not be null");
        }
        String trimmed = raw.trim();
        if (trimmed.indexOf('/') >= 0 || trimmed.endsWith(".txt")) {
            return trimmed;
        }
        return "data/MVPRP/" + trimmed + ".txt";
    }

    private static String instanceName(String instancePath) {
        String fileName = Paths.get(instancePath).getFileName().toString();
        if (fileName.endsWith(".txt")) {
            return fileName.substring(0, fileName.length() - 4);
        }
        return fileName;
    }

    private static ArrayList<String> buildAblationTargetKeys(String[] instancePaths, VariantSpec[] variants) {
        ArrayList<String> keys = new ArrayList<String>();
        for (int idx = 0; idx < instancePaths.length; idx++) {
            String name = instanceName(normalizeInstancePath(instancePaths[idx]));
            for (int v = 0; v < variants.length; v++) {
                keys.add(buildKey(name, variants[v].name));
            }
            keys.add(buildKey(name, "full"));
        }
        return keys;
    }

    private static ArrayList<String> buildSensitivityQTargetKeys(String[] instancePaths) {
        ArrayList<String> keys = new ArrayList<String>();
        for (int idx = 0; idx < instancePaths.length; idx++) {
            String name = instanceName(normalizeInstancePath(instancePaths[idx]));
            for (int qIdx = 0; qIdx < Q_RATIOS.length; qIdx++) {
                keys.add(buildKey(name, fmtRatio(Q_RATIOS[qIdx])));
            }
        }
        return keys;
    }

    private static ArrayList<String> buildSensitivityKTargetKeys(String[] instancePaths) {
        ArrayList<String> keys = new ArrayList<String>();
        for (int idx = 0; idx < instancePaths.length; idx++) {
            String name = instanceName(normalizeInstancePath(instancePaths[idx]));
            for (int kIdx = 0; kIdx < K_VALUES.length; kIdx++) {
                keys.add(buildKey(name, Integer.toString(K_VALUES[kIdx])));
            }
        }
        return keys;
    }

    private static ArrayList<String> buildSensitivityFTargetKeys(String[] instancePaths) {
        ArrayList<String> keys = new ArrayList<String>();
        for (int idx = 0; idx < instancePaths.length; idx++) {
            String name = instanceName(normalizeInstancePath(instancePaths[idx]));
            for (int fIdx = 0; fIdx < F_RATIOS.length; fIdx++) {
                keys.add(buildKey(name, fmtRatio(F_RATIOS[fIdx])));
            }
        }
        return keys;
    }

    private static ArrayList<String> buildSensitivityHTargetKeys(String[] instancePaths) {
        ArrayList<String> keys = new ArrayList<String>();
        for (int idx = 0; idx < instancePaths.length; idx++) {
            String name = instanceName(normalizeInstancePath(instancePaths[idx]));
            for (int hIdx = 0; hIdx < H_RATIOS.length; hIdx++) {
                keys.add(buildKey(name, fmtRatio(H_RATIOS[hIdx])));
            }
        }
        return keys;
    }

    private static int countCompletedTargets(ExistingOutput existing, List<String> targetKeys) {
        int count = 0;
        for (int i = 0; i < targetKeys.size(); i++) {
            if (existing.rowsByKey.containsKey(targetKeys.get(i))) {
                count++;
            }
        }
        return count;
    }

    private static boolean containsNonOptimal(ExistingOutput existing, List<String> targetKeys) {
        for (int i = 0; i < targetKeys.size(); i++) {
            ExistingRow row = existing.rowsByKey.get(targetKeys.get(i));
            if (row != null && !row.optimal) {
                return true;
            }
        }
        return false;
    }

    private static void validateExistingKeys(ExistingOutput existing, List<String> targetKeys, Path outputPath) {
        HashSet<String> allowed = new HashSet<String>(targetKeys);
        for (String key : existing.rowsByKey.keySet()) {
            if (!allowed.contains(key)) {
                throw new IllegalStateException("Unexpected existing row in " + outputPath + ": " + key);
            }
        }
    }

    private static ExistingFamily detectSensitivityFamily(
            ExistingOutput existing,
            String[] primaryInstances,
            String[] fallbackInstances
    ) {
        Set<String> primary = toNameSet(primaryInstances);
        Set<String> fallback = toNameSet(fallbackInstances);
        boolean sawPrimary = false;
        boolean sawFallback = false;

        for (ExistingRow row : existing.rowsByKey.values()) {
            String instance = row.values.length == 0 ? "" : row.values[0];
            if (primary.contains(instance)) {
                sawPrimary = true;
            } else if (fallback.contains(instance)) {
                sawFallback = true;
            } else {
                throw new IllegalStateException("Unexpected instance in sensitivity output: " + instance);
            }
        }
        if (sawPrimary && sawFallback) {
            throw new IllegalStateException("Mixed primary/fallback rows in sensitivity output");
        }
        if (sawFallback) {
            return ExistingFamily.FALLBACK;
        }
        if (sawPrimary) {
            return ExistingFamily.PRIMARY;
        }
        return ExistingFamily.NONE;
    }

    private static Set<String> toNameSet(String[] instanceNames) {
        HashSet<String> names = new HashSet<String>();
        for (int i = 0; i < instanceNames.length; i++) {
            names.add(instanceName(normalizeInstancePath(instanceNames[i])));
        }
        return names;
    }

    private static String buildKey(String... parts) {
        StringBuilder key = new StringBuilder();
        for (int i = 0; i < parts.length; i++) {
            if (i > 0) {
                key.append('\u0001');
            }
            key.append(parts[i] == null ? "" : parts[i]);
        }
        return key.toString();
    }

    private static String buildKey(String[] values, int[] indices) {
        String[] parts = new String[indices.length];
        for (int i = 0; i < indices.length; i++) {
            int index = indices[i];
            if (index < 0 || index >= values.length) {
                throw new IllegalStateException("Missing key column at index " + index);
            }
            parts[i] = values[index];
        }
        return buildKey(parts);
    }

    private static void writeCsvRow(BufferedWriter writer, String[] values) throws IOException {
        StringBuilder sb = new StringBuilder(256);
        for (int i = 0; i < values.length; i++) {
            if (i > 0) {
                sb.append(',');
            }
            appendCsvField(sb, values[i]);
        }
        sb.append(System.lineSeparator());
        writer.write(sb.toString());
        writer.flush();
    }

    private static void appendCsvField(StringBuilder sb, String value) {
        String safe = (value == null) ? "" : value;
        sb.append('"');
        for (int i = 0; i < safe.length(); i++) {
            char ch = safe.charAt(i);
            if (ch == '"') {
                sb.append('"');
            }
            sb.append(ch);
        }
        sb.append('"');
    }

    private static List<String> parseCsvLine(String line) {
        ArrayList<String> fields = new ArrayList<String>();
        StringBuilder current = new StringBuilder();
        boolean inQuotes = false;
        for (int i = 0; i < line.length(); i++) {
            char ch = line.charAt(i);
            if (ch == '"') {
                if (inQuotes && i + 1 < line.length() && line.charAt(i + 1) == '"') {
                    current.append('"');
                    i++;
                } else {
                    inQuotes = !inQuotes;
                }
            } else if (ch == ',' && !inQuotes) {
                fields.add(current.toString());
                current.setLength(0);
            } else {
                current.append(ch);
            }
        }
        fields.add(current.toString());
        return fields;
    }

    private static boolean parseBoolean(String text) {
        return Boolean.parseBoolean(text);
    }

    private static double parseDouble(String text) {
        if (text == null || text.trim().isEmpty() || "NaN".equalsIgnoreCase(text.trim())) {
            return Double.NaN;
        }
        return Double.parseDouble(text.trim());
    }

    private static String fmt(double value) {
        if (Double.isNaN(value)) {
            return "NaN";
        }
        if (Double.isInfinite(value)) {
            return value > 0.0 ? "Infinity" : "-Infinity";
        }
        return String.format(Locale.US, "%.6f", value);
    }

    private static String fmtRatio(double value) {
        return String.format(Locale.US, "%.2f", value);
    }

    private static String fmtBoolean(boolean value) {
        return value ? "True" : "False";
    }

    private static <T> T runQuietly(CheckedSupplier<T> supplier) throws Exception {
        PrintStream originalOut = System.out;
        PrintStream originalErr = System.err;
        PrintStream silent = new PrintStream(DEV_NULL);
        try {
            System.setOut(silent);
            System.setErr(silent);
            return supplier.get();
        } finally {
            System.setOut(originalOut);
            System.setErr(originalErr);
            silent.close();
        }
    }

    private interface CheckedSupplier<T> {
        T get() throws Exception;
    }

    private static final class VariantSpec {
        final String name;
        final LbbdReformulationSolver.Config config;

        VariantSpec(String name, LbbdReformulationSolver.Config config) {
            this.name = name;
            this.config = config;
        }
    }

    private static final class SolveWithDetails {
        final SolveResult summary;
        final LbbdDetailedResult details;

        SolveWithDetails(SolveResult summary, LbbdDetailedResult details) {
            this.summary = summary;
            this.details = details;
        }

        int iterations() {
            return details == null ? 0 : details.iterations;
        }

        int feasibilityCuts() {
            return details == null ? 0 : details.feasibilityCuts;
        }

        int optimalityCuts() {
            return details == null ? 0 : details.optimalityCuts;
        }
    }

    private static final class FullResultRow {
        final boolean feasible;
        final boolean optimal;
        final double objective;
        final double bestBound;
        final double gap;
        final double timeSec;
        final int iterations;
        final int feasibilityCuts;
        final int optimalityCuts;

        FullResultRow(
                boolean feasible,
                boolean optimal,
                double objective,
                double bestBound,
                double gap,
                double timeSec,
                int iterations,
                int feasibilityCuts,
                int optimalityCuts
        ) {
            this.feasible = feasible;
            this.optimal = optimal;
            this.objective = objective;
            this.bestBound = bestBound;
            this.gap = gap;
            this.timeSec = timeSec;
            this.iterations = iterations;
            this.feasibilityCuts = feasibilityCuts;
            this.optimalityCuts = optimalityCuts;
        }
    }

    private static final class FamilyRun {
        final ArrayList<String[]> rows;
        final boolean allOptimal;

        FamilyRun(ArrayList<String[]> rows, boolean allOptimal) {
            this.rows = rows;
            this.allOptimal = allOptimal;
        }
    }

    private static final class ExistingOutput {
        final HashMap<String, ExistingRow> rowsByKey = new HashMap<String, ExistingRow>();
    }

    private static final class ExistingRow {
        final String[] values;
        final boolean optimal;

        ExistingRow(String[] values, int optimalColumn) {
            this.values = values;
            this.optimal = optimalColumn >= 0
                    && optimalColumn < values.length
                    && parseBoolean(values[optimalColumn]);
        }
    }

    private enum ExistingFamily {
        NONE,
        PRIMARY,
        FALLBACK
    }

    private static final class Metrics {
        final double routeCost;
        final double inventoryCost;
        final double productionCost;
        final double avgVisitFrequency;
        final double avgPickupPerVisit;
        final int numProductionPeriods;
        final double avgProductionBatch;
        final double rawInventoryAvg;
        final double finishedInventoryAvg;
        final double supplierInventoryAvg;

        Metrics(
                double routeCost,
                double inventoryCost,
                double productionCost,
                double avgVisitFrequency,
                double avgPickupPerVisit,
                int numProductionPeriods,
                double avgProductionBatch,
                double rawInventoryAvg,
                double finishedInventoryAvg,
                double supplierInventoryAvg
        ) {
            this.routeCost = routeCost;
            this.inventoryCost = inventoryCost;
            this.productionCost = productionCost;
            this.avgVisitFrequency = avgVisitFrequency;
            this.avgPickupPerVisit = avgPickupPerVisit;
            this.numProductionPeriods = numProductionPeriods;
            this.avgProductionBatch = avgProductionBatch;
            this.rawInventoryAvg = rawInventoryAvg;
            this.finishedInventoryAvg = finishedInventoryAvg;
            this.supplierInventoryAvg = supplierInventoryAvg;
        }

        static Metrics from(Instance ins, LbbdDetailedResult details) {
            if (details == null) {
                return new Metrics(Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
            }
            double[][] supplierEndInventory = LbbdSolutionAnalytics.computeSupplierEndInventory(ins, details.qBar, details.zBar);
            return new Metrics(
                    details.routeCost,
                    details.inventoryCost,
                    details.productionCost,
                    LbbdSolutionAnalytics.computeAverageVisitFrequency(ins, details.zBar),
                    LbbdSolutionAnalytics.computeAveragePickupPerVisit(ins, details.qBar, details.zBar),
                    LbbdSolutionAnalytics.countProductionPeriods(details.y),
                    LbbdSolutionAnalytics.computeAverageProductionBatch(details.p, details.y),
                    LbbdSolutionAnalytics.computeRawInventoryAverage(ins, details.i0),
                    LbbdSolutionAnalytics.computeFinishedInventoryAverage(ins, details.p0),
                    LbbdSolutionAnalytics.computeSupplierInventoryAverage(ins, supplierEndInventory)
            );
        }
    }

    private static final class CliOptions {
        final String mode;
        final double timeLimitSec;
        final String[] explicitInstances;

        CliOptions(String mode, double timeLimitSec, String[] explicitInstances) {
            this.mode = mode;
            this.timeLimitSec = timeLimitSec;
            this.explicitInstances = explicitInstances;
        }

        boolean hasExplicitInstances() {
            return explicitInstances != null && explicitInstances.length > 0;
        }

        static CliOptions parse(String[] args) {
            String mode = "all";
            double timeLimitSec = DEFAULT_TIME_LIMIT_SEC;
            String[] explicitInstances = null;
            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (arg.startsWith("--time-limit=")) {
                    timeLimitSec = Double.parseDouble(arg.substring("--time-limit=".length()));
                } else if (arg.startsWith("--instances=")) {
                    explicitInstances = splitCsvValues(arg.substring("--instances=".length()));
                } else if (!arg.startsWith("--")) {
                    mode = normalizeMode(arg);
                } else {
                    throw new IllegalArgumentException("Unsupported option: " + arg);
                }
            }
            return new CliOptions(mode, timeLimitSec, explicitInstances);
        }

        private static String normalizeMode(String raw) {
            String normalized = raw.trim().toLowerCase(Locale.ROOT);
            if ("sensitivity_q".equals(normalized) || "sensitivity-q".equals(normalized)) {
                return "sensitivity_q";
            }
            if ("sensitivity_k".equals(normalized) || "sensitivity-k".equals(normalized)) {
                return "sensitivity_k";
            }
            if ("sensitivity_f".equals(normalized) || "sensitivity-f".equals(normalized)) {
                return "sensitivity_f";
            }
            if ("sensitivity_h".equals(normalized) || "sensitivity-h".equals(normalized)) {
                return "sensitivity_h";
            }
            return normalized;
        }

        private static String[] splitCsvValues(String raw) {
            if (raw == null || raw.trim().isEmpty()) {
                return new String[0];
            }
            String[] tokens = raw.split(",");
            ArrayList<String> values = new ArrayList<String>();
            for (int i = 0; i < tokens.length; i++) {
                String token = tokens[i].trim();
                if (!token.isEmpty()) {
                    values.add(token);
                }
            }
            return values.toArray(new String[0]);
        }
    }
}
