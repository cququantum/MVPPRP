import instance.Instance;
import lbbdModel.LbbdReformulationSolver;
import model.CplexConfig;
import model.SolveResult;
import originalModel.OriginalModelSolver;
import reformulationModel.ReformulationModelSolver;

import java.io.IOException;
import java.util.Locale;

public class Main {
    private enum SolverMode { ALL, ORIGINAL, REFORM, LBBD_NO_INIT, LBBD }
    private static final double RESULT_TOL = CplexConfig.MIP_GAP;

    public static void main(String[] args) {
        SolverMode solverMode = SolverMode.LBBD_NO_INIT;
        int startArg = 0;
        if (args.length > 0 && args[0].startsWith("--solver=")) {
            solverMode = parseSolverMode(args[0].substring("--solver=".length()));
            startArg = 1;
        }


        String[] instancePaths;
        if (args.length > startArg) {
            instancePaths = new String[args.length - startArg];
            System.arraycopy(args, startArg, instancePaths, 0, instancePaths.length);
        } else {
            instancePaths = new String[]{"data/MVPRP/MVPRP1_15_6_2.txt"};
        }

        for (int idx = 0; idx < instancePaths.length; idx++) {
            if (idx > 0) {
                System.out.println();
            }
            runSingleInstance(instancePaths[idx], solverMode);
        }
    }

    private static void runSingleInstance(String instancePath, SolverMode solverMode) {
        Instance.Options options = Instance.Options.defaults();
        options.distanceMode = Instance.Options.DistanceMode.EUCLIDEAN_FLOAT;
        options.autoSetDt = true;

        try {
            Instance ins = Instance.fromFile(instancePath, options);

            SolveResult originalResult = null;
            SolveResult reformulationResult = null;
            SolveResult lbbdNoInitResult = null;
            SolveResult lbbdResult = null;

            if (solverMode == SolverMode.ALL || solverMode == SolverMode.ORIGINAL) {
                originalResult = new OriginalModelSolver().solve(ins);
            }
            if (solverMode == SolverMode.ALL || solverMode == SolverMode.REFORM) {
                ReformulationModelSolver reformulationSolver = new ReformulationModelSolver();
                reformulationResult = reformulationSolver.solve(ins);
            }
            if (solverMode == SolverMode.ALL || solverMode == SolverMode.LBBD_NO_INIT) {
                LbbdReformulationSolver lbbdNoInitSolver = new LbbdReformulationSolver(ins, false);
                lbbdNoInitResult = lbbdNoInitSolver.solve(Double.NaN, RESULT_TOL);
            }
            if (solverMode == SolverMode.ALL || solverMode == SolverMode.LBBD) {
                LbbdReformulationSolver lbbdSolver = new LbbdReformulationSolver(ins);
                lbbdResult = lbbdSolver.solve(Double.NaN, RESULT_TOL);
            }

            System.out.println("Instance: " + instancePath);
            System.out.println("n=" + ins.n + ", l=" + ins.l + ", K=" + ins.K + ", Q=" + format(ins.Q));
            if (originalResult != null) {
                System.out.println(originalResult.toSummaryLine());
            }
            if (reformulationResult != null) {
                System.out.println(reformulationResult.toSummaryLine());
            }
            if (lbbdNoInitResult != null) {
                System.out.println(lbbdNoInitResult.toSummaryLine());
            }
            if (lbbdResult != null) {
                System.out.println(lbbdResult.toSummaryLine());
            }
            if (originalResult != null && reformulationResult != null) {
                printComparison(originalResult, reformulationResult);
            }
            if (reformulationResult != null && lbbdNoInitResult != null) {
                printValidation("reform-vs-lbbd_no_init", reformulationResult, lbbdNoInitResult, CplexConfig.MIP_GAP);
            }
            if (reformulationResult != null && lbbdResult != null) {
                printValidation("reform-vs-lbbd", reformulationResult, lbbdResult, CplexConfig.MIP_GAP);
            }
        } catch (IOException e) {
            throw new RuntimeException("Failed to load instance file: " + instancePath, e);
        }
    }

    private static SolverMode parseSolverMode(String s) {
        if (s == null) {
            return SolverMode.ALL;
        }
        String x = s.trim().toLowerCase(Locale.ROOT);
        if ("all".equals(x)) {
            return SolverMode.ALL;
        }
        if ("original".equals(x)) {
            return SolverMode.ORIGINAL;
        }
        if ("reform".equals(x) || "reformulation".equals(x)) {
            return SolverMode.REFORM;
        }
        if ("lbbd_no_init".equals(x) || "lbbd-no-init".equals(x) || "lbbd_noinit".equals(x)) {
            return SolverMode.LBBD_NO_INIT;
        }
        if ("lbbd".equals(x)) {
            return SolverMode.LBBD;
        }
        throw new IllegalArgumentException("Unsupported --solver mode: " + s);
    }

    private static void printComparison(SolveResult original, SolveResult reformulation) {
        if (!hasObjective(original) || !hasObjective(reformulation)) {
            System.out.println("Comparison: at least one model has no feasible incumbent objective.");
            return;
        }

        double objDiff = reformulation.objective - original.objective;
        double timeDiff = reformulation.solveTimeSec - original.solveTimeSec;

        String betterObj;
        if (Math.abs(objDiff) <= 1e-6) {
            betterObj = "Tie";
        } else {
            betterObj = (objDiff < 0.0) ? reformulation.modelName : original.modelName;
        }

        String faster;
        if (Math.abs(timeDiff) <= 1e-6) {
            faster = "Tie";
        } else {
            faster = (timeDiff < 0.0) ? reformulation.modelName : original.modelName;
        }

        System.out.println(
                "Comparison: betterObj=" + betterObj
                        + ", objDelta(reform-original)=" + format(objDiff)
                        + ", faster=" + faster
                        + ", timeDeltaSec(reform-original)=" + format(timeDiff)
        );
    }

    private static boolean hasObjective(SolveResult r) {
        return r.feasible && !Double.isNaN(r.objective);
    }

    private static void printValidation(String label, SolveResult baseline, SolveResult candidate, double tol) {
        if (!hasObjective(baseline) || !hasObjective(candidate)) {
            System.out.println("Validation(" + label + "): SKIP | missing feasible objective.");
            return;
        }

        double diff = candidate.objective - baseline.objective;
        double absDiff = Math.abs(diff);
        String verdict = absDiff <= tol ? "PASS" : "FAIL";

        System.out.println(
                "Validation(" + label + "): " + verdict
                        + " | tol=" + format(tol)
                        + " | objDelta(candidate-baseline)=" + format(diff)
                        + " | absDelta=" + format(absDiff)
        );
    }

    private static String format(double v) {
        return String.format(Locale.US, "%.6f", v);
    }
}
