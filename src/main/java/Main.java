import instance.Instance;
import lbbdModel.LbbdReformulationSolver;
import model.SolveResult;
import originalModel.OriginalModelSolver;
import reformulationModel.ReformulationModelSolver;

import java.io.IOException;
import java.util.Locale;

public class Main {
    public static void main(String[] args) {
        String instancePath = (args.length > 0) ? args[0] : "data/MVPRP/MVPRP2_10_3_2.txt";

        Instance.Options options = Instance.Options.defaults();
        options.distanceMode = Instance.Options.DistanceMode.EUCLIDEAN_FLOAT;
        options.autoSetDt = true;

        try {
            Instance ins = Instance.fromFile(instancePath, options);

            SolveResult originalResult = new OriginalModelSolver().solve(ins);
            SolveResult reformulationResult = new ReformulationModelSolver().solve(ins);
            SolveResult lbbdResult = new LbbdReformulationSolver().solve(ins);

            System.out.println("Instance: " + instancePath);
            System.out.println("n=" + ins.n + ", l=" + ins.l + ", K=" + ins.K + ", Q=" + format(ins.Q));
            System.out.println(originalResult.toSummaryLine());
            System.out.println(reformulationResult.toSummaryLine());
            System.out.println(lbbdResult.toSummaryLine());
            printComparison(originalResult, reformulationResult);
            printValidation("reform-vs-lbbd", reformulationResult, lbbdResult, 1e-4);
        } catch (IOException e) {
            throw new RuntimeException("Failed to load instance file: " + instancePath, e);
        }
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
