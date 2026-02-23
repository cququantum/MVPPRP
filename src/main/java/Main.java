import instance.Instance;
import model.SolveResult;
import originalModel.OriginalModelSolver;
import reformulationModel.ReformulationModelSolver;

import java.io.IOException;
import java.util.Locale;

public class Main {
    public static void main(String[] args) {
        String instancePath = (args.length > 0) ? args[0] : "data/MVPRP/MVPRP1_10_6_2.txt";

        Instance.Options options = Instance.Options.defaults();
        options.distanceMode = Instance.Options.DistanceMode.EUCLIDEAN_FLOAT;
        options.autoSetDt = true;

        try {
            Instance ins = Instance.fromFile(instancePath, options);

            SolveResult originalResult = new OriginalModelSolver().solve(ins);
            SolveResult reformulationResult = new ReformulationModelSolver().solve(ins);

            System.out.println("Instance: " + instancePath);
            System.out.println("n=" + ins.n + ", l=" + ins.l + ", K=" + ins.K + ", Q=" + format(ins.Q));
            System.out.println(originalResult.toSummaryLine());
            System.out.println(reformulationResult.toSummaryLine());
            printComparison(originalResult, reformulationResult);
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

    private static String format(double v) {
        return String.format(Locale.US, "%.6f", v);
    }
}
