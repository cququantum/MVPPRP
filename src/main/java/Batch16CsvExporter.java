import instance.Instance;
import lbbdModel.LbbdReformulationSolver;
import model.CplexConfig;
import model.SolveResult;
import originalModel.OriginalModelSolver;
import reformulationModel.ReformulationModelSolver;

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
import java.util.List;
import java.util.Locale;
import java.util.Set;

public final class Batch16CsvExporter {
    private static final double RESULT_TOL = CplexConfig.MIP_GAP;
    private static final int METHODS_PER_INSTANCE = 4;
    private static final String CSV_HEADER =
            "instance,method,status,feasible,optimal,objective,best_bound,gap,time_sec";
    private static final OutputStream DEV_NULL = new OutputStream() {
        @Override
        public void write(int b) {
            // discard
        }
    };
    private static final String[] INSTANCE_PATHS = new String[]{
            "data/MVPRP/MVPRP1_10_3_2.txt",
            "data/MVPRP/MVPRP1_10_3_3.txt",
            "data/MVPRP/MVPRP1_10_6_2.txt",
            "data/MVPRP/MVPRP1_10_6_3.txt",
            "data/MVPRP/MVPRP1_10_9_2.txt",
            "data/MVPRP/MVPRP1_10_9_3.txt",
            "data/MVPRP/MVPRP2_10_3_2.txt",
            "data/MVPRP/MVPRP2_10_3_3.txt",
            "data/MVPRP/MVPRP2_10_6_2.txt",
            "data/MVPRP/MVPRP2_10_6_3.txt",
            "data/MVPRP/MVPRP2_10_9_2.txt",
            "data/MVPRP/MVPRP2_10_9_3.txt",
            "data/MVPRP/MVPRP3_10_3_2.txt",
            "data/MVPRP/MVPRP3_10_3_3.txt",
            "data/MVPRP/MVPRP3_10_6_2.txt",
            "data/MVPRP/MVPRP3_10_6_3.txt",
            "data/MVPRP/MVPRP3_10_9_2.txt",
            "data/MVPRP/MVPRP3_10_9_3.txt",
            "data/MVPRP/MVPRP4_10_3_2.txt",
            "data/MVPRP/MVPRP4_10_3_3.txt",
            "data/MVPRP/MVPRP4_10_6_2.txt",
            "data/MVPRP/MVPRP4_10_6_3.txt",
            "data/MVPRP/MVPRP4_10_9_2.txt",
            "data/MVPRP/MVPRP4_10_9_3.txt"
    };
    private static final Path OUTPUT_CSV = Paths.get("results_24cases_600s.csv");

    private Batch16CsvExporter() {
    }

    public static void main(String[] args) throws Exception {
        Set<String> completedKeys = new HashSet<String>();
        ensureOutputFile(completedKeys);

        try (BufferedWriter writer = Files.newBufferedWriter(
                OUTPUT_CSV,
                StandardCharsets.UTF_8,
                StandardOpenOption.APPEND
        )) {
            int totalRows = INSTANCE_PATHS.length * METHODS_PER_INSTANCE;
            int progress = completedKeys.size();

            for (int idx = 0; idx < INSTANCE_PATHS.length; idx++) {
                String instancePath = INSTANCE_PATHS[idx];
                String instanceName = instanceName(instancePath);
                final Instance ins = loadInstance(instancePath);

                progress = runAndWriteIfNeeded(
                        writer,
                        ins,
                        instanceName,
                        "origin",
                        completedKeys,
                        progress,
                        totalRows,
                        new SolveSupplier() {
                            @Override
                            public SolveResult get() {
                                return new OriginalModelSolver().solve(ins);
                            }
                        }
                );

                progress = runAndWriteIfNeeded(
                        writer,
                        ins,
                        instanceName,
                        "reform",
                        completedKeys,
                        progress,
                        totalRows,
                        new SolveSupplier() {
                            @Override
                            public SolveResult get() {
                                return new ReformulationModelSolver().solve(ins);
                            }
                        }
                );

                progress = runAndWriteIfNeeded(
                        writer,
                        ins,
                        instanceName,
                        "lbbd_no_init",
                        completedKeys,
                        progress,
                        totalRows,
                        new SolveSupplier() {
                            @Override
                            public SolveResult get() {
                                return new LbbdReformulationSolver(ins, false).solve(Double.NaN, RESULT_TOL);
                            }
                        }
                );

                progress = runAndWriteIfNeeded(
                        writer,
                        ins,
                        instanceName,
                        "lbbd",
                        completedKeys,
                        progress,
                        totalRows,
                        new SolveSupplier() {
                            @Override
                            public SolveResult get() {
                                return new LbbdReformulationSolver(ins).solve(Double.NaN, RESULT_TOL);
                            }
                        }
                );
            }
        }
    }

    private static Instance loadInstance(String instancePath) throws IOException {
        Instance.Options options = Instance.Options.defaults();
        options.distanceMode = Instance.Options.DistanceMode.EUCLIDEAN_FLOAT;
        options.autoSetDt = true;
        return Instance.fromFile(instancePath, options);
    }

    private static int runAndWriteIfNeeded(
            BufferedWriter writer,
            Instance ins,
            String instanceName,
            String method,
            Set<String> completedKeys,
            int progress,
            int totalRows,
            SolveSupplier solveSupplier
    ) throws Exception {
        String key = buildKey(instanceName, method);
        if (completedKeys.contains(key)) {
            System.out.println("[" + progress + "/" + totalRows + "] skip: " + instanceName + " " + method);
            return progress;
        }

        SolveResult result = runQuietly(solveSupplier);
        writeRow(writer, instanceName, method, result);
        completedKeys.add(key);
        int newProgress = progress + 1;
        System.out.println("[" + newProgress + "/" + totalRows + "] done: " + instanceName + " " + method);
        return newProgress;
    }

    private static SolveResult runQuietly(SolveSupplier solveSupplier) throws Exception {
        PrintStream originalOut = System.out;
        PrintStream originalErr = System.err;
        PrintStream silent = new PrintStream(DEV_NULL);
        try {
            System.setOut(silent);
            System.setErr(silent);
            return solveSupplier.get();
        } finally {
            System.setOut(originalOut);
            System.setErr(originalErr);
            silent.close();
        }
    }

    private static void ensureOutputFile(Set<String> completedKeys) throws IOException {
        if (!Files.exists(OUTPUT_CSV)) {
            writeHeader(OUTPUT_CSV);
            return;
        }

        try (BufferedReader reader = Files.newBufferedReader(OUTPUT_CSV, StandardCharsets.UTF_8)) {
            String header = reader.readLine();
            if (header == null) {
                writeHeader(OUTPUT_CSV);
                return;
            }
            if (!CSV_HEADER.equals(header)) {
                throw new IllegalStateException("Unexpected CSV header in " + OUTPUT_CSV + ": " + header);
            }

            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    continue;
                }
                List<String> fields = parseCsvLine(line);
                if (fields.size() < 2) {
                    throw new IllegalStateException("Malformed CSV row in " + OUTPUT_CSV + ": " + line);
                }
                completedKeys.add(buildKey(fields.get(0), fields.get(1)));
            }
        }
    }

    private static void writeHeader(Path outputCsv) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(
                outputCsv,
                StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING
        )) {
            writer.write(CSV_HEADER);
            writer.write(System.lineSeparator());
            writer.flush();
        }
    }

    private static void writeRow(BufferedWriter writer, String instanceName, String method, SolveResult result) throws IOException {
        StringBuilder row = new StringBuilder(256);
        appendCsvField(row, instanceName);
        row.append(',');
        appendCsvField(row, method);
        row.append(',');
        appendCsvField(row, safeString(result.status));
        row.append(',');
        row.append(result.feasible);
        row.append(',');
        row.append(result.optimal);
        row.append(',');
        row.append(fmt(result.objective));
        row.append(',');
        row.append(fmt(result.bestBound));
        row.append(',');
        row.append(fmt(result.mipGap));
        row.append(',');
        row.append(fmt(result.solveTimeSec));
        row.append(System.lineSeparator());

        writer.write(row.toString());
        writer.flush();
    }

    private static void appendCsvField(StringBuilder sb, String value) {
        sb.append('"');
        for (int i = 0; i < value.length(); i++) {
            char ch = value.charAt(i);
            if (ch == '"') {
                sb.append('"');
            }
            sb.append(ch);
        }
        sb.append('"');
    }

    private static List<String> parseCsvLine(String line) {
        List<String> fields = new ArrayList<String>();
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

    private static String instanceName(String instancePath) {
        String fileName = Paths.get(instancePath).getFileName().toString();
        if (fileName.endsWith(".txt")) {
            return fileName.substring(0, fileName.length() - 4);
        }
        return fileName;
    }

    private static String buildKey(String instanceName, String method) {
        return instanceName + "|" + method;
    }

    private static String safeString(String s) {
        return s == null ? "" : s;
    }

    private static String fmt(double value) {
        if (Double.isNaN(value)) {
            return "NaN";
        }
        return String.format(Locale.US, "%.6f", value);
    }

    private interface SolveSupplier {
        SolveResult get() throws Exception;
    }
}
