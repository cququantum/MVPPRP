package model;

public final class CplexConfig {
    public static final double TIME_LIMIT_SEC = 300.0;
    public static final double MIP_GAP = 1e-4;
    public static final boolean LOG_TO_CONSOLE = true;

    private CplexConfig() {
    }
}
