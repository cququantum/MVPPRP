package model;

public final class CplexConfig {
    public static final double TIME_LIMIT_SEC = 6000.0;
    public static final double MIP_GAP = 1e-6;
    public static final boolean LOG_TO_CONSOLE = true;

    private CplexConfig() {
    }
}
