package lbbdModel.rmp;

import java.util.Arrays;

public final class ScenarioRouteColumn {
    public final int[] scenarioIndicesInOrder;
    public final double cost;
    public final double load;

    public ScenarioRouteColumn(int[] scenarioIndicesInOrder, double cost, double load) {
        this.scenarioIndicesInOrder = scenarioIndicesInOrder;
        this.cost = cost;
        this.load = load;
    }

    public String key() {
        if (scenarioIndicesInOrder.length == 0) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int k = 0; k < scenarioIndicesInOrder.length; k++) {
            if (k > 0) {
                sb.append('-');
            }
            sb.append('s').append(scenarioIndicesInOrder[k]);
        }
        return sb.toString();
    }

    @Override
    public String toString() {
        return "ScenarioRouteColumn{cost=" + cost + ", load=" + load + ", scenarios="
                + Arrays.toString(scenarioIndicesInOrder) + "}";
    }
}
