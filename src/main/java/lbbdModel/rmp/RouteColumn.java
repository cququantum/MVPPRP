package lbbdModel.rmp;

import java.util.Arrays;

public final class RouteColumn {
    public final int[] globalCustomersInOrder;
    public final double cost;
    public final double load;

    public RouteColumn(int[] globalCustomersInOrder, double cost, double load) {
        this.globalCustomersInOrder = globalCustomersInOrder;
        this.cost = cost;
        this.load = load;
    }

    public String key() {
        if (globalCustomersInOrder.length == 0) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int k = 0; k < globalCustomersInOrder.length; k++) {
            if (k > 0) {
                sb.append('-');
            }
            sb.append(globalCustomersInOrder[k]);
        }
        return sb.toString();
    }

    @Override
    public String toString() {
        return "RouteColumn{cost=" + cost + ", load=" + load + ", path=" + Arrays.toString(globalCustomersInOrder) + "}";
    }
}
