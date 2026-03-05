# speed.tex 完整集成计划（分步实现 + 逐步验证）

> **验证方法**：每完成一步后运行 `--solver=lbbd` 测试实例（如 `MVPRP1_10_3_2.txt`），确认：
> 1. 最终目标值与 `--solver=reform` 一致（tolerance 1e-4 内）
> 2. 不出现新的 infeasible/crash
> 3. 观察迭代次数变化趋势

---

## 全局分析：现有代码结构要点

### MasterModel 字段（其它方法可访问）
```
omega[l+1], m[l+1], lambda[n+1][l+2][l+2], routingLambdaRefs
```

### MasterModel.build() 局部变量（仅 build 内可用）
```
y[t], p[t], p0[t], i0[t]
```

### 关键事实
- **BMP 中无显式 `z_{it}` 变量**：`z_{it} = sum_{v=pi(i,t)}^{t-1} lambda_{ivt}`，所有涉及 z 的约束必须用 lambda 和来表达。
- **现有 VehicleCap**：`sum g*lambda <= Q*m[t]` + `m[t] <= K` 在 LP 松弛中已蕴含 A.2，但显式添加 `sum g*lambda <= KQ` 可帮助 CPLEX 更快传播界。
- **现有代码在 build() 中已跳过 `g > Q` 的弧**（不创建 lambda 变量），A.3 只需处理 `Q/2 < g <= Q`。
- **T2 是独立 LP**：因为 T2 需要自己的 y/p/P0/I0 等连续变量 + 路由列变量 + 跨期约束，应作为独立的 `T2InitialCutSolver` 类实现，不复用 MasterModel。

---

## Step 1：有效不等式 A.1 — 最小收集启动时间约束

### 修改文件
`LbbdReformulationSolver.java` → `MasterModel.build()` 末尾

### 做什么
计算 `tPrime`：满足 `P00 + I00/k < sum_{tau=1}^{t} d_tau` 的最小 t。若存在，添加：
```
sum_{tau=1}^{tPrime} sum_{i=1}^{n} sum_{v=pi(i,tau)}^{tau-1} lambda[i][v][tau] >= 1
```

### 代码骨架
```java
// === Valid Inequality A.1: Minimum collection start time (eq 4-1/4-2) ===
{
    double initialCap = ins.P00 + ins.I00 / ins.k;
    double cumDemand = 0.0;
    int tPrime = -1;
    for (int t = 1; t <= l; t++) {
        cumDemand += ins.dt[t];
        if (initialCap < cumDemand - 1e-9) {  // 严格 < 
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
```

### 正确性检查
1. **严格不等号**：条件是 `P00 + I00/k < cumDemand`，不是 `<=`。如果等于，初始库存刚好够，不需要提前收集。用 `initialCap < cumDemand - 1e-9` 来实现严格 `<`（考虑浮点误差）。
2. **z 用 lambda 替代**：公式中的 `z_{i,tau}` 等价于 `sum_{v} lambda_{i,v,tau}`，上面代码直接枚举所有 (i,v,tau) 的 lambda，正确。
3. **边界情况**：若 `ins.k == 0` 则 `I00/k` 为无穷，此时 initialCap 为无穷，不存在 tPrime，不添加约束——安全。但 Instance 默认 k=1.0，实际不会出现 k=0。
4. **若所有 `ins.dt[t]` 为 0**：cumDemand 始终为 0，initialCap >= 0，不存在 tPrime——安全。

### 验证
```bash
mvn compile exec:java -Dexec.mainClass="Main" -Dexec.args="--solver=lbbd data/MVPRP/MVPRP1_10_3_2.txt"
# 预期：目标值不变（Validation PASS），迭代次数可能微减
```

---

## Step 2：有效不等式 A.2 — 周期车辆总容量约束

### 修改文件
同上，`MasterModel.build()` 末尾（Step 1 之后）

### 做什么
对每个 `t ∈ {1..l}`，添加：
```
sum_{i=1}^{n} sum_{v=pi(i,t)}^{t-1} g(i,v,t) * lambda[i][v][t] <= K * Q
```

### 代码骨架
```java
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
```

### 正确性检查
1. **RHS 必须是 `K * Q`**（不是 `Q` 也不是 `K`）。
2. **不会切掉可行解**：在任何可行的整数解中，所有收集量都必须被 K 辆容量 Q 的车装走，总量自然 <= KQ。

### 验证
运行对比目标值不变。

---

## Step 3：有效不等式 A.3 — 大收集量客户约束

### 修改文件
同上，`MasterModel.build()` 末尾（Step 2 之后）

### 做什么
对每个 `t ∈ {1..l}`，添加：
```
sum_{(i,v): g(i,v,t) > Q/2} lambda[i][v][t] <= K
```

### 代码骨架
```java
// === Valid Inequality A.3: Large-item constraint (eq 4-4) ===
for (int t = 1; t <= l; t++) {
    IloLinearNumExpr expr = cplex.linearNumExpr();
    boolean hasTerms = false;
    for (int i = 1; i <= n; i++) {
        for (int v = ins.pi[i][t]; v <= t - 1; v++) {
            if (lambda[i][v][t] != null) {
                double gVal = ins.g(i, v, t);
                if (gVal > ins.Q / 2.0 + 1e-9) {  // 严格 >
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
```

### 正确性检查
1. **严格大于**：条件是 `g > Q/2`，用 `> Q/2 + epsilon` 实现。
2. **现有代码已不为 `g > Q` 创建 lambda 变量**，所以这里只处理 `Q/2 < g <= Q` 的弧。
3. **hasTerms 检查**：如果某期没有大客户弧，跳过。

### 验证
```bash
mvn compile exec:java -Dexec.mainClass="Main" -Dexec.args="--solver=lbbd data/MVPRP/MVPRP1_10_3_2.txt"
# 三个有效不等式全加后对比：目标值不变，期望迭代次数有所减少
```

---

## Step 4：新建 T1InitialCutSolver

### 新建文件
`src/main/java/lbbdModel/rmp/T1InitialCutSolver.java`

### T1 与现有 PeriodRouteMasterLpSolver 的精确差异

| 方面 | PeriodRouteMasterLpSolver | T1InitialCutSolver |
|------|--------------------------|-------------------|
| 覆盖约束 RHS | 活跃弧 = 1.0，其余 = 0.0 | **全部 = 1.0** |
| 约束类型 | `addEq`（=） | **`addGe`（>=）** |
| 对偶 w 符号 | 无约束（可正可负） | **w >= 0**（由 >= 保证） |
| 依赖 BMP 解 | 是（需要 prevVisitBySupplier） | **否**（全局有效） |
| 何时调用 | LBBD 迭代中按需调用 | **迭代前一次性调用** |

### 接口设计
```java
package lbbdModel.rmp;

public final class T1InitialCutSolver {
    
    public static final class T1Result {
        public final boolean feasible;
        public final boolean optimal;
        public final boolean pricingProvedOptimal;
        public final boolean artificialClean;
        public final double lpObjective;
        public final double[][] dualW;    // [n+1][l+2]: dualW[i][v]
        public final double dualU0;       // 车辆限制对偶 (<= 0)
        public final int generatedColumns;
    }
    
    public T1InitialCutSolver(Instance ins) { ... }
    
    /** 对 period t 求解 T1 LP，返回对偶变量 */
    public T1Result solve(int t) { ... }
}
```

### 核心求解逻辑

```
1. 枚举情景：
   对所有 i ∈ {1..n}, v ∈ {pi(i,t)..t-1}:
     如果 g(i,v,t) <= Q: 创建情景 (supplier=i, prev=v, demand=g(i,v,t))
   → 情景数组: scenarioSupplier[], scenarioPrev[], scenarioDemand[]

2. 所有情景 RHS 均为 1.0 (区别于 PeriodRouteMasterLpSolver 的选择性 0/1)

3. 建 LP:
   - 目标: min sum b_r * xi_r
   - 覆盖约束: sum_r a_{s,r} * xi_r >= 1.0 对每个情景 s  ← addGe, 不是 addEq
   - 车辆约束: sum_r xi_r <= K
   - 人工变量: a_s >= 0 带大罚项, 加在覆盖约束上
   
4. 初始列: 每个情景的单客户路线 (depot → supplier_s → depot)

5. 列生成循环:
   a. 解 LP
   b. 提取对偶: dualByScenario[s] 从覆盖约束, dualU0 从车辆约束
   c. 调用 PricingEspprcSolver.findNegativeScenarioRoutes() 找负 RC 路线
      → 复用现有定价器, 无需改动
   d. 添加新列, 重复直到无负 RC 列

6. 收敛后:
   a. 检查人工变量 = 0
   b. 检查 dualU0 <= 1e-7
   c. 将 dualByScenario 映射回 dualW[i][v]:
      for s in scenarios: dualW[scenarioSupplier[s]][scenarioPrev[s]] = dual[s]
```

### 正确性关键点

1. **列生成必须完全收敛**：pricing 确认 `foundNegativeColumn == false` 才能取对偶。未收敛的对偶可能不满足 DualCVRP(t) 的全局对偶约束，导致割切掉最优解（参见 lbbd.tex 注释 3）。

2. **>= 约束保证 w >= 0**：对偶理论——最小化问题中 `>= b` 约束的对偶变量非负。

3. **情景中同一 supplier i 可有多个 v 值**：例如 supplier 3 可能有 (3,0,t), (3,1,t), (3,2,t) 三个情景。一条路线只能覆盖同一 supplier 的至多一个情景。`PricingEspprcSolver` 的 `ScenarioLongSearch` 已通过 `supplierMask` 处理此约束——复用时无需改动。

4. **T1 的 RHS 全 1 意味着什么**：T1 要求每个情景 (i,v) 至少被路线覆盖一次。在实际路由中，对于给定的 lambda 配置，只有 lambda_{ivt}=1 的 (i,v) 才需要被覆盖。T1 的"过度覆盖"使得 T1 的最优值通常 >= 实际 CVRP(t) 最优值（因为要服务更多客户），但其对偶因 w>=0 而给出全局有效的下界。

### 验证（独立验证 T1 求解器，暂不接入 LBBD）
在 T1 的 solve() 返回前打印：
```
[T1] period=t, lpObj=..., dualU0=..., artificialClean=..., cols=...
```
确认每一期：lpObj >= 0，dualU0 <= 0，artificialClean == true，所有 dualW[i][v] >= -1e-9。

---

## Step 5：MasterModel 添加初始割方法 + 接入 T1 到 LBBD

### 修改文件
`LbbdReformulationSolver.java`

### 5a. MasterModel 中添加方法

```java
/**
 * 添加对偶初始割: omega[t] >= sum_{i,v} w[i][v] * lambda[i][v][t] + K * u0
 * 此方法同时适用于 T1 和 T2 的按期割（数学形式相同，只是 w/u0 数值不同）。
 */
IloRange addDualInitialCut(int t, double[][] dualW, double dualU0, String tag) 
        throws IloException {
    IloLinearNumExpr expr = cplex.linearNumExpr();
    expr.addTerm(1.0, omega[t]);
    double rhs = ins.K * dualU0;

    for (int i = 1; i <= ins.n; i++) {
        for (int v = ins.pi[i][t]; v <= t - 1; v++) {
            if (lambda[i][v][t] == null) continue;
            double w = 0.0;
            if (dualW != null && i < dualW.length && v < dualW[i].length) {
                w = dualW[i][v];
            }
            if (Math.abs(w) > 1e-9) {
                expr.addTerm(-w, lambda[i][v][t]);
            }
        }
    }
    return cplex.addGe(expr, rhs, tag);
}
```

**公式推导验证**：
- 割：`omega[t] >= sum_{i,v} w*_{iv} * lambda_{ivt} + K * u0*`
- 移项：`omega[t] - sum_{i,v} w*_{iv} * lambda_{ivt} >= K * u0*`
- 代码：`expr = omega[t] + sum(-w * lambda)`，`rhs = K * u0`
- `u0 <= 0` → `rhs <= 0` → 不等式方向合理
- 当所有 lambda=0 时：`omega[t] >= K*u0 <= 0`，这是平凡的（omega >= 0 已知），不会切掉可行解
- 当某些 lambda=1 时：`omega[t] >= sum(active w) + K*u0`，这正是弱对偶给出的 CVRP(t) 下界 ✓

### 5b. solve() 方法中调用 T1

在 `LbbdReformulationSolver.solve()` 中，位于 `lbRResult` 打印之后、`while (iteration < MAX_ITERATIONS)` 之前：

```java
// === T1 Initial Cuts (speed.tex Section B.2, eq 4-9) ===
T1InitialCutSolver t1Solver = new T1InitialCutSolver(ins);
int t1CutsAdded = 0;
for (int t = 1; t <= ins.l; t++) {
    T1InitialCutSolver.T1Result t1 = t1Solver.solve(t);
    if (t1.feasible && t1.optimal && t1.pricingProvedOptimal 
            && t1.artificialClean && t1.dualU0 <= 1e-7) {
        master.addDualInitialCut(t, t1.dualW, t1.dualU0, "T1_InitCut_t" + t);
        t1CutsAdded++;
    }
}
System.out.println("[LBBD] T1 initial cuts added: " + t1CutsAdded + "/" + ins.l);
```

**注意**：
- T1 割永远不需要被删除（全局有效），不加入 `activeRmpDualCuts` 等管理队列。
- T1 求解是一次性的，不在主循环中重复。
- 若某期 T1 求解失败，仅跳过该期（不影响正确性）。

### 验证
```bash
mvn compile exec:java -Dexec.mainClass="Main" -Dexec.args="--solver=lbbd data/MVPRP/MVPRP1_10_3_2.txt"
# 预期：
# 1. Validation(reform-vs-lbbd): PASS
# 2. [LBBD] T1 initial cuts added: l/l (全部成功)
# 3. 首几次迭代的 delta(phi-omega) 应比无割时更小（omega 不再从 0 开始）
# 4. 总迭代次数减少
```

---

## Step 6：新建 T2InitialCutSolver

### 新建文件
`src/main/java/lbbdModel/rmp/T2InitialCutSolver.java`

### T2 的本质
T2 是**改写模型（reformulationModel.tex）的完整 LP 松弛**，其中路由部分用集合划分列表示。T2 **跨期耦合**（通过生产/库存平衡），因此比 T1（单期独立）给出更紧的对偶下界。

### 接口设计

```java
package lbbdModel.rmp;

public final class T2InitialCutSolver {

    public static final class T2Result {
        public final boolean feasible;
        public final boolean optimal;
        public final boolean pricingProvedOptimal;
        public final boolean artificialClean;
        public final double lpObjective;
        public final double[][][] dualWByPeriod;  // [l+1][n+1][l+2]
        public final double[] dualU0ByPeriod;     // [l+1]
        public final int totalGeneratedColumns;
    }

    public T2InitialCutSolver(Instance ins) { ... }
    
    /** 求解全局 T2 LP，返回各期的对偶变量 */
    public T2Result solve() { ... }
}
```

### T2 模型完整结构

#### 决策变量（全部连续）
```
y[t] ∈ [0,1]              t=1..l      生产开关
p[t] >= 0                 t=1..l      生产量
P0[t] >= 0                t=1..l      成品库存
I0[t] >= 0                t=1..l      原料库存
lambda[i][v][t] ∈ [0,1]   有效(i,v,t) 访问弧（t=1..l+1）
xi[r][t] >= 0              列生成添加  路由列变量（t=1..l）
```

**注意：此处 lambda 是 T2 内部 LP 的连续变量，不是 MasterModel 的 lambda（那是整数变量）。T2 完全独立于 MasterModel，在自己的 IloCplex 实例中建模。**

#### 目标函数（对照 reformulationModel.tex 2.4 节）
```
min  sum_{t=1}^{l} [ sum_r b_r * xi_{rt}             // 路由成本
                    + u * p[t] + f * y[t]               // 生产成本
                    + h0 * I0[t] + hp * P0[t] ]         // 工厂/成品库存
   + sum_{t=1}^{l+1} sum_i sum_v e(i,v,t) * lambda[i][v][t]  // 供应商库存
```

#### 约束清单（与 reformulationModel.tex 逐条对应）

**① 成品平衡（2-37）**: `P0[t-1] + p[t] = d[t] + P0[t]`，t=1 时 `P0[0]=P00`
**② 生产能力（2-38）**: `p[t] <= C * y[t]`
**③ 成品容量（2-39）**: `P0[t] <= Lp`
**④ 工厂原料平衡（2-40）**: `I0[t-1] + sum_{i,v} g(i,v,t)*lambda[i][v][t] = k*p[t] + I0[t]`，t=1 时 `I0[0]=I00`
**⑤ 工厂容量（2-41）**: `I0[t] <= L0`
**⑥ 起始流（2-43）**: `sum_{t=1}^{mu(i,0)} lambda[i][0][t] = 1`
**⑦ 流平衡（2-44）**: `sum_v lambda[i][v][t] - sum_tau lambda[i][t][tau] = 0`，对 i∈V_S, t∈T
**⑧ 终止流（2-45）**: `sum_{t=pi(i,l+1)}^{l} lambda[i][t][l+1] = 1`
**⑨★ 路由覆盖（生成 w\*）**: `sum_r a_{ivr}*xi[r][t] = lambda[i][v][t]`，对每个有效 (i,v) 且 t∈T
**⑩★ 车辆限制（生成 u0\*）**: `sum_r xi[r][t] <= K`，对每个 t∈T

**T2 中不含的约束**（与弧变量相关的，已被集合划分列替代）：
- 无 MTZ 子回路消除
- 无弧变量 x_{ijt}
- 无 m[t] 变量（车辆数由 sum xi_{rt} 隐含）
- 无 VehicleCap: `sum g*lambda <= Q*m[t]`（容量由列的可行性保证）
- 无 z_{it} 变量（由 lambda 流约束隐含）
- 无 VisitLink: `sum x = z`（无 x 变量）

### 覆盖约束的建模细节（关键）

覆盖约束 ⑨ 将 lambda 与 xi 链接：
```
sum_r a_{ivr} * xi_{rt} = lambda_{ivt}
```

建模时，lambda 是 T2 LP 的变量（系数 -1），xi 通过列生成添加（系数 +1）：
```java
// 建模: sum_r a*xi - lambda = 0
IloLinearNumExpr lhs = cplex.linearNumExpr();
lhs.addTerm(-1.0, lambdaVar[i][v][t]);
coverEq[t][scenarioIndex] = cplex.addEq(lhs, 0.0, name);

// 列生成添加 xi 列时:
IloColumn col = cplex.column(obj, routeCost)
    .and(cplex.column(vehicleLimit[t], 1.0));
for (int s : coveredScenarios) {
    col = col.and(cplex.column(coverEq[t][s], 1.0));
}
cplex.numVar(col, 0.0, Double.MAX_VALUE, varName);
```

### 列生成流程

```
初始化:
  1. 建所有非路由变量和约束 (①~⑧)
  2. 对每个 period t, 枚举情景 (i,v), 建覆盖约束 ⑨ 和车辆限制 ⑩
  3. 对每个 period t 的每个情景, 添加单客户初始列
  4. 对每个 period t 的每个覆盖约束, 添加人工变量（大罚项）

列生成循环:
  while (true):
    a. 解 T2 主 LP
    b. foundAnyNewColumn = false
    c. for t = 1..l:
         提取 dualByScenario_t[] 从 coverEq[t][*]
         提取 dualU0_t 从 vehicleLimit[t]
         调用 PricingEspprcSolver.findNegativeScenarioRoutes(
             ins, scenarioSupplier_t, scenarioDemand_t,
             dualByScenario_t, dualU0_t,
             existingRouteKeys_t, maxCols
         )
         如果找到负 RC 列: 添加列, foundAnyNewColumn = true
    d. if (!foundAnyNewColumn): break

提取结果:
  for t = 1..l:
    dualU0ByPeriod[t] = cplex.getDual(vehicleLimit[t])
    for each scenario s of period t:
      dualWByPeriod[t][supplier[s]][prev[s]] = cplex.getDual(coverEq[t][s])
```

### 正确性关键点

1. **为什么 T2 的 (w\*_{ivt}, u0t\*) 是 DualCVRP(t) 的可行解？**
   `xi[r][t]` 仅出现在两类约束中：覆盖约束（系数 a_{ivr}）和车辆约束（系数 1）。
   由 LP 对偶理论，`xi[r][t]` 的对偶约束为：`sum_{(i,v)} a_{ivr} * w_{ivt} + u0t <= b_r`，对所有路线 r。
   这与 DualCVRP(t) 的对偶约束**形式完全相同**。✓

2. **T2 的 w 可以为负**：因为覆盖约束是等式（= lambda），对偶 w 无符号约束。与 T1 不同。但割仍然安全。

3. **列生成必须完全收敛**（同 T1）。

4. **Lambda 的 t=l+1 弧**：lambda[i][v][l+1] 出现在终止流约束和目标函数（e 成本），但**不出现在路由覆盖约束中**（t=l+1 不是实际路由期）。覆盖约束只对 t=1..l 建立。

5. **lambda 在 T2 中是连续的 [0,1]**：这是 LP 松弛。得到的对偶解用于割，割在 BMP 中对整数 lambda 也成立（因为 DualCVRP 可行性不依赖 lambda 的值）。

6. **T2 比 T1 更紧的原因**：T2 的覆盖约束是 `= lambda`（RHS 由 LP 优化决定），而 T1 是 `>= 1`（固定 RHS）。T2 还包含跨期耦合约束（生产/库存），使得 lambda 的取值受到全局生产计划的约束，对偶值更精确地反映了路由成本的下界。

### 数值验证方法
对偶割的一致性检查：求解后，对每个 t，计算
```
cutVal_t = sum_{i,v} w*_{ivt} * lambdaSolVal_{ivt} + K * u0t*
```
应当近似等于 T2 LP 中第 t 期的路由成本贡献 `sum_r b_r * xiSolVal_{rt}`。若不一致，说明对偶提取有符号错误。

### 验证（独立验证 T2 求解器）
打印：
```
[T2] lpObj=..., totalCols=..., artificialClean=..., converged=...
[T2] period t: dualU0=..., maxW=..., minW=...
```
确认：
- lpObj >= 0 且 lpObj <= reform 最优值（T2 是 LP 松弛，应 <= MIP 最优值）
- 每期 dualU0[t] <= 0
- artificialClean == true

---

## Step 7：接入 T2 初始割到 LBBD

### 修改文件
`LbbdReformulationSolver.java` → `solve()` 方法

### 在 T1 初始割之后添加

```java
// === T2 Initial Cuts (speed.tex Section B.3, eq 4-12, per-period form) ===
T2InitialCutSolver t2Solver = new T2InitialCutSolver(ins);
T2InitialCutSolver.T2Result t2 = t2Solver.solve();
int t2CutsAdded = 0;
if (t2.feasible && t2.optimal && t2.pricingProvedOptimal && t2.artificialClean) {
    for (int t = 1; t <= ins.l; t++) {
        if (t2.dualU0ByPeriod[t] <= 1e-7) {
            master.addDualInitialCut(t, t2.dualWByPeriod[t], 
                    t2.dualU0ByPeriod[t], "T2_InitCut_t" + t);
            t2CutsAdded++;
        }
    }
}
System.out.println("[LBBD] T2 initial cuts added: " + t2CutsAdded + "/" + ins.l);
```

### 关于割的形式选择

speed.tex 给出的原始形式（4-12）是**总和割**（一条约束跨所有期）：
```
sum_t omega[t] >= sum_t (sum_{i,v} w*_{ivt} * lambda_{ivt} + K * u0t*)
```

上面的代码用的是**按期割**（每期一条，共 l 条）：
```
omega[t] >= sum_{i,v} w*_{ivt} * lambda_{ivt} + K * u0t*,  for each t
```

**按期割严格强于总和割**：按期割的总和蕴含总和割，但反之不成立。speed.tex 也提到"可以拆成每期一条（更紧）"。使用按期割是更好的选择。

### T1 与 T2 初始割共存

- 同一期 t 会有 T1 割和 T2 割，两者系数不同但形式相同。
- 两者同时存在是安全的（都是有效不等式），CPLEX 会自动处理。
- T2 割通常更紧（因为包含跨期信息），但 T1 割的 w >= 0 性质有时在某些 lambda 配置下更有效。两者互补。
- **不要因为有 T2 就跳过 T1**：T1 计算快且稳定，是安全的保底。

### 验证
```bash
# 单实例测试
mvn compile exec:java -Dexec.mainClass="Main" -Dexec.args="--solver=lbbd data/MVPRP/MVPRP1_10_3_2.txt"

# 多实例测试
mvn compile exec:java -Dexec.mainClass="Main" \
  -Dexec.args="data/MVPRP/MVPRP1_10_3_2.txt data/MVPRP/MVPRP2_10_3_2.txt data/MVPRP/MVPRP3_10_3_2.txt data/MVPRP/MVPRP4_10_3_2.txt"

# 预期:
# 1. 所有 Validation(reform-vs-lbbd): PASS
# 2. [LBBD] T1 initial cuts added: l/l
# 3. [LBBD] T2 initial cuts added: l/l
# 4. 迭代次数进一步减少（相比仅 T1）
# 5. 首次迭代的 delta(phi-omega) 应明显更小
```

---

## 完整文件清单

### 新建文件
1. `src/main/java/lbbdModel/rmp/T1InitialCutSolver.java`
2. `src/main/java/lbbdModel/rmp/T2InitialCutSolver.java`

### 修改文件
3. `src/main/java/lbbdModel/LbbdReformulationSolver.java`
   - `MasterModel.build()` 末尾：添加 A.1 / A.2 / A.3 有效不等式
   - `MasterModel` 类：添加 `addDualInitialCut()` 方法
   - `solve()` 方法：在主循环前依次调用 T1 和 T2 初始割

### 不修改的文件
- `PeriodRouteMasterLpSolver.java`、`PricingEspprcSolver.java`、`RoutingLowerBoundSolver.java`
- `Instance.java`、`Main.java`
- 所有 `originalModel/` 和 `reformulationModel/` 代码

---

## 各步预期效果与失败信号

| 步骤 | 新增内容 | 正常预期 | 失败信号 |
|------|---------|---------|---------|
| Step 1 | A.1 约束 | obj 不变，LP 微收紧 | obj 变大 → 约束过强/错误 |
| Step 2 | A.2 约束 | obj 不变，LP 收紧 | obj 变大 → RHS 错误 |
| Step 3 | A.3 约束 | obj 不变，LP 微收紧 | obj 变大 → 条件判断错误 |
| Step 4 | T1 求解器 | 独立可运行，w>=0 | T1 不可行 / w<0 / 未收敛 |
| Step 5 | T1 接入 | 迭代数减少，obj 不变 | obj 偏移 → 割方向/系数错 |
| Step 6 | T2 求解器 | 独立可运行，lpObj <= reform | 人工变量残留 / 不收敛 |
| Step 7 | T2 接入 | 迭代数进一步减少 | obj 偏移 → 对偶提取符号错 |

---

## 附录 A：T2 中 lambda 变量的双重角色

在 T2 LP 中，`lambda[i][v][t]` 同时出现在：

| 约束 | lambda 系数 | 约束中的对偶变量 |
|------|------------|----------------|
| 覆盖约束 ⑨ (t∈T) | -1 | w_{ivt} ← **我们需要的** |
| 工厂原料平衡 ④ | +g(i,v,t) | 工厂平衡对偶 |
| 流平衡 ⑦ (入弧) | +1 | 流对偶 |
| 流平衡 ⑦ (出弧) | -1 | 流对偶 |
| 起始流 ⑥ / 终止流 ⑧ | +1 | 端点对偶 |
| 目标函数 | e(i,v,t) | — |

T2 的覆盖约束对偶 w\*_{ivt} 已经通过 LP 对偶性自动包含了所有交互效应——这是 T2 比 T1 更紧的根本原因。

但在提取割系数时，我们**只需要覆盖约束 ⑨ 和车辆约束 ⑩ 的对偶**（因为 xi 只出现在这两类约束中）。

## 附录 B：对偶提取的符号验证方法

实现完 T2 后，可用以下方法验证对偶符号是否正确：

```java
// 求解后，对每个 period t：
double routingCostT = 0.0; // sum_r b_r * xiVal_r (已知)
double cutValT = dualU0[t] * ins.K;
for (int s = 0; s < scenarioCount_t; s++) {
    cutValT += cplex.getDual(coverEq[t][s]) * cplex.getValue(lambdaVar[...]);
}
// routingCostT 应近似等于 cutValT（互补松弛）
// 如果差距大 → 符号有问题
```

## 附录 C：T1 复用 PeriodRouteMasterLpSolver 代码的指南

T1InitialCutSolver 的核心逻辑与 PeriodRouteMasterLpSolver.solve() 高度相似。建议的复用策略：

**直接复制 PeriodRouteMasterLpSolver.solve() 后做以下修改**：

1. **情景 RHS**：
   ```java
   // 原: rhsList.add((vBar == v) ? 1.0 : 0.0);
   // 改: 
   rhsList.add(1.0);  // T1: 所有情景 RHS 均为 1.0
   ```

2. **约束类型**：
   ```java
   // 原: coverEq[s] = cplex.addEq(..., scenarioRhs[s], ...);
   // 改:
   coverEq[s] = cplex.addGe(..., scenarioRhs[s], ...);  // T1: >= 而非 =
   ```

3. **移除 prevVisitBySupplier 参数**：T1 不依赖当前 BMP 解。

4. **返回类型**：改为 T1Result，dualW 映射到 [i][v] 而非 [i][prevVisit]。

5. **其余代码不变**：初始列、人工变量、列生成循环、定价器调用全部复用。
