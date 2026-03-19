# 实验补充方案 — 供 Claude Code 执行

本文档描述第五章两组待补充实验的完整方案。请按顺序执行。

---

## 前置信息

- 项目路径：应该在用户的 Java 项目中，入口类和 Instance 加载逻辑参考 `src/main/java/instance/Instance.java`
- 求解器：CPLEX（Java API）
- 已有结果文件：`merged_results.csv`，包含 origin / reform / lbbd / alns 四种方法在 48 个实例上的结果
- LBBD 的实现中，status 字段已经记录了 `feasCuts`, `optCuts`, `lpDualCuts` 等统计量

---

## 实验一：5.3 消融实验（Ablation Study）

### 目标

验证 LBBD 中三项加速策略的增量贡献：
1. 加强割（加强可行性割 + 加强最优性割）
2. 有效不等式（3.6 节的三类）
3. 路径成本下界（用于最优性割中的 θ_t 下界项）

### 实验变体定义

需要运行以下 **4 个 LBBD 变体**（Full-LBBD 已有数据，无需重跑）：

| 变体名称 | 加强割 | 有效不等式 | 路径成本下界 | 说明 |
|---------|:-----:|:---------:|:----------:|------|
| `base`  | ✗ | ✗ | ✗ | 仅使用基础可行性割(3-14)和基础最优性割(3-18) |
| `+cuts` | ✓ | ✗ | ✗ | 替换为加强可行性割(3-16)和加强最优性割(3-19) |
| `+vi`   | ✓ | ✓ | ✗ | 在 +cuts 基础上添加三类有效不等式 |
| `full`  | ✓ | ✓ | ✓ | 完整算法（= 已有的 lbbd 结果，无需重跑） |

**实现方式**：在 LBBD 求解器中添加配置开关（布尔参数），控制上述三项策略的启用/禁用。具体而言：

1. **加强割开关** (`useStrengthenedCuts`):
   - `false`：可行性割使用公式 (3-14)（基础组合割，仅排除当前不可行解），最优性割使用公式 (3-18)（对偶系数仅绑定在当前访问弧 λ_{iv_i,t} 上）
   - `true`：可行性割使用公式 (3-16)（利用单调性加强，将约束扩展到 s ∈ [ε_i(t), v_i] 区间），最优性割使用公式 (3-19)（对偶系数扩展到 s ∈ [ε_i(t), v_i] 的所有可行访问弧）

2. **有效不等式开关** (`useValidInequalities`):
   - `false`：主问题不添加 3.6 节的三类有效不等式
   - `true`：向主问题添加最小收集启动时间约束 (3-29)、周期车辆总容量约束 (3-30)、大收集量客户约束 (3-31)

3. **路径成本下界开关** (`useRoutingLowerBound`):
   - `false`：最优性割中的 θ_t 下界项设为 0（即公式 (3-12) 中 LB_t = 0）
   - `true`：通过求解松弛路由问题 (3-27) 计算路径成本下界 LB_t

### 测试实例

**不需要在全部 48 个实例上跑**。选取 LBBD 优势最显著的子集，即 Full-LBBD 全部最优但 CPLEX 吃力的 16 个实例：

```
# 参数组合 10_6_3 (4个实例)
MVPRP1_10_6_3, MVPRP2_10_6_3, MVPRP3_10_6_3, MVPRP4_10_6_3

# 参数组合 10_9_2 (4个实例)
MVPRP1_10_9_2, MVPRP2_10_9_2, MVPRP3_10_9_2, MVPRP4_10_9_2

# 参数组合 10_9_3 (4个实例)
MVPRP1_10_9_3, MVPRP2_10_9_3, MVPRP3_10_9_3, MVPRP4_10_9_3

# 参数组合 15_3_3 (4个实例)
MVPRP1_15_3_3, MVPRP2_15_3_3, MVPRP3_15_3_3, MVPRP4_15_3_3
```

### 时间限制

每个变体每个实例的全局时间限制设为 **600 秒**（与 full 一致）。

### 输出格式

每个变体运行完后，输出一个 CSV 文件 `ablation_results.csv`，格式如下：

```csv
instance,variant,feasible,optimal,objective,best_bound,gap,time_sec,iterations,feasCuts,optCuts
MVPRP1_10_6_3,base,True,False,45300.12,44800.55,0.011,600.0,8,15,6
MVPRP1_10_6_3,+cuts,True,True,45223.94,45223.94,0.0,12.3,21,2,19
MVPRP1_10_6_3,+vi,True,True,45223.94,45223.94,0.0,8.1,18,2,16
...
```

字段说明：
- `variant`：变体名称（`base` / `+cuts` / `+vi`），`full` 的数据从已有 `merged_results.csv` 中提取即可
- `iterations`：LBBD 迭代次数
- `feasCuts`：添加的可行性割总数
- `optCuts`：添加的最优性割总数

### 实现步骤

1. 先阅读 LBBD 求解器的代码，找到：
   - 可行性割的生成位置（找公式 3-14 / 3-16 对应的代码）
   - 最优性割的生成位置（找公式 3-18 / 3-19 对应的代码）
   - 有效不等式的添加位置（找公式 3-29 / 3-30 / 3-31 对应的代码）
   - 路径成本下界的计算和使用位置（找公式 3-27 / 3-12 对应的代码）

2. 为每项策略添加布尔开关参数，并确保 `false` 时回退到基础版本

3. 写一个运行脚本/main方法，按 3 个变体 × 16 个实例 = 48 次运行，收集结果

4. 如果代码中基础割和加强割已经是用 if-else 或策略模式分开的，直接用开关控制；如果基础割的代码已经被删除只保留了加强版，需要根据论文公式 (3-14) 和 (3-18) 补写基础版本

---

## 实验二：5.4 参数敏感性分析

### 目标

分析关键运营参数对最优决策结构的影响，提炼管理洞察。

### 基准实例选取

从 48 个实例中选取 **2 个代表性实例** 作为基准，要求 LBBD 能在合理时间内求解到最优：

```
MVPRP1_10_6_2  （标准成本结构，中等规模）
MVPRP3_10_6_2  （高运输成本结构，中等规模）
```

如果这两个实例在参数变化后某些方案 LBBD 超时，可改用 `MVPRP1_10_3_2` 和 `MVPRP3_10_3_2`（更小规模，确保可解）。

### 实验 2A：车辆容量 Q 的影响

在基准实例上，保持其他参数不变，将车辆容量 Q 按以下比例缩放：

| 方案 | Q 倍率 | n=10,K=2 时的 Q 值 |
|------|:------:|:-----------------:|
| Q-80%  | 0.80 | 158 |
| Q-100% | 1.00 | 198（基准） |
| Q-120% | 1.20 | 238 |
| Q-150% | 1.50 | 297 |

**实现方式**：在 Instance 加载后、求解前，直接修改 `Q` 的值（乘以倍率后向下取整）。

**求解方法**：使用 LBBD 求解每个方案。时限 600 秒。

**输出格式**：`sensitivity_Q.csv`

```csv
instance,Q_ratio,Q_value,feasible,optimal,objective,time_sec,route_cost,inventory_cost,production_cost,avg_visit_frequency,avg_pickup_per_visit
MVPRP1_10_6_2,0.80,158,True,True,46500.12,25.3,12000.5,8500.3,26000.0,3.2,45.6
...
```

关键：除了目标值外，还需要从最优解中提取以下决策结构信息：
- `route_cost`：总运输成本（所有时期路径成本之和）
- `inventory_cost`：总库存成本（收集点库存 + 工厂原料库存 + 成品库存）
- `production_cost`：总生产成本（固定设置成本 + 可变生产成本）
- `avg_visit_frequency`：平均每个收集点在规划期内被访问的次数 = Σ_i Σ_t z_{it} / n
- `avg_pickup_per_visit`：平均每次访问的取货量 = Σ_i Σ_t q_{it} / Σ_i Σ_t z_{it}

### 实验 2B：车辆数 K 的影响

在基准实例上，保持其他参数不变，变化 K：

| 方案 | K 值 |
|------|:----:|
| K=1  | 1    |
| K=2  | 2（基准） |
| K=3  | 3    |
| K=4  | 4    |

**注意**：K=1 可能导致某些实例不可行（单车容量不够），如果不可行则记录 `feasible=False`。

**输出格式**：`sensitivity_K.csv`，字段同上。

### 实验 2C：生产设置成本 f 的影响

在基准实例上，保持其他参数不变，将生产固定设置成本 f 按以下比例缩放：

| 方案 | f 倍率 |
|------|:------:|
| f-50%   | 0.50 |
| f-100%  | 1.00（基准） |
| f-200%  | 2.00 |
| f-500%  | 5.00 |

**实现方式**：在 Instance 加载后、求解前，直接修改 `f`（设置成本）的值。

**输出格式**：`sensitivity_f.csv`

```csv
instance,f_ratio,f_value,feasible,optimal,objective,time_sec,route_cost,inventory_cost,production_cost,num_production_periods,avg_production_batch,raw_inv_avg,finished_inv_avg
MVPRP1_10_6_2,0.50,125.0,True,True,43000.12,18.5,12000.5,7500.3,23500.0,4,75.2,30.5,20.1
...
```

额外字段：
- `num_production_periods`：有生产活动的时期数 = Σ_t y_t
- `avg_production_batch`：平均生产批量 = Σ_t p_t / Σ_t y_t
- `raw_inv_avg`：工厂原料库存平均水平 = Σ_t I^0_t / T
- `finished_inv_avg`：工厂成品库存平均水平 = Σ_t P^0_t / T

### 实验 2D：库存持有成本的影响

在基准实例上，保持其他参数不变，同比例缩放**所有三级库存持有成本**（收集点 h_i、工厂原料 h_0、成品 h_p）：

| 方案 | h 倍率 |
|------|:------:|
| h-50%  | 0.50 |
| h-100% | 1.00（基准） |
| h-200% | 2.00 |
| h-500% | 5.00 |

**输出格式**：`sensitivity_h.csv`

```csv
instance,h_ratio,feasible,optimal,objective,time_sec,route_cost,inventory_cost,production_cost,supplier_inv_avg,raw_inv_avg,finished_inv_avg,avg_visit_frequency
...
```

额外字段：
- `supplier_inv_avg`：收集点库存平均水平 = Σ_i Σ_t I_{it} / (n × T)

### 实现步骤

1. 阅读 Instance.java，找到 Q、K、f、h_i、h_0、h_p 这些参数的存储位置

2. 写一个参数扫描的 main 方法或脚本：
   - 加载基准实例
   - 修改目标参数
   - 调用 LBBD 求解
   - 从最优解中提取决策结构信息（成本分解、访问频率、库存水平等）
   - 写入 CSV

3. 如果从最优解中提取成本分解信息的接口不存在，需要在 LBBD 的解输出部分添加：遍历最优解的 z_{it}、q_{it}、路径方案、y_t、p_t、I^0_t、P^0_t，分别求和

4. 每个实验（2A/2B/2C/2D）× 2 个基准实例 × 4 个参数水平 = 8 次运行，四组共 32 次，时限 600 秒/次，总计最多 ~5.3 小时

---

## 输出文件汇总

实验完成后，应产出以下 5 个 CSV 文件：

```
ablation_results.csv      # 实验一：3 变体 × 16 实例 = 48 行
sensitivity_Q.csv         # 实验 2A：4 水平 × 2 实例 = 8 行
sensitivity_K.csv         # 实验 2B：4 水平 × 2 实例 = 8 行
sensitivity_f.csv         # 实验 2C：4 水平 × 2 实例 = 8 行
sensitivity_h.csv         # 实验 2D：4 水平 × 2 实例 = 8 行
```

请将这些文件放在项目根目录或指定的输出目录中。
