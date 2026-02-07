# MVPPRP

项目执行需要用到cplex，我安装的cplex版本是CPLEX_Studio2211
执行代码的指令是export MAVEN_OPTS="-Djava.library.path=/Applications/CPLEX_Studio2211/cplex/bin/arm64_osx"
mvn compile exec:java -Dexec.mainClass="Main"

## model.tex 修改记录（为保证原始模型与改写模型一致）

### 1) 修改位置：`model.tex` 第 2.4 节，`e_{ivt}` 定义（`v=0` 情况）

- 修改前：
```tex
h_i \sum_{j=0}^{t-1} \left(I_{i0} + \sum_{r=1}^{j} s_{ir}\right), \quad \text{如果 } v = 0
```
- 修改后：
```tex
h_i \sum_{j=1}^{t-1} \left(I_{i0} + \sum_{r=1}^{j} s_{ir}\right), \quad \text{如果 } v = 0
```
- 修改原因：原式会额外计入初始库存 `I_{i0}` 一期持有成本，导致改写模型目标值与原始模型相差常数项。

### 2) 修改位置：`model.tex` 第 2.4 节，改写目标函数中的生产成本项

- 修改前：
```tex
\underbrace{u p_t + f y_t)}_{\text{生产成本}}
```
- 修改后：
```tex
\underbrace{u p_t + f y_t}_{\text{生产成本}}
```
- 修改原因：修正括号笔误，避免公式书写不规范。

### 3) 修改位置：`model.tex` 第 2.4 节，车辆路径约束中的流量平衡约束索引域

- 修改前：
```tex
\sum_{j \in V} x_{ijt} = \sum_{j \in V} x_{jit}, \quad \forall i \in V, t \in T
```
- 修改后：
```tex
\sum_{j \in V} x_{ijt} = \sum_{j \in V} x_{jit}, \quad \forall i \in V_S, t \in T
```
- 修改原因：与原始模型一致，流量平衡约束仅对收集点集合 `V_S` 施加。
