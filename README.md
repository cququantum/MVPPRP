# MVPPRP

项目执行需要用到cplex，我安装的cplex版本是CPLEX_Studio2211
执行代码的指令是export MAVEN_OPTS="-Djava.library.path=/Applications/CPLEX_Studio2211/cplex/bin/arm64_osx"
mvn compile exec:java -Dexec.mainClass="Main"

按照下面这个测试
测试 1/2/3/4_10_3_2 的指令是
export MAVEN_OPTS="-Djava.library.path=/Applications/CPLEX_Studio2211/cplex/bin/arm64_osx"
mvn compile exec:java -Dexec.mainClass="Main" -Dexec.args="data/MVPRP/MVPRP1_10_3_2.txt data/MVPRP/MVPRP2_10_3_2.txt data/MVPRP/MVPRP3_10_3_2.txt data/MVPRP/MVPRP4_10_3_2.txt"

## 备注

当前 `CVRP` 子问题的求解已经基本对齐 `LRP` 中的 branch-and-price 思路，但 `period-CVRP` 里还没有移植 `subset-row cuts`。

这一点目前不是正确性问题。按照 `reformulationModel.tex`、`lbbd.tex`、`speed.tex` 的模型语义，当前实现仍然可以通过 `exact pricing fallback + branching` 保证 `CVRP` 子问题精确求解；`subset-row cuts` 只属于进一步强化 LP 下界、改善大实例性能的增强项。

如果后续在更大的实例上发现 `CVRP` 子问题重新成为主要瓶颈，优先考虑把 `LRP` 中的 `subset-row cuts` 作为性能增强移植进 `period-CVRP`。现阶段 README 中的 `1/2/3/4_10_3_2` 测试已通过，不需要为了正确性补这一层。
