# MVPPRP

项目执行需要用到cplex，我安装的cplex版本是CPLEX_Studio2211
执行代码的指令是export MAVEN_OPTS="-Djava.library.path=/Applications/CPLEX_Studio2211/cplex/bin/arm64_osx"
mvn compile exec:java -Dexec.mainClass="Main"

测试 1/2/3/4_10_3_2 的指令是
export MAVEN_OPTS="-Djava.library.path=/Applications/CPLEX_Studio2211/cplex/bin/arm64_osx"
mvn compile exec:java -Dexec.mainClass="Main" -Dexec.args="data/MVPRP/MVPRP1_10_3_2.txt data/MVPRP/MVPRP2_10_3_2.txt data/MVPRP/MVPRP3_10_3_2.txt data/MVPRP/MVPRP4_10_3_2.txt"
