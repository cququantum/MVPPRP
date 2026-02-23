# MVPPRP

项目执行需要用到cplex，我安装的cplex版本是CPLEX_Studio2211
执行代码的指令是export MAVEN_OPTS="-Djava.library.path=/Applications/CPLEX_Studio2211/cplex/bin/arm64_osx"
mvn compile exec:java -Dexec.mainClass="Main"
