file (GLOB cur *.cpp) # GLOB选项将会为所有匹配查询表达式(*.cpp)的文件生成一个文件list, 并将该list存储进变量cur里.
add_library (deom ${cur}) # 将源文件(cur)生成static链接文件, 然后添加到工程中去.
add_dependencies(deom json11) # 令libdeom.a依赖于libjson11.a.
