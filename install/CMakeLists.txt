cmake_minimum_required (VERSION 2.8) # 指定需要的最小的cmake版本.
project (DEOM) # 工程名(DEOM).
set (DEOM_VERSION_MAJOR 0) # 主版本号.
set (DEOM_VERSION_MINOR 0) # 次版本号.

set (CMAKE_VERBOSE_MAKEFILE ON) # 不要显示执行构建过程中详细的信息.

file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/bin) # 存放生成的可执行文件.
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/include) # 存放需要include的头文件.
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib) # 存放静态static链接(.a)文件.
#file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/share) # 存放动态shared链接(.so)文件(本project中暂不需要).

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin) # 可执行文件输出目录.
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib) # Module库文件输出目录.
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib) # Static库文件输出目录.

set (CMAKE_CXX_FLAGS "-std=c++11 -g -Wall -Ofast -march=ivybridge")

include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/src/deom
                    ${PROJECT_SOURCE_DIR}/src/thirdparty/trie
                   )

link_directories(${PROJECT_SOURCE_DIR}/lib)

add_subdirectory(${PROJECT_SOURCE_DIR}/src/deom)

FILE(GLOB apps_src ${PROJECT_SOURCE_DIR}/src/apps/*.cpp)  # GLOB选项将会为所有匹配查询表达式(*.cpp)的文件生成一个文件list, 并将该list存储进变量apps_src里.

foreach (p ${apps_src})
    get_filename_component(p_bin ${p} NAME_WE)
    add_executable(${p_bin} ${p})
    target_link_libraries(${p_bin} deom json11 armadillo lapack blas)
endforeach ()

include(ExternalProject)

if (EXISTS ${PROJECT_SOURCE_DIR}/lib/libjson11.a)
else()
  ExternalProject_Add(json11
    GIT_REPOSITORY https://github.com/dropbox/json11.git
    PREFIX ${PROJECT_BINARY_DIR}/thirdparty
    INSTALL_DIR ${PROJECT_SOURCE_DIR}
    INSTALL_COMMAND cp <SOURCE_DIR>/json11.hpp <INSTALL_DIR>/include && cp <BINARY_DIR>/libjson11.a <INSTALL_DIR>/lib
  )
endif()

#if (EXISTS ${PROJECT_SOURCE_DIR}/include/armadillo)
#else()
#  ExternalProject_Add(armadillo
#    URL http://sourceforge.net/projects/arma/files/armadillo-9.880.1.tar.xz
#    PREFIX ${PROJECT_BINARY_DIR}/thirdparty
#    INSTALL_DIR ${PROJECT_SOURCE_DIR}
#    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#  )
#endif()
