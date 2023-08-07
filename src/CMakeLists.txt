cmake_minimum_required(VERSION 3.15)
# 指定工程名字
project(test_fem)
# set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
# 指定C++标准，CMAKE_CXX_STANDARD是一个宏
# set(CMAKE_CXX_STANDARD 17)
# 指定要编译的源文件
file(GLOB SRC_PATH ./*.cpp)
# 指定可执行文件的输出路径
set(EXECUTABLE_OUTPUT_PATH ${EXEC_PATH})
# 指定头文件的搜索路径
include_directories(
${HEAD_PATH}
) 

add_executable(test_fem ${SRC_PATH})