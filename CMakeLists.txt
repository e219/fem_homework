cmake_minimum_required(VERSION 3.15)
# 指定工程名字
project(Fem)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
# 指定C++标准，CMAKE_CXX_STANDARD是一个宏
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

set(EXEC_PATH ${PROJECT_SOURCE_DIR}/bin)

set(INTERN_HEAD ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(EXTERN_HEAD /usr/include/eigen3)

set(HEAD_PATH ${INTERN_HEAD} ${EXTERN_HEAD})

add_subdirectory(src)
add_subdirectory(test)