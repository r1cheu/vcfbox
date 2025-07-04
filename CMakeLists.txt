cmake_minimum_required(VERSION 3.18)
project(
  f1_conbine
  VERSION 0.0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_BUILD_TYPE Release)

set(HTSLIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.pixi/envs/default")
include_directories(${HTSLIB_ROOT}/include)
link_directories(${HTSLIB_ROOT}/lib)

find_package(Threads REQUIRED)
add_executable(vcf_combine src/main.cpp src/vcf.cpp src/utils.cpp)
add_executable(test src/tester.cpp)
target_link_libraries(vcf_combine PRIVATE ${HTSLIB_ROOT}/lib/libhts.so
                                          Threads::Threads)
target_link_libraries(test PRIVATE ${HTSLIB_ROOT}/lib/libhts.so
                                   Threads::Threads)
