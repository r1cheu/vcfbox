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
link_directories(${HTSLIB_ROOT}/lib)

find_package(Threads REQUIRED)

add_subdirectory(ext/CLI11 EXCLUDE_FROM_ALL)
add_subdirectory(ext/eigen EXCLUDE_FROM_ALL)

add_library(barkeep INTERFACE)
target_include_directories(
  barkeep INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/ext/barkeep)

add_executable(
  vcfbox
  src/main.cpp
  src/app/combine_command.cpp
  src/app/cli.cpp
  src/app/parentage_command.cpp
  src/combine/combine_genotypes.cpp
  src/common/line_writer.cpp
  src/common/lines.cpp
  src/common/path.cpp
  src/common/progress.cpp
  src/common/strings.cpp
  src/parentage/bitmatrix.cpp
  src/parentage/build_matrix.cpp
  src/parentage/likelihood.cpp
  src/parentage/parentage.cpp
  src/parentage/pileup.cpp
  src/parentage/sites.cpp)

target_include_directories(
  vcfbox PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${HTSLIB_ROOT}/include)
target_link_libraries(vcfbox PRIVATE ${HTSLIB_ROOT}/lib/libhts.so
                                     Threads::Threads CLI11::CLI11
                                     Eigen3::Eigen barkeep)
