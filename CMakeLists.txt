cmake_minimum_required(VERSION 3.18)
project(
  f1_conbine
  VERSION 0.0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_BUILD_TYPE Release)
add_compile_options(-march=native)

set(HTSLIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.pixi/envs/default")
link_directories(${HTSLIB_ROOT}/lib)
list(APPEND CMAKE_PREFIX_PATH ${HTSLIB_ROOT})

find_package(Threads REQUIRED)

find_package(CLI11 CONFIG REQUIRED)

add_library(barkeep INTERFACE)
target_include_directories(
  barkeep INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/ext/barkeep)

# Business logic without CLI wiring; shared by the binary and the tests.
add_library(
  vcfbox_core STATIC
  src/combine/combine_genotypes.cpp
  src/parentage/counts.cpp
  src/parentage/counts_io.cpp
  src/parentage/panel.cpp
  src/parentage/parentage.cpp
  src/parentage/pileup.cpp
  src/parentage/report.cpp
  src/parentage/score.cpp
  src/parentage/sites.cpp
  src/common/line_writer.cpp
  src/common/lines.cpp
  src/common/path.cpp
  src/common/progress.cpp
  src/common/strings.cpp)

target_include_directories(
  vcfbox_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src ${HTSLIB_ROOT}/include)
target_link_libraries(vcfbox_core PUBLIC ${HTSLIB_ROOT}/lib/libhts.so
                                         Threads::Threads barkeep)

add_executable(
  vcfbox src/main.cpp src/app/combine_command.cpp
         src/app/parentage_command.cpp src/app/cli.cpp)

target_link_libraries(vcfbox PRIVATE vcfbox_core CLI11::CLI11)

enable_testing()
add_subdirectory(tests)
