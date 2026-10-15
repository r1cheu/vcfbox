find_package(Catch2 3 CONFIG REQUIRED)

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "test_*.cpp")
list(APPEND TEST_SOURCES file_fixture.cpp vcf_fixture.cpp)

add_executable(vcfbox_tests ${TEST_SOURCES})

target_include_directories(vcfbox_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(vcfbox_tests PRIVATE vcfbox_core Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(vcfbox_tests)
