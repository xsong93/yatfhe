message(STATUS "Building tests...")
project(google_tests)

add_subdirectory(lib)

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

file(GLOB SOURCE_FILES *.cpp)

add_executable(google_tests_run ${SOURCE_FILES}) # adding the google_tests_run target

target_link_libraries(google_tests_run yatfhe_lib) # linking google_tests_run with yatfhe_lib which will be tested
target_link_libraries(google_tests_run gtest gtest_main)

message(STATUS "Building tests...DONE!!")