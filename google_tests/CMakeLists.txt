project(google_tests)

add_subdirectory(lib)

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

#set(SOURCE_FILES ntt_test.cpp)

file(GLOB SOURCE_FILES *.cpp)

add_executable(google_tests_run ${SOURCE_FILES}) # adding the google_tests_run target

target_link_libraries(google_tests_run yatfhe_lib) # linking google_tests_run with yatfhe_lib which will be tested
target_link_libraries(google_tests_run gtest gtest_main)

# Build GMP
if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    message("Building on macOS")
    set(GMP_DIR /opt/homebrew/Cellar/gmp/6.3.0)
    include_directories(${GMP_DIR}/include)
    link_directories(${GMP_DIR}/lib)
    target_link_libraries(yatfhe_lib "${GMP_DIR}/lib/libgmp.a")
    target_link_libraries(yatfhe_lib "${GMP_DIR}/lib/libgmpxx.a")
else()
    message("Building on Linux")
    find_library(GMP_LIB gmp)
    target_link_libraries(yatfhe_lib libgmp)
endif ()