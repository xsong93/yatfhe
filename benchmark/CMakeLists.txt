project(benchmark)

set(SOURCE_FILES bench.cpp)

add_executable(benchmark ${SOURCE_FILES})
add_executable(swap_copy swap_copy.cpp)

target_link_libraries(benchmark LINK_PRIVATE yatfhe_lib)
target_link_libraries(swap_copy LINK_PRIVATE yatfhe_lib)