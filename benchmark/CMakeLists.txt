project(benchmark)

add_executable(benchmark bench.cpp)
add_executable(swap_copy swap_copy.cpp)
add_executable(mod_mult mod_mult.cpp)
add_executable(crt crt.cpp)

target_link_libraries(benchmark PRIVATE yatfhe_lib)
target_link_libraries(swap_copy PRIVATE yatfhe_lib)
target_link_libraries(mod_mult PRIVATE yatfhe_lib)
target_link_libraries(crt PRIVATE yatfhe_lib)