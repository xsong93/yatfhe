message(STATUS "Building benches...")

file(GLOB YA_BENCHMARK_SOURCES *.cpp)

foreach(bench_source ${YA_BENCHMARK_SOURCES})
    get_filename_component(target_name ${bench_source} NAME_WE)
    add_executable(${target_name} ${bench_source})
    target_include_directories(${target_name} PRIVATE ${CMAKE_SOURCE_DIR}/yatfhe/include)
    target_link_libraries(${target_name} PRIVATE yatfhe_lib)

endforeach()

message(STATUS "Building benches...DONE!!")
