message(STATUS "Building google benches...")
project(google_benchmark)

find_package(benchmark REQUIRED)

file(GLOB BENCHMARK_SOURCES *.cpp)

get_filename_component(PROJECT_ROOT "../" ABSOLUTE)

foreach(bench_source ${BENCHMARK_SOURCES})
    get_filename_component(target_name ${bench_source} NAME_WE)

    message(STATUS "Adding benchmark target: ${target_name}")

    add_executable(${target_name} ${bench_source})

    target_include_directories(${target_name} PRIVATE
            ${CMAKE_SOURCE_DIR}/yatfhe/include
    )

    target_link_libraries(${target_name} PRIVATE
            benchmark::benchmark
            yatfhe_lib
            pthread
    )

endforeach()

message(STATUS "Building google benches...DONE!!")
