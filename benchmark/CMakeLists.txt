project(benchmark)

set(SOURCE_FILES bench.cpp)

add_executable(benchmark ${SOURCE_FILES})

target_link_libraries(benchmark LINK_PRIVATE yatfhe_lib)