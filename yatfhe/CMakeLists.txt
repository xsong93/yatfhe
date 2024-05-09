project(yatfhe)

file(GLOB SOURCES
        ./*.cpp
        ./src/*.cpp
)

set(SOURCE_FILES ${SOURCES})
find_library(GMP_LIB gmp)
add_library(yatfhe_lib STATIC ${SOURCE_FILES})

target_include_directories(yatfhe_lib PRIVATE
        include/yautil
        include/yatfhe
)
target_include_directories(yatfhe_lib PUBLIC include)

target_link_libraries(yatfhe_lib gmp)