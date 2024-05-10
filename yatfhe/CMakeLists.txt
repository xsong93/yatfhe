project(yatfhe)

file(GLOB SOURCES
        ./*.cpp
        ./src/*.cpp
)

set(SOURCE_FILES ${SOURCES})
add_library(yatfhe_lib STATIC ${SOURCE_FILES})

target_include_directories(yatfhe_lib PRIVATE
        include/yautil
        include/yatfhe
)
target_include_directories(yatfhe_lib PUBLIC include)

# Build GMP
message("Linking GMP")
find_package(PkgConfig REQUIRED)
pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)
target_link_libraries(yatfhe_lib PkgConfig::GMPXX)
