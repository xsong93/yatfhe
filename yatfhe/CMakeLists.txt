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