message(STATUS "Building yatfhe...")
project(yatfhe)

set(HEXL_REQUIRED TRUE)

option(ENABLE_TIMER "Enable execution timers" OFF)
option(PRINTER_ON "Enable debug printing" OFF)
option(DEBUG_MODE "Enable debug mode" OFF)
option(TERNARY "Use ternary secret instead of binary" OFF)
option(USE_CRT "Use CRT" OFF)

set(TORUS_TYPE "32" CACHE STRING "Torus type (32, 56)")

file(GLOB SOURCES
        ./*.cpp
        ./src/*.cpp
)

set(SOURCE_FILES ${SOURCES})
add_library(yatfhe_lib STATIC ${SOURCE_FILES})

#target_include_directories(yatfhe_lib PRIVATE
#        include/yautil
#        include/yatfhe
#)
target_include_directories(yatfhe_lib PUBLIC include)

target_compile_definitions(yatfhe_lib PUBLIC
        TORUS${TORUS_TYPE}
        $<$<BOOL:${ENABLE_TIMER}>:ENABLE_TIMER>
        $<$<BOOL:${PRINTER_ON}>:PRINTER_ON>
        $<$<BOOL:${DEBUG_MODE}>:DEBUG_MODE>
        $<$<BOOL:${TERNARY}>:TERNARY>
        $<$<BOOL:${USE_CRT}>:USE_CRT>
)

# Build GMP
message(STATUS "Linking GMP...")
find_package(PkgConfig REQUIRED)
pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)
target_link_libraries(yatfhe_lib PkgConfig::GMPXX)
message(STATUS "Linking GMP...DONE!!")

# Link HEXL
message(STATUS "Linking HEXL...")
if(HEXL_REQUIRED)
    find_package(HEXL REQUIRED)
    message(STATUS "HEXL found, setting USING_NTT=HEXL")
    target_link_libraries(yatfhe_lib HEXL::hexl)
    add_compile_definitions(USING_NTT=HEXL)
else()
    message(STATUS "HEXL not required, setting USING_NTT=NORMAL")
    add_compile_definitions(USING_NTT=NORMAL)
endif()
message(STATUS "Linking HEXL...DONE!!")

# Link the threading library (pthread on Unix)
message(STATUS "Linking Threads...")
find_package(Threads REQUIRED)
target_link_libraries(yatfhe_lib Threads::Threads)
message(STATUS "Linking Threads...DONE!!")

# Add installation command
if(INSTALL)
    message(STATUS "Installing headers to /usr/local/include/yatfhe /usr/local/include/yautil")

    # Specify where to install the headers and library
    install(TARGETS yatfhe_lib
            ARCHIVE DESTINATION /usr/local/lib  # Static library installation path
            LIBRARY DESTINATION /usr/local/lib  # Shared library installation path
            RUNTIME DESTINATION /usr/local/bin  # Executable installation path
    )
    add_custom_command(
            TARGET yatfhe_lib PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E remove_directory /usr/local/include/yatfhe /usr/local/include/yautil
    )
    install(DIRECTORY include/yatfhe DESTINATION /usr/local/include)
    install(DIRECTORY include/yautil DESTINATION /usr/local/include)
endif ()

message(STATUS "Building yatfhe...DONE!!")