message(STATUS "Building yatfhe...")
project(yatfhe)

file(GLOB SOURCES
        ./*.cpp
        ./src/*.cpp
)

set(SOURCE_FILES ${SOURCES}
        include/yautil/multi_threading.h
        src/multi_threading.cpp)
add_library(yatfhe_lib STATIC ${SOURCE_FILES})

#target_include_directories(yatfhe_lib PRIVATE
#        include/yautil
#        include/yatfhe
#)
target_include_directories(yatfhe_lib PUBLIC include)

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
    message(STATUS "HEXL found, setting NTT=HEXL")
    target_link_libraries(yatfhe_lib HEXL::hexl)
    add_compile_definitions(NTT=HEXL)
else()
    message(STATUS "HEXL not required, setting NTT=NORMAL")
    add_compile_definitions(NTT=NORMAL)
endif()
message(STATUS "Linking HEXL...DONE!!")

# Link the threading library (pthread on Unix)
message(STATUS "Linking Threads...")
find_package(Threads REQUIRED)
target_link_libraries(yatfhe_lib Threads::Threads)
message(STATUS "Linking Threads...DONE!!")

# Specify where to install the headers and library
install(TARGETS yatfhe_lib
        ARCHIVE DESTINATION /usr/local/lib  # Static library installation path
        LIBRARY DESTINATION /usr/local/lib  # Shared library installation path
        RUNTIME DESTINATION /usr/local/bin  # Executable installation path (if needed)
)

# Add installation command for headers
if(INSTALL_HEADERS)
    message(STATUS "Installing headers to /usr/local/include/yatfhe /usr/local/include/yautil")
    add_custom_command(
            TARGET yatfhe_lib PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E remove_directory /usr/local/include/yatfhe /usr/local/include/yautil
    )
    install(DIRECTORY include/yatfhe DESTINATION /usr/local/include)
    install(DIRECTORY include/yautil DESTINATION /usr/local/include)
endif ()

message(STATUS "Building yatfhe...DONE!!")