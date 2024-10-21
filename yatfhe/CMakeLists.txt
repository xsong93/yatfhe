project(yatfhe)

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

# Build GMP
message("Linking GMP")
find_package(PkgConfig REQUIRED)
pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)
target_link_libraries(yatfhe_lib PkgConfig::GMPXX)

# Specify where to install the headers and library
install(TARGETS yatfhe_lib
        ARCHIVE DESTINATION /usr/local/lib  # Static library installation path
        LIBRARY DESTINATION /usr/local/lib  # Shared library installation path
        RUNTIME DESTINATION /usr/local/bin  # Executable installation path (if needed)
)

# Add installation command for headers
install(DIRECTORY include/yatfhe DESTINATION /usr/local/include)
install(DIRECTORY include/yautil DESTINATION /usr/local/include)