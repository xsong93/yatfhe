cmake_minimum_required(VERSION 3.10)
message(STATUS "Building base...")
project(base)

set(CMAKE_CXX_STANDARD 17)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()
    message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

    # Tune for the build host (AVX512). Turn OFF for portable binaries.
    option(NATIVE_ARCH "Compile with -march=native" ON)
    if(NATIVE_ARCH)
        include(CheckCXXCompilerFlag)
        check_cxx_compiler_flag(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
        if(COMPILER_SUPPORTS_MARCH_NATIVE)
            add_compile_options(-march=native)
        endif()
    endif()

    # OFF by default: measured a wash on this tree, since the hot leaf functions
    # (longModP, subTorus, signedGadgetDecomposition) are already header-defined.
    option(ENABLE_LTO "Build with link-time optimization" OFF)
    if(ENABLE_LTO)
        include(CheckIPOSupported)
        check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)
        if(IPO_SUPPORTED)
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
            message(STATUS "LTO: enabled")
        else()
            message(STATUS "LTO: unsupported (${IPO_ERROR})")
        endif()
    endif()
endif()

add_subdirectory(yatfhe)
add_subdirectory(ya_benchmark)
add_subdirectory(google_tests)
