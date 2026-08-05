cmake_minimum_required(VERSION 3.10)
message(STATUS "Building base...")
project(base)

set(CMAKE_CXX_STANDARD 17)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()
    message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

    # --- ISA baseline ---------------------------------------------------------
    # Pick an explicit level for reproducible or portable builds.
    #
    #   native  - tune for the build host (default; fastest, least portable)
    #   avx512  - require AVX-512F/DQ/VL/BW; runs only on AVX-512 hardware
    #   avx2    - AVX2 + FMA; runs on anything from Haswell onwards
    #   generic - no ISA flags at all
    #
    # This governs our own code only -- the HEXL you link has to match separately.
    # Point CMAKE_PREFIX_PATH at a matching HEXL install.
    set(TARGET_ARCH "native" CACHE STRING "ISA baseline: native, avx512, avx2, generic")
    set_property(CACHE TARGET_ARCH PROPERTY STRINGS native avx512 avx2 generic)

    option(NATIVE_ARCH "Deprecated; use TARGET_ARCH=generic" ON)
    if(NOT NATIVE_ARCH AND TARGET_ARCH STREQUAL "native")
        message(WARNING "NATIVE_ARCH=OFF is deprecated; use -DTARGET_ARCH=generic")
        set(TARGET_ARCH "generic")
    endif()

    include(CheckCXXCompilerFlag)

    # Add the flags only if the compiler actually accepts them.
    function(add_arch_flags label)
        string(REPLACE ";" " " flags_str "${ARGN}")
        string(MAKE_C_IDENTIFIER "COMPILER_SUPPORTS_ARCH_${label}" cache_var)
        check_cxx_compiler_flag("${flags_str}" ${cache_var})
        if(${cache_var})
            add_compile_options(${ARGN})
        else()
            message(WARNING "TARGET_ARCH=${label}: compiler rejects '${flags_str}'; "
                            "building with no ISA flags")
        endif()
    endfunction()

    if(TARGET_ARCH STREQUAL "native")
        add_arch_flags(native -march=native)
    elseif(TARGET_ARCH STREQUAL "avx512")
        add_arch_flags(avx512 -mavx512f -mavx512dq -mavx512vl -mavx512bw -mfma)
    elseif(TARGET_ARCH STREQUAL "avx2")
        add_arch_flags(avx2 -mavx2 -mfma)
    elseif(NOT TARGET_ARCH STREQUAL "generic")
        message(FATAL_ERROR
                "TARGET_ARCH must be one of: native, avx512, avx2, generic (got '${TARGET_ARCH}')")
    endif()

    if(TARGET_ARCH STREQUAL "avx512" AND EXISTS /proc/cpuinfo)
        file(READ /proc/cpuinfo CPUINFO)
        if(NOT CPUINFO MATCHES "avx512f")
            message(WARNING "TARGET_ARCH=avx512 but this host reports no AVX-512; "
                            "the binaries will not run here.")
        endif()
    endif()
    message(STATUS "TARGET_ARCH: ${TARGET_ARCH}")

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
