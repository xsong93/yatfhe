#! /bin/bash
set -e

# HEXL must be built with AVX512 enabled -- a scalar-only HEXL costs ~4x on every
# NTT, which would silently invalidate any blind-rotate timing measured here.
# Point HEXL_ROOT at such an install prefix; the default is the sibling hexl
# checkout built with -DCMAKE_BUILD_TYPE=Release.
HEXL_ROOT="${HEXL_ROOT:-$HOME/Github/hexl/install-avx512}"

# Build type and -march=native come from CMakeLists.txt, which applies them when
# this project is configured standalone, as it is here.
CMAKE_ARGS=(-DINSTALL="OFF" -DTORUS_TYPE="32" -DPRINTER_ON="ON" -DENABLE_TIMER="ON")
if [ -d "$HEXL_ROOT" ]; then
    CMAKE_ARGS+=(-DCMAKE_PREFIX_PATH="$HEXL_ROOT")
else
    echo "warning: HEXL_ROOT=$HEXL_ROOT not found, falling back to system HEXL" >&2
    echo "         that build has no AVX512 -- timings will be ~4x pessimistic" >&2
fi

cmake -B CMAKE_BUILD "${CMAKE_ARGS[@]}"
cmake --build CMAKE_BUILD -j$(nproc)
