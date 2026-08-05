#! /bin/bash
set -e

# ISA baseline for our own code: native (default), avx512, avx2 or generic.
# Use avx2/generic to produce binaries that run on machines without AVX-512.
TARGET_ARCH="${TARGET_ARCH:-native}"

# The default install follows TARGET_ARCH:
# AVX-512-tuned HEXL for native/avx512, the system HEXL for avx2/generic.
# Set HEXL_ROOT explicitly to override, or HEXL_ROOT= to force the system one.
case "$TARGET_ARCH" in
    native|avx512) DEFAULT_HEXL_ROOT="$HOME/Github/hexl/install-avx512" ;;
    *)             DEFAULT_HEXL_ROOT="" ;;
esac
HEXL_ROOT="${HEXL_ROOT-$DEFAULT_HEXL_ROOT}"

# Build type and ISA flags come from CMakeLists.txt, which applies them when
# this project is configured standalone, as it is here.
CMAKE_ARGS=(-DINSTALL="OFF" -DTORUS_TYPE="32" -DPRINTER_ON="ON" -DENABLE_TIMER="ON"
            -DTARGET_ARCH="$TARGET_ARCH")
if [ -n "$HEXL_ROOT" ] && [ -d "$HEXL_ROOT" ]; then
    CMAKE_ARGS+=(-DCMAKE_PREFIX_PATH="$HEXL_ROOT")
elif [ -n "$HEXL_ROOT" ]; then
    echo "warning: HEXL_ROOT=$HEXL_ROOT not found, falling back to system HEXL" >&2
    echo "         that build has no AVX512 -- timings will be ~4x pessimistic" >&2
fi

cmake -B CMAKE_BUILD "${CMAKE_ARGS[@]}"
cmake --build CMAKE_BUILD -j$(nproc)
