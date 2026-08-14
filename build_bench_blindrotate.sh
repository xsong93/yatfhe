#! /bin/bash
set -e

# 32, 56
TORUS="${TORUS:-32}"

# ISA baseline for our own code: native (default), avx512, avx2 or generic.
TARGET_ARCH="${TARGET_ARCH:-native}"

# Install prefixes for different HEXL build.
# Set HEXL_ROOT to override.
case "$TARGET_ARCH" in
    native|avx512) DEFAULT_HEXL_ROOT="/usr/local" ;;
    avx2)          DEFAULT_HEXL_ROOT="$HOME/.hexl/avx2" ;;
    *)             DEFAULT_HEXL_ROOT="$HOME/.hexl/generic" ;;
esac
HEXL_ROOT="${HEXL_ROOT-$DEFAULT_HEXL_ROOT}"

if [ ! -d "$HEXL_ROOT" ]; then
    echo "error: HEXL_ROOT=$HEXL_ROOT not found. Run ./hexl_build_all.sh" >&2
    exit 1
fi

BUILD_DIR="CMAKE_BUILD"

# You can disable debug messages by setting -DPRINTER_ON to OFF.
# Likewise, disable timing by setting -DENABLE_TIMER to OFF.
cmake -B "$BUILD_DIR" \
      -UHEXL_DIR \
      -DINSTALL="OFF" \
      -DTORUS_TYPE="$TORUS" \
      -DPRINTER_ON="ON" \
      -DENABLE_TIMER="ON" \
      -DTARGET_ARCH="$TARGET_ARCH" \
      -DCMAKE_PREFIX_PATH="$HEXL_ROOT"

cmake --build "$BUILD_DIR" -j$(nproc)
