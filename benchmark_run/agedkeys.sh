#!/bin/bash
# Point the ablation's five key names at the most aged file of each key type.
#
# Why this exists.
# A key file written in the same session reads faster than an aged one.
# Only aged keys simulate the serving condition.
#
# These have a pool counterpart of the same key type and size:
#   A  BootstrappingKeyMP            <- BSK_GINX_*.bin
#   B  BootstrappingKeyWWL24         <- BSK_WWL+24_*.bin
#   E  BootstrappingKeyMPLazyPipeAlt <- BSK_LAZY_*.bin
# The other two key types exist only inside the ablation,
# which have aged for the length of the suite by the time the measurement runs.

set -eu
B="$(cd "$(dirname "$0")" && pwd)"
SRV="$B/server"
DIR="$SRV/agedkeys"
TENANT="${AGED_TENANT:-1}"      # which pool tenant is the aged reference
mkdir -p "$DIR"
rm -f "$DIR"/*.bin

link() {
    if [ ! -f "$SRV/$2" ]; then
        echo "agedkeys: $2 is missing; run stage 0 (gen_benchkeys) first" >&2
        exit 1
    fi
    ln -s "../$2" "$DIR/$1"
}

link COMP_A_GINX.bin        "BSK_GINX_${TENANT}.bin"
link COMP_C_WWL24.bin       "BSK_WWL+24_${TENANT}.bin"
link COMP_G_PIPE_ALT.bin    "BSK_LAZY_${TENANT}.bin"
link COMP_D_WWL24_ALT.bin   "COMP_D_WWL24_ALT.bin"
link COMP_E_PAR_LAZY.bin    "COMP_E_PAR_LAZY.bin"

# The mapping is only sound if the type really matches
# a silent mismatch would be timed as a valid row
check() {
    local got
    got=$(stat -Lc%s "$DIR/$1")
    if [ "$got" != "$2" ]; then
        echo "agedkeys: $1 is $got bytes, expected $2 ($3): wrong key type for that row" >&2
        exit 1
    fi
}
check COMP_A_GINX.bin      178372849 "GINX key, analogous to BSK_GINX"
check COMP_C_WWL24.bin     89446336  "succinct WWL+24 key, analogous to BSK_WWL+24"
check COMP_D_WWL24_ALT.bin 44890016  "plaintext-domain WWL+24 key, ablation only"
check COMP_E_PAR_LAZY.bin  89330040  "parallel-expansion key, ablation only"
check COMP_G_PIPE_ALT.bin  44890016  "pipelined lazy key, analogous to BSK_LAZY"

echo "agedkeys: server/agedkeys ready (3 pool symlinks + 2 stage-2 files, tenant $TENANT)"
for f in "$DIR"/*.bin; do
    printf "  %-24s -> %s\n" "$(basename "$f")" "$(readlink "$f")"
done
