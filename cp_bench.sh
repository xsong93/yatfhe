#! /bin/bash

set -eu
for b in CMAKE_BUILD/ya_benchmark/*; do
    if [ -f "$b" ] && [ -x "$b" ]; then
        cp "$b" benchmark_run/server/.
        echo "copied $(basename "$b")"
    fi
done
