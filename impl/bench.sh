#!/bin/sh

echo "-------------------------------------"
echo "$1: $2"
echo "Running the benchmarks..."
# Reduced criterion budget by default: 20 samples of ~1s with 0.5s warmup
# keeps a full sweep down to minutes. Set FLATTEN_BENCH_QUALITY=high for the
# full 100-sample, 3s warmup + 5s measurement protocol. --bench bench keeps
# the arguments away from the lib unittest target, whose libtest harness
# rejects criterion's options.
if [ "$FLATTEN_BENCH_QUALITY" = "high" ]; then
    FLATTEN_INPUT="$2" cargo bench --bench bench $1
else
    FLATTEN_INPUT="$2" cargo bench --bench bench $1 -- --sample-size 20 --measurement-time 1 --warm-up-time 1
fi
echo "Extracting benchmark results..."
flatten-helper criterion $1 -i . -o ../notes/results/bench-$1-$2-$3.md
echo "Making graphs..."
flatten-helper graph -i ../notes/results/bench-$1-$2-$3.md -o ../notes/results/bench-$1-$2-$3.svg -t "$2 dataset: $1"
echo "$2 done."
