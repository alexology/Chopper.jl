#!/usr/bin/env bash
# Outputs of the original programs for the validation scripts, on Linux (or
# WSL), with the release binaries of chopper 0.14.1 and minimap2 2.30:
#
#     bash validation/make_oracles.sh chopper-linux minimap2 reads.fastq.gz contaminants.fa outdir
#
# Writes, from `reads.fastq.gz`: q20_t1.fq, split_t1.fq, best_t1.fq,
# gcinv_t1.fq (compare_outputs.jl); and, for the first 20,000 reads
# (outdir/r20k.fq): contam.fq (chopper -c) and contam.paf (minimap2), for
# compare_contaminants.jl. The output directory must have room for about
# 4 times the uncompressed input; do not use a small /tmp (tmpfs).
set -eu
chopper=$1 minimap2=$2 reads=$3 contam=$4 out=$5
mkdir -p "$out"
run() { name=$1; shift; "$chopper" -t 1 -i "$reads" "$@" > "$out/${name}_t1.fq"; }
run q20 -q 20
run split --trim-approach split-by-low-quality --cutoff 15 -l 200
run best -q 15 --trim-approach best-read-segment --cutoff 18 -l 500
run gcinv --mingc 0.35 --maxgc 0.45 --inverse -q 12 --maxlength 3000
gzip -cd "$reads" | head -n 80000 > "$out/r20k.fq"
"$chopper" -t 1 -c "$contam" -i "$out/r20k.fq" > "$out/contam.fq"
"$minimap2" -x lr:hq -c --no-hash-name -t 4 "$contam" "$out/r20k.fq" > "$out/contam.paf" 2> /dev/null
md5sum "$out"/*.fq
