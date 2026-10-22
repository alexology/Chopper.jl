# Validation against chopper and minimap2

The unit and reference tests (`test/`) use chopper's small test files. These
scripts compare Chopper.jl with the original programs on real data, which is
not distributed with the package.

1. On Linux (or WSL), with the release binaries of chopper 0.14.1
   (`chopper-linux`) and minimap2 2.30 (`minimap2-2.30_x64-linux`):

   ```sh
   bash validation/make_oracles.sh chopper-linux minimap2 reads.fastq.gz contaminants.fa oracles/
   ```

2. Then, with Julia:

   ```sh
   julia -t 4 --project=. validation/compare_outputs.jl reads.fastq.gz oracles/
   julia -t 4 --project=. validation/compare_contaminants.jl oracles/r20k.fq contaminants.fa oracles/contam.paf oracles/contam.fq
   ```

`compare_outputs.jl` runs four settings (quality filter; split by low
quality; best read segment with a quality filter; GC window with `--inverse`)
with 1 thread and with all threads, and reports times and whether the output
is byte-identical. `compare_contaminants.jl` compares match lengths with
minimap2, the classification at 90%, and the output with chopper's.

Results on 2026-10-21 (200,000 R10.4.1 reads, dorado hac 5.2.0, of an ND2–COI
metabarcoding run; contaminants: savont amplicon sequences of the same and of
another sample): every output identical, see the README and
[docs/CONTAMINANTS.md](../docs/CONTAMINANTS.md).
