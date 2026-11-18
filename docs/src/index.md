# Chopper.jl

Filtering and trimming of long-read FASTQ files (Oxford Nanopore, PacBio) in
native Julia: a port of [chopper](https://github.com/wdecoster/chopper) 0.14.1
by Wouter De Coster (Rust; the successor of NanoFilt and NanoLyse). No
external programs: the contamination screen, which chopper delegates to
minimap2, uses [Minimap2.jl](https://github.com/alexology/Minimap2.jl), a native port of the same
minimap2 2.30.

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.

Reads are filtered on mean quality, length, GC content and similarity to
contaminant sequences, then optionally trimmed with one of four strategies.
With the same options the output is **byte-identical to chopper's** (see
[Validation](@ref "Validation")), and it does not depend on the number of threads.

Julia ≥ 1.10 (tested with 1.12). Dependencies: FASTX.jl (FASTQ/FASTA parsing),
CodecZlibNG.jl (gzip, via zlib-ng), TranscodingStreams.jl, StringViews.jl,
PrecompileTools.jl. Optional: CodecBzip2, CodecXz, CodecZstd for bzip2, xz and
zstd input.

## Installation

This package is not registered in the Julia General registry, and neither
is its dependency Minimap2.jl. Install both from GitHub, Minimap2.jl first:

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/Minimap2.jl")
Pkg.add(url="https://github.com/alexology/Chopper.jl")
```

or, from local copies, `Pkg.develop(path = "path/to/Minimap2.jl")` and then
`Pkg.develop(path = "path/to/Chopper.jl")`. With Julia ≥ 1.11 `Project.toml`
also points to `../Minimap2.jl` (`[sources]`), so a checkout next to
Minimap2.jl works directly.

## Command line

`bin/chopper.jl` takes the options of chopper 0.14.1 (`--help` lists them),
plus `-o/--output` (a file, gzip-compressed if it ends in `.gz`) and
`--summary` (counts per filter on stderr). Start Julia with the threads to use:

```sh
julia -t 4 --project=path/to/Chopper.jl path/to/Chopper.jl/bin/chopper.jl -q 10 -l 500 -i reads.fastq.gz > filtered.fastq
julia -t 4 --project=... bin/chopper.jl -q 20 -i reads.fastq.gz -o filtered.fastq.gz --summary
gunzip -c reads.fastq.gz | julia --project=... bin/chopper.jl --trim-approach trim-by-quality --cutoff 15
julia --project=... bin/chopper.jl -c lambda.fasta --trim-approach split-by-low-quality --cutoff 15 --split-window 5 -l 50 -i reads.fastq
```

Messages, warnings, errors and exit codes are those of chopper ("Kept N reads
out of M reads" on stderr). Julia's start-up takes a few seconds; for many
files, call `chopper` from a Julia session or script instead.

## Julia API

```julia
using Chopper

stats = chopper("reads.fastq.gz", "filtered.fastq.gz"; min_quality = 10, min_length = 500)
stats.segments, stats.reads                  # chopper's "Kept N reads out of M reads"
stats                                        # counts per filter (ChopperStats)

chopper("reads.fastq", "trimmed.fastq"; trim = SplitByLowQuality(15; window = 5), min_length = 50)
chopper("reads.fastq", "clean.fastq"; contaminants = "lambda.fasta", trim = FixedCrop(20, 20),
        update_mods = true)

# The same decisions for one read (raw Phred+33 quality line):
f = ReadFilter(min_quality = 12, max_length = 30_000)
read_segments(seq, qual; filter = f, trim = TrimByQuality(10))   # -> Vector{UnitRange{Int}}
read_status(f, seq, qual)                                       # PASSED or the first failed filter
mean_quality(codeunits(qual)), gc_content(codeunits(seq))
trim_ranges(BestReadSegment(15), qual)

idx = ContaminantIndex("lambda.fasta")       # reuse it across files
is_contaminant(idx, seq), best_match_length(idx, seq)
rewrite_mod_tags(["MM:Z:C+m?,1,1;", "ML:B:C,200,10"], seq, 11:19)
```

| chopper option | Julia |
|---|---|
| `-q`, `--maxqual` | `min_quality`, `max_quality` |
| `-l`, `--maxlength` | `min_length`, `max_length` |
| `--mingc`, `--maxgc` | `min_gc`, `max_gc` |
| `-c` | `contaminants` (path or `ContaminantIndex`), `contaminant_fraction` (0.9) |
| `--trim-approach fixed-crop --headcrop H --tailcrop T` | `trim = FixedCrop(H, T)` |
| `--trim-approach trim-by-quality --cutoff C` | `trim = TrimByQuality(C)` |
| `--trim-approach best-read-segment --cutoff C` | `trim = BestReadSegment(C)` |
| `--trim-approach split-by-low-quality --cutoff C --split-window W` | `trim = SplitByLowQuality(C; window = W)` |
| `--update-mods`, `--inverse` | `update_mods`, `inverse` |
| `-t` | `threads` (at most `Threads.nthreads()`) |

A trimming strategy is a subtype of `TrimStrategy` implementing
`trim_ranges!(out, strategy, qual)`, so new strategies plug into `chopper`
without other changes (the tests define one).

Semantics, as in chopper: the filters judge the whole, untrimmed read, in the
order length, quality, GC, contaminants; trimming is applied afterwards to
the reads they select, and only the minimum length is checked again on each
segment (the mean quality of a segment is not recomputed). With `inverse`,
the reads failing a filter are written (and trimmed) instead.

## Citation and license

Please cite chopper: De Coster W., Rademakers R. (2023) NanoPack2: population
scale evaluation of long-read sequencing data. *Bioinformatics* 39(5):
btad311. Chopper.jl is released under the MIT license ([LICENSE](https://github.com/alexology/Chopper.jl/blob/main/LICENSE)). chopper
is also MIT-licensed; its copyright notice and minimap2's are in
[NOTICE](https://github.com/alexology/Chopper.jl/blob/main/NOTICE). The upstream version is recorded in [UPSTREAM.toml](https://github.com/alexology/Chopper.jl/blob/main/UPSTREAM.toml).
