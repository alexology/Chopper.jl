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
[Validation](#validation)), and it does not depend on the number of threads.

Julia ≥ 1.10 (tested with 1.12). Dependencies: FASTX.jl (FASTQ/FASTA parsing),
CodecZlibNG.jl (gzip, via zlib-ng), TranscodingStreams.jl, StringViews.jl,
PrecompileTools.jl. Optional: CodecBzip2, CodecXz, CodecZstd for bzip2, xz and
zstd input.

## Install from GitHub

This package is not registered in the Julia General registry. Install it
from GitHub together with its unregistered dependency:

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/Minimap2.jl")
Pkg.add(url="https://github.com/alexology/Chopper.jl")
```

## Installation

The package is not registered. From a local copy:

```julia
using Pkg
Pkg.develop(path = "path/to/Chopper.jl")
```

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

## Differences from chopper

Details in [docs/DIFFERENCES.md](docs/DIFFERENCES.md). In short:

- **deterministic output**: records are written in input order with any
  number of threads; chopper writes them in the order its threads finish
  whenever `-t > 1`;
- **library API**: strategies as types (invalid combinations of options
  cannot be expressed), 1-based `UnitRange` segments, statistics per filter
  returned as `ChopperStats`, errors as exceptions (`ModTagError`) instead
  of `process::exit`, no global state;
- **robustness**: a read with an empty sequence is discarded by the length
  filter (chopper panics and stops writing, exit status 101); input
  compression recognised from the content also on stdin, including zstd;
  output to a file, optionally gzip-compressed;
- `--inverse` stops at the first failed filter as normal mode does (the read
  is kept either way); chopper evaluates every filter, including the
  alignment against contaminants;
- **contaminants without the minimap2 binary**: Minimap2.jl maps the reads
  with the `lr:hq` preset as chopper's minimap2 crate does
  ([docs/CONTAMINANTS.md](docs/CONTAMINANTS.md)).

The mean quality uses Julia's `exp10` and `log10`, which can differ from the
C library used by chopper in the last bit: only a read whose mean quality is
exactly on a threshold could be classified differently (none in the
validation data).

## Validation

- `test/runtests.jl` (`julia -t 4 --project=. -e 'using Pkg; Pkg.test()'`,
  762 tests, about 4 minutes): chopper's own Rust unit tests, ported;
  byte-for-byte comparison with the output, stderr and exit status of the
  chopper 0.14.1 Linux binary on 35 cases (`test/reference/`, made by
  `tools/make_reference.sh`: every filter and trimming strategy, inverse,
  contaminants, warnings and errors, header edge cases, MM/ML/MN rewriting,
  broken tags), each with 1 and 4 threads; match lengths of the contaminant
  screen identical to minimap2 2.30's on chopper's test reads;
  compressed inputs (gzip, multi-member gzip, bzip2, xz, zstd, stdin).
- Real data (`validation/`), 200,000 ONT reads of the user's metabarcoding
  run (R10.4.1, dorado hac 5.2.0, 1.09 GB, tags in the header), against
  chopper 0.14.1 with one thread: output identical for `-q 20`;
  `split-by-low-quality --cutoff 15 -l 200` (799,589 records);
  `-q 15 best-read-segment --cutoff 18 -l 500`; `--mingc 0.35 --maxgc 0.45
  --inverse -q 12 --maxlength 3000`; with 1 and 4 threads.
- Contamination screen, against chopper and minimap2 2.30 (the version in
  chopper 0.14.1): identical classification and output on chopper's test data and on 20,000 of
  those reads screened against amplicon sequences of the same sample (16,070
  contaminants) and of another sample (none); see
  [docs/CONTAMINANTS.md](docs/CONTAMINANTS.md).

Speed on the 200,000 reads (Windows 11, 4 cores; chopper in WSL on the same
machine):

| | chopper 0.14.1 | Chopper.jl |
|---|---|---|
| `-q 20`, gzip input, 1 / 4 threads | 4.5 / 5.3 s | 5.4 / 4.3 s |
| `-q 20`, plain input, 4 threads | 1.6 s | 0.84 s |
| `-c`, 20,000 reads, 32 amplicons (16,070 contaminants) | 46 s (1 thread) | 24 s (4 threads) |
| `-c`, 20,000 reads, 3 amplicons of another species | 63 s (1 thread) | 24 s (4 threads) |

With gzip input both are limited by decompression (zlib-rs in chopper,
zlib-ng here, in a separate task). The contamination screen is limited by
Minimap2.jl, about 1.6 times slower than minimap2 per thread; reads are
mapped in parallel.

## Citation and license

Please cite chopper: De Coster W., Rademakers R. (2023) NanoPack2: population
scale evaluation of long-read sequencing data. *Bioinformatics* 39(5):
btad311. Chopper.jl is released under the MIT license ([LICENSE](LICENSE)). chopper
is also MIT-licensed; its copyright notice and minimap2's are in
[NOTICE](NOTICE). The upstream version is recorded in [UPSTREAM.toml](UPSTREAM.toml).
