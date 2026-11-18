# Differences and validation

## Differences from chopper

Details in [docs/src/DIFFERENCES.md](DIFFERENCES.md). In short:

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
  ([docs/src/CONTAMINANTS.md](CONTAMINANTS.md)).

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
  [docs/src/CONTAMINANTS.md](CONTAMINANTS.md).

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
