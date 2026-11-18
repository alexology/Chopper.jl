# Chopper.jl and chopper 0.14.1

What is the same, what differs on purpose, and what was found in chopper while
porting it.

## Same behaviour

With the same options and one thread in chopper, the output (records, bytes),
the stderr messages and the exit status are identical; `test/reference/`
checks this on 35 cases and `validation/` on real data. In particular:

- **Mean quality**: error probabilities from a table of `10^(-q/10)`, summed
  left to right, `-10 log10(mean)`, with Julia's `exp10` and `log10` (see
  below for the last bit).
- **Filter order and activation**: length, quality (only if `-q > 0` or
  `--maxqual < 1000`), GC (only if `--mingc` or `--maxgc` is set),
  contaminants. GC counts `G`/`C` in either case over all bases, `N`
  included in the denominator.
- **Trimming**: the four strategies with chopper's tie rules (the best read
  segment prefers, among equal scores, the longest and then the first; a
  running sum of exactly zero is not reset); filters judge the untrimmed read
  and only `-l` is checked again on each segment, not with
  `--trim-approach` unset.
- **Headers**: the trailing whitespace of the header line is removed, as
  chopper's FASTQ parser does; split reads get `_segment_<n>` on the read
  name, before any description or tab-separated tag.
- **Base modification tags** (`--update-mods`): MM/ML/MN rewritten with the
  same rules (htslib semantics, `N` counts every base, lower-case tag names,
  two-code lists, MM without ML, `+` before a skip count accepted like Rust's
  integer parser) and the same error messages; a broken tag is fatal and the
  output incomplete.
- **Command line**: same options, defaults, warnings for ignored trimming
  options, errors for missing ones.

## Differences on purpose

| | chopper | Chopper.jl |
|---|---|---|
| Output order with `-t > 1` | order in which threads finish (varies between runs) | input order, for any number of threads |
| Use from code | command line only | `chopper()`, `read_segments`, `read_status`, `trim_ranges`, `ContaminantIndex`, ... |
| Statistics | "Kept N reads out of M reads" | `ChopperStats`: reads/bases read and written, count per failed filter, reads emptied by trimming; `--summary` |
| Trimming options | flags validated at run time (`--cutoff` with `fixed-crop` is a warning) | one type per strategy; invalid combinations cannot be built in the API (the CLI keeps chopper's warnings) |
| Segments | 0-based `(start, end)` | 1-based `UnitRange` |
| Errors | `eprintln!` + `process::exit` from any thread | exceptions (`ModTagError` with the read name); the CLI turns them into chopper's messages and exit codes |
| Global state | counter of reads with MM tags | none (`stats.reads_with_mods`) |
| `--inverse` | evaluates every filter, including the alignment against contaminants | stops at the first failed filter: the read is kept anyway, so the result is the same |
| Empty reads (`@r` + empty sequence) | the FASTQ parser fails; chopper panics (exit 101) and stops writing | read normally; discarded by `-l` (≥ 1 by default) |
| Input compression | gzip, bzip2, xz for files (magic bytes or extension); stdin uncompressed | gzip (zlib-ng), bzip2, xz, zstd, from magic bytes, files and stdin |
| Output | stdout | stdout, or `-o` file (gzip if `.gz`) |
| Contaminant screen | minimap2 2.30 (C library, Rust crate) | Minimap2.jl, a Julia port of minimap2 2.30 ([CONTAMINANTS.md](CONTAMINANTS.md)) |
| Fraction of matches for contaminants | fixed 0.9 | `contaminant_fraction` keyword (0.9) |
| `-t` default | 4 | the Julia threads (`julia -t N`) |

Not supported (as in chopper): sequencing summary files. FASTQ records must
fit on four lines (FASTX.jl); chopper's parser also accepts sequences split
over several lines, which no current basecaller writes.

## Found in chopper

- **Empty reads crash the run.** A record with an empty sequence and quality
  makes the `bio` parser return `IncompleteRecord`; chopper panics ("problem
  parsing fastq record"), exits with status 101 and the records after it are
  lost (checked with the 0.14.1 binary).
- **Non-deterministic output** with more than one thread (documented
  upstream): on 200,000 reads `-t 1` and `-t 4` give different files with the
  same sorted lines.
- **The comment on `--inverse`** in `get_valid_segments` says that all filters
  must be evaluated in inverse mode; since a read is kept as soon as one
  filter fails, the early exit is also valid there.
- **Floating point and the platform**: chopper's mean quality depends on the
  C library (`powf` for the table, `log10`), so chopper itself may classify a
  read exactly on a threshold differently on Linux and macOS. Against glibc,
  Julia's `exp10` differs in 4 of the 256 table entries and `log10` in about
  2.5% of the arguments, always in the last bit; a port of glibc's functions
  was tried and dropped as unnecessary (no read of the validation data is
  affected).
- Not reproduced: chopper's comment on `par_bridge`/`for_each_init` (more
  initialisations than channels would drop reads with "failed to send the
  read for writing"); 1–32 threads kept the same 199,961 of 200,000 reads.
