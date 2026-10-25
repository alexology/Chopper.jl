# Contamination screen

chopper's `-c contaminants.fasta` removes reads that align to a contaminant
(typically the lambda phage DNA control strand of ONT kits) with at least 90%
of the read length as matching bases. chopper maps every read with minimap2
(Rust binding `minimap2` 0.1.31, minimap2 2.30) configured as
`Aligner::builder().preset(LrHq).with_cigar().with_index(fasta)`, and checks
the "match length" (PAF column 10) of every reported alignment.

Chopper.jl does the mapping with [Minimap2.jl](../../Minimap2.jl), a native
Julia port of minimap2 2.30, configured in the same way (`src/contaminants.jl`):

1. **Index** (`ContaminantIndex`): `MinimapIndex` of the contaminant
   sequences with the `lr:hq` minimizers (k = 19, w = 19); the mapping
   options are then updated from the index, as `with_index` does (occurrence
   threshold from the 2·10⁻⁴ quantile, clamped to 50–500).
2. **Mapping**: `map_query` with base-level alignment (`F_CIGAR`), which
   runs minimap2's seeding, chaining, chain selection and ksw2 alignment.
3. **Decision** (Chopper): a read is a contaminant if the largest match
   length over its alignments is at least `contaminant_fraction` (0.9) times
   its length.

The index is only read while mapping, so the batches of `chopper()` share it
across threads.

Until 2026-10-21 Chopper.jl had its own reimplementation of the same mapping
(minimizers, chaining, banded alignment by anti-diagonals). Its match
lengths differed from minimap2's for 0.4–61% of the aligned reads (co-optimal
paths, narrower bands, Z-drop points), although the classification agreed on
the validation data. It was replaced by Minimap2.jl so that there is one
mapper in the suite.

## Agreement with minimap2 and chopper

`test/test_contaminants.jl` compares the match lengths with minimap2 2.30
(`-x lr:hq -c`, largest per read; `tools/make_contaminant_reference.jl`) on
chopper's `test.fastq` (250 reads up to 77 kb, ~10% errors) against
`random_contam.fa` (one 75 kb read): identical for all 250 reads (182
aligned).

`validation/compare_contaminants.jl`, 20,000 ONT reads (ND2–COI amplicons,
R10.4.1), minimap2 2.30 release binary (`-x lr:hq -c --no-hash-name`) and
chopper 0.14.1 (`-t 1 -c`):

| Contaminants | Reads aligned (minimap2 / Julia) | Match length identical | Contaminants (chopper / Julia) | Output identical |
|---|---|---|---|---|
| 3 amplicon sequences of another sample | 15,761 / 15,761 | 20,000 / 20,000 | 0 / 0 | yes |
| 32 amplicon sequences of the same and the other sample | 19,060 / 19,060 | 20,000 / 20,000 | 16,070 / 16,070 | yes |

`--no-hash-name` matters: minimap2 breaks ties between equal chains with a
hash seeded by the read name, and chopper, like Chopper.jl, maps reads
without their names. Without the option one read (`95b66141…`, 2,833 bases,
fraction 0.88 against two nearly tied references, below the threshold either
way) has a different match length; it was first reported here as a
Minimap2.jl discrepancy (2026-10-22: it is not one).

Speed (20,000 reads, `chopper()` with 4 threads): 24 s with the 32
amplicons and 24 s with the 3 of another species; chopper took 46 s and
63 s with one thread.
