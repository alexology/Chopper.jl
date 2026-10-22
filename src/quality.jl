# Quality and composition statistics of a read.
#
# Qualities are handled as the raw bytes of the FASTQ quality line (Phred+33),
# so that no decoded copy of the line is needed.

"""Offset of the Phred+33 (Sanger, Illumina 1.8+, ONT, PacBio) quality encoding."""
const PHRED_OFFSET = 0x21

# Error probability 10^(-q/10) of every Phred score 0..255. A table lookup
# replaces a `^` call per base in the hot loop; covering the whole `UInt8`
# range means that any byte indexes safely. (`exp10` equals the C library's
# `pow(10, x)` used by chopper in all but 4 entries, in the last bit; `^`
# differs in 13.)
const PHRED_PROBABILITY = ntuple(i -> exp10((i - 1) / -10.0), 256)

"""
    phred_to_probability(q::Integer) -> Float64

Error probability `10^(-q/10)` of the Phred score `q` (0–255).
"""
@inline phred_to_probability(q::Integer) = @inbounds PHRED_PROBABILITY[(q % UInt8) + 1]

"""
    mean_quality(qual::AbstractVector{UInt8}) -> Float64

Mean quality of a read, from the raw Phred+33 bytes of its quality line,
computed as chopper (and NanoFilt) do: the Phred scores are converted to error
probabilities, averaged, and the mean probability is converted back to a Phred
score, `-10 log10(mean(p))`. This is lower than the arithmetic mean of the
scores whenever the scores vary.

The sum runs left to right, as in chopper. The result can differ from
chopper's in the last bit (Julia's `exp10` and `log10` against the C
library's), which matters only for a read whose mean quality lies exactly
on a threshold. An empty read gives `NaN`.
"""
function mean_quality(qual::AbstractVector{UInt8})
    total = 0.0
    @inbounds for q in qual
        total += phred_to_probability(q - PHRED_OFFSET)
    end
    return -10.0 * log10(total / length(qual))
end

"""
    gc_content(seq::AbstractVector{UInt8}) -> Float64

Fraction of `G`/`C` bases (either case) among all bases of `seq`, including
`N` and other symbols in the denominator, as in chopper. An empty sequence
gives `NaN`.
"""
function gc_content(seq::AbstractVector{UInt8})
    n = 0
    @inbounds for b in seq
        c = b | 0x20                     # lower case
        n += (c == UInt8('g')) | (c == UInt8('c'))
    end
    return n / length(seq)
end
