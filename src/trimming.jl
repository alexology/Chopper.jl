# Trimming strategies.
#
# A strategy never copies or modifies a read: it returns the ranges of the
# read to keep (1-based, inclusive `UnitRange`s); more than one range means
# that the read was split. Only the writer slices the record.

"""
    TrimStrategy

Abstract type of the trimming strategies. A strategy `s` implements
`trim_ranges!(out, s, qual)`, which empties `out`, fills it with the ranges
(`UnitRange{Int}`, 1-based) of the read to keep, in increasing order, and
returns it. `qual` holds the raw Phred+33 bytes of the quality line.
"""
abstract type TrimStrategy end

"""
    trim_ranges(s::TrimStrategy, qual::AbstractVector{UInt8}) -> Vector{UnitRange{Int}}

Ranges of the read to keep under strategy `s`, from the raw Phred+33 bytes of
its quality line (`codeunits` of the line, or `Vector{UInt8}`). An empty
vector means that nothing is kept. See [`trim_ranges!`](@ref) for a version
that reuses the output vector.
"""
trim_ranges(s::TrimStrategy, qual::AbstractVector{UInt8}) =
    trim_ranges!(UnitRange{Int}[], s, qual)
trim_ranges(s::TrimStrategy, qual::AbstractString) = trim_ranges(s, codeunits(qual))

"""
    trim_ranges!(out::Vector{UnitRange{Int}}, s::TrimStrategy, qual) -> out

In-place version of [`trim_ranges`](@ref).
"""
function trim_ranges! end

"""
    NoTrim()

Keep the whole read. With this strategy (the default) the minimum length is
not checked a second time after trimming, as in chopper without
`--trim-approach`.
"""
struct NoTrim <: TrimStrategy end

function trim_ranges!(out::Vector{UnitRange{Int}}, ::NoTrim, qual::AbstractVector{UInt8})
    empty!(out)
    push!(out, 1:length(qual))
    return out
end

function check_cutoff(cutoff::Integer)
    0 <= cutoff <= 255 || throw(ArgumentError("cutoff must be a Phred score in 0:255, got $cutoff"))
    return UInt8(cutoff)
end

"""
    FixedCrop(head, tail)

Remove `head` bases from the start and `tail` bases from the end of every
read (`--trim-approach fixed-crop --headcrop head --tailcrop tail`). Reads
not longer than `head + tail` are discarded. At least one of the two must be
positive.
"""
struct FixedCrop <: TrimStrategy
    head::Int
    tail::Int
    function FixedCrop(head::Integer, tail::Integer)
        (head >= 0 && tail >= 0) || throw(ArgumentError("head and tail crop must be non-negative"))
        (head > 0 || tail > 0) ||
            throw(ArgumentError("at least one of head and tail crop must be greater than 0"))
        return new(head, tail)
    end
end
FixedCrop(; head::Integer = 0, tail::Integer = 0) = FixedCrop(head, tail)

function trim_ranges!(out::Vector{UnitRange{Int}}, s::FixedCrop, qual::AbstractVector{UInt8})
    empty!(out)
    n = length(qual)
    max(n - s.tail, 0) <= s.head && return out
    push!(out, (s.head + 1):(n - s.tail))
    return out
end

"""
    TrimByQuality(cutoff)

Trim bases with a Phred score below `cutoff` from both ends of the read,
stopping at the first base with score `>= cutoff` on each side
(`--trim-approach trim-by-quality --cutoff cutoff`). Low-quality bases inside
the read are kept. A read without any base `>= cutoff` is discarded.
"""
struct TrimByQuality <: TrimStrategy
    cutoff::UInt8
    TrimByQuality(cutoff::Integer) = new(check_cutoff(cutoff))
end

function trim_ranges!(out::Vector{UnitRange{Int}}, s::TrimByQuality, qual::AbstractVector{UInt8})
    empty!(out)
    n = length(qual)
    first = 1
    @inbounds while first <= n && qual[first] - PHRED_OFFSET < s.cutoff
        first += 1
    end
    last = n
    @inbounds while last > first && qual[last] - PHRED_OFFSET < s.cutoff
        last -= 1
    end
    last >= first && push!(out, first:last)
    return out
end

"""
    BestReadSegment(cutoff)

Keep the single segment of the read with the highest score, where every base
scores `p(cutoff) - p(q)` and `p` is the error probability of a Phred score
(`--trim-approach best-read-segment --cutoff cutoff`). This is the modified
Mott algorithm of Phred quality trimming: bases better than `cutoff`
contribute positively, worse bases negatively. Among segments with the same
score the longest (then the first) is kept. A read whose bases are all worse
than `cutoff` is discarded.
"""
struct BestReadSegment <: TrimStrategy
    cutoff::UInt8
    probability::Float64
    BestReadSegment(cutoff::Integer) =
        new(check_cutoff(cutoff), phred_to_probability(check_cutoff(cutoff)))
end

function trim_ranges!(out::Vector{UnitRange{Int}}, s::BestReadSegment, qual::AbstractVector{UInt8})
    empty!(out)
    # Kadane's maximum subarray, with the tie rules of chopper: an equal score
    # replaces the best segment only if longer, and a running sum of exactly
    # zero is not reset.
    best_start = 0
    best_stop = 0
    best_score = 0.0
    best_length = 0
    start = 1
    score = -1.0
    @inbounds for i in eachindex(qual)
        gain = s.probability - phred_to_probability(qual[i] - PHRED_OFFSET)
        if score < 0.0
            score = 0.0
            start = i
        end
        score += gain
        len = i - start + 1
        if best_score < score || (best_score == score && best_length < len)
            best_start = start
            best_stop = i
            best_score = score
            best_length = len
        end
    end
    best_start > 0 && push!(out, best_start:best_stop)
    return out
end

"""
    SplitByLowQuality(cutoff; window = 1)

Split reads at runs of at least `window` consecutive bases with a Phred score
below `cutoff`, and keep the high-quality parts between them
(`--trim-approach split-by-low-quality --cutoff cutoff --split-window
window`). Shorter low-quality runs stay inside their segment; low-quality
bases at the ends of a segment are removed. Each part must still reach the
minimum read length to be written; parts are named `<read>_segment_<n>`.
A `window` of 0 is treated as 1.
"""
struct SplitByLowQuality <: TrimStrategy
    cutoff::UInt8
    window::Int
    function SplitByLowQuality(cutoff::Integer; window::Integer = 1)
        window >= 0 || throw(ArgumentError("window must be non-negative"))
        return new(check_cutoff(cutoff), max(window, 1))
    end
end

function trim_ranges!(out::Vector{UnitRange{Int}}, s::SplitByLowQuality, qual::AbstractVector{UInt8})
    empty!(out)
    # `start` is the first good base of the open segment (0: none open),
    # `last_good` the latest good base, `bad_run` the bad bases since then.
    start = 0
    last_good = 0
    bad_run = 0
    @inbounds for i in eachindex(qual)
        if qual[i] - PHRED_OFFSET >= s.cutoff
            start == 0 && (start = i)
            last_good = i
            bad_run = 0
        else
            bad_run += 1
            if bad_run >= s.window
                start > 0 && push!(out, start:last_good)
                start = 0
                last_good = 0
            end
        end
    end
    start > 0 && push!(out, start:last_good)
    return out
end
