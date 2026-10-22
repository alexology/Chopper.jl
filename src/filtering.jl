# Read filters and the decision for one read.

"""
    ReadStatus

Outcome of the filters for one read: `PASSED`, or the first filter it failed
(`FAILED_LENGTH`, `FAILED_QUALITY`, `FAILED_GC`, `CONTAMINANT`). Filters run
in that order, cheapest first.
"""
@enum ReadStatus::UInt8 PASSED FAILED_LENGTH FAILED_QUALITY FAILED_GC CONTAMINANT

"""
    ReadFilter(; min_quality = 0, max_quality = 1000, min_length = 1,
               max_length = typemax(Int), min_gc = nothing, max_gc = nothing,
               contaminants = nothing, contaminant_fraction = 0.9)

Filters applied to whole, untrimmed reads, as chopper's filtering options:

| keyword | chopper | |
|---|---|---|
| `min_quality`, `max_quality` | `-q`, `--maxqual` | mean quality ([`mean_quality`](@ref)) |
| `min_length`, `max_length` | `-l`, `--maxlength` | read length; `min_length` is checked again on trimmed segments |
| `min_gc`, `max_gc` | `--mingc`, `--maxgc` | GC fraction, 0–1 ([`gc_content`](@ref)) |
| `contaminants` | `-c` | a [`ContaminantIndex`](@ref) or the path of a FASTA file |
| `contaminant_fraction` | (fixed 0.9) | see [`is_contaminant`](@ref) |

As in chopper, the quality filter is evaluated only if `min_quality > 0` or
`max_quality < 1000`, and the GC filter only if `min_gc` or `max_gc` is set:
with the defaults, reads of any quality (or GC content) pass.
"""
struct ReadFilter
    min_quality::Float64
    max_quality::Float64
    min_length::Int
    max_length::Int
    min_gc::Float64
    max_gc::Float64
    check_quality::Bool
    check_gc::Bool
    contaminants::Union{Nothing,ContaminantIndex}
    contaminant_fraction::Float64
end

function ReadFilter(; min_quality::Real = 0.0, max_quality::Real = 1000.0,
                    min_length::Integer = 1, max_length::Integer = typemax(Int),
                    min_gc::Union{Nothing,Real} = nothing, max_gc::Union{Nothing,Real} = nothing,
                    contaminants::Union{Nothing,AbstractString,ContaminantIndex} = nothing,
                    contaminant_fraction::Real = 0.9)
    min_length >= 0 || throw(ArgumentError("min_length must be non-negative"))
    max_length >= 0 || throw(ArgumentError("max_length must be non-negative"))
    for (name, v) in (("min_gc", min_gc), ("max_gc", max_gc))
        v === nothing || 0 <= v <= 1 ||
            throw(ArgumentError("$name must be a fraction from 0 to 1, got $v"))
    end
    0 < contaminant_fraction <= 1 ||
        throw(ArgumentError("contaminant_fraction must be in (0, 1], got $contaminant_fraction"))
    index = contaminants isa AbstractString ? ContaminantIndex(contaminants) : contaminants
    return ReadFilter(min_quality, max_quality, min_length, max_length,
                      something(min_gc, 0.0), something(max_gc, 1.0),
                      min_quality > 0 || max_quality < 1000, min_gc !== nothing || max_gc !== nothing,
                      index, contaminant_fraction)
end

"""
    read_status(filter::ReadFilter, seq, qual) -> ReadStatus

First filter of `filter` that the read fails, or `PASSED`. `seq` and `qual`
are the sequence and the raw Phred+33 quality line (bytes or strings).
"""
read_status(f::ReadFilter, seq, qual) =
    read_status(f, as_bytes(seq), as_bytes(qual))

as_bytes(x::AbstractString) = codeunits(x)
as_bytes(x::AbstractVector{UInt8}) = x

function read_status(f::ReadFilter, seq::AbstractVector{UInt8}, qual::AbstractVector{UInt8})
    n = length(seq)
    f.min_length <= n <= f.max_length || return FAILED_LENGTH
    if f.check_quality
        q = mean_quality(qual)
        f.min_quality <= q <= f.max_quality || return FAILED_QUALITY
    end
    if f.check_gc
        gc = gc_content(seq)
        f.min_gc <= gc <= f.max_gc || return FAILED_GC
    end
    if f.contaminants !== nothing
        is_contaminant(f.contaminants, seq; min_fraction = f.contaminant_fraction) && return CONTAMINANT
    end
    return PASSED
end

# Everything that decides the fate of a read.
struct ChopperConfig{T<:TrimStrategy}
    filter::ReadFilter
    trim::T
    inverse::Bool
    update_mods::Bool
end

"""
    read_segments(seq, qual; filter = ReadFilter(), trim = NoTrim(), inverse = false)
        -> Vector{UnitRange{Int}}

Ranges of the read that chopper would write: empty if the read is discarded,
one range per written record otherwise (several when `trim` splits it).

The filters judge the whole read; then the trimming strategy is applied and
each trimmed segment shorter than `filter.min_length` is dropped (not with
`NoTrim`). With `inverse = true` the reads that fail a filter are kept and
trimmed instead, as chopper's `--inverse`.
"""
function read_segments(seq, qual; filter::ReadFilter = ReadFilter(), trim::TrimStrategy = NoTrim(),
                       inverse::Bool = false)
    s, q = as_bytes(seq), as_bytes(qual)
    length(s) == length(q) || throw(DimensionMismatch("sequence and quality differ in length"))
    out = UnitRange{Int}[]
    read_segments!(out, ChopperConfig(filter, trim, inverse, false), s, q)
    return out
end

# Fills `out` with the segments to write and returns the filter status.
function read_segments!(out::Vector{UnitRange{Int}}, cfg::ChopperConfig, seq::AbstractVector{UInt8},
                        qual::AbstractVector{UInt8})
    # Filters stop at the first failure also with `inverse`: the read is
    # then kept whatever the remaining filters say. (chopper evaluates all of
    # them in that case, including the alignment against contaminants.)
    status = read_status(cfg.filter, seq, qual)
    if (status == PASSED) == cfg.inverse
        empty!(out)
        return status
    end
    trim_ranges!(out, cfg.trim, qual)
    drop_short_segments!(out, cfg.trim, cfg.filter.min_length)
    return status
end

drop_short_segments!(out, ::NoTrim, min_length) = out
drop_short_segments!(out, ::TrimStrategy, min_length) = filter!(r -> length(r) >= min_length, out)
