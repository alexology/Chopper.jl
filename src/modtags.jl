# Rewriting of base modification tags (MM/ML/MN) for trimmed reads.
#
# `samtools fastq -T MM,ML,MN` writes SAM tags as tab-separated fields after
# the read name. MM encodes the positions of modified bases as counts of
# canonical bases to skip, so trimming or splitting a read invalidates it;
# this file recomputes MM, ML and MN for a subsequence. The rules follow
# htslib (`sam_mods.c`), as in chopper:
#
# - FASTQ is always in the original read orientation, so canonical bases are
#   counted left to right and the `+`/`-` strand of MM is carried through;
# - the canonical base `N` matches every base, not only literal `N`s;
# - ML is optional; a list with several codes (`C+mh`) has several ML values
#   per call, and the ML values of all lists are concatenated in list order.
#
# Tags that cannot be rewritten raise a `ModTagError`: writing them through
# would attach modification coordinates to a sequence they do not describe.

"""
    ModTagError

Base modification tags of a read that cannot be rewritten: malformed `MM`,
`ML` or `MN` tags, `MM` and `ML` with different numbers of calls, an `MM` tag
that refers to more bases than the read has, an `MN` tag with a different
read length, or `ML` without `MM`. Fields: `kind`, one of `:malformed_mm`,
`:malformed_ml`, `:malformed_mn`, `:count_mismatch`, `:beyond_sequence`,
`:stale_mn`, `:ml_without_mm`; `msg`; `read`, the read name (empty if not
known).
"""
struct ModTagError <: Exception
    kind::Symbol
    msg::String
    read::String
end
ModTagError(kind::Symbol, msg::AbstractString) = ModTagError(kind, msg, "")

Base.showerror(io::IO, e::ModTagError) =
    print(io, isempty(e.read) ? e.msg : "read $(e.read): $(e.msg)")

# One MM sublist, e.g. `C+mh?,5,12`.
struct ModList
    prefix::SubString{String}    # canonical base, strand, codes, optional flag
    base::UInt8                  # canonical base to count (U counted as T)
    stride::Int                  # ML values per modified position
    deltas::Vector{Int}
end

# Unsigned decimal integer as parsed by Rust (`str::parse::<usize>`): an
# optional `+`, then at least one digit; no spaces, no sign `-`.
function parse_unsigned(s::AbstractString, max::Integer = typemax(Int))
    b = codeunits(s)
    i = 1
    i <= length(b) && b[i] == UInt8('+') && (i += 1)
    i <= length(b) || return nothing
    v = 0
    for j in i:length(b)
        d = b[j] - UInt8('0')
        d <= 9 || return nothing
        v > (max - d) ÷ 10 && return nothing
        v = 10v + d
    end
    return v
end

function parse_mm(mm::AbstractString)
    lists = ModList[]
    for chunk in eachsplit(String(mm), ';')
        isempty(chunk) && continue
        comma = findfirst(==(','), chunk)
        head = comma === nothing ? chunk : SubString(chunk, 1, prevind(chunk, comma))
        deltas_str = comma === nothing ? "" : SubString(chunk, nextind(chunk, comma))
        bytes = codeunits(head)
        length(bytes) < 3 && throw(ModTagError(:malformed_mm, "malformed MM tag: '$chunk' is too short"))
        bytes[2] in (UInt8('+'), UInt8('-')) ||
            throw(ModTagError(:malformed_mm, "malformed MM tag: '$head' has no '+' or '-' strand"))
        base = bytes[1] == UInt8('U') ? UInt8('T') : bytes[1]
        base in (UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T'), UInt8('N')) ||
            throw(ModTagError(:malformed_mm, "malformed MM tag: '$(Char(bytes[1]))' is not a canonical base"))
        # The codes run to the end of the head, minus the optional '.'/'?' flag.
        codes = SubString(head, 3)
        if endswith(codes, '.') || endswith(codes, '?')
            codes = chop(codes)
        end
        # A numeric ChEBI id is a single modification; letters are one each.
        stride = if !isempty(codes) && isdigit(first(codes))
            all(isdigit, codes) || throw(ModTagError(:malformed_mm,
                "malformed MM tag: '$codes' is neither a ChEBI id nor a list of modification codes"))
            1
        else
            (!isempty(codes) && all(c -> isascii(c) && isletter(c), codes)) ||
                throw(ModTagError(:malformed_mm, "malformed MM tag: '$head' has no modification code"))
            ncodeunits(codes)
        end
        deltas = Int[]
        if !isempty(deltas_str)
            for d in eachsplit(deltas_str, ',')
                v = parse_unsigned(d)
                v === nothing && throw(ModTagError(:malformed_mm, "malformed MM tag: '$d' is not a valid skip count"))
                push!(deltas, v)
            end
        end
        push!(lists, ModList(head, base, stride, deltas))
    end
    return lists
end

function parse_ml(ml::AbstractString)
    values = UInt8[]
    isempty(ml) && return values
    for v in eachsplit(ml, ',')
        x = parse_unsigned(v, 255)
        x === nothing && throw(ModTagError(:malformed_ml, "malformed ML tag: '$v' is not a value in 0-255"))
        push!(values, x)
    end
    return values
end

# MM and ML for the subsequence seq[start+1:stop] (0-based half-open
# [start, stop) internally, as in chopper). Lists that lose all their calls
# are kept as empty lists, so that "this modification was called and none
# was found here" survives.
function subset_mods(mm::AbstractString, ml::Union{Nothing,Vector{UInt8}},
                     seq::AbstractVector{UInt8}, start::Int, stop::Int)
    lists = parse_mm(mm)
    calls = sum((l.stride * length(l.deltas) for l in lists); init = 0)
    if ml !== nothing && calls != length(ml)
        throw(ModTagError(:count_mismatch,
            "MM lists $calls modification call(s) but ML holds $(length(ml)) probability value(s)"))
    end
    new_mm = IOBuffer()
    new_ml = UInt8[]
    ml_cursor = 0
    n = length(seq)
    for list in lists
        print(new_mm, list.prefix)
        pos = 0               # next base to examine (0-based)
        window_occ = 0        # canonical bases seen so far inside [start, stop)
        last_kept_occ = 0     # window_occ at the last kept call
        for delta in list.deltas
            remaining = delta
            p = -1
            i = pos
            @inbounds while i < n
                b = seq[i + 1]
                if list.base == UInt8('N') || uppercase_byte(b) == list.base
                    start <= i < stop && (window_occ += 1)
                    if remaining == 0
                        p = i
                        break
                    end
                    remaining -= 1
                end
                i += 1
            end
            p < 0 && throw(ModTagError(:beyond_sequence,
                "MM tag refers to more '$(Char(list.base))' bases than the sequence contains; " *
                "the tags do not describe this read"))
            pos = p + 1
            if start <= p < stop
                # window_occ counts p itself, so the new skip count is the
                # number of occurrences between the previous kept call and p.
                print(new_mm, ',', window_occ - 1 - last_kept_occ)
                last_kept_occ = window_occ
                ml === nothing || append!(new_ml, view(ml, (ml_cursor + 1):(ml_cursor + list.stride)))
            end
            ml_cursor += list.stride
        end
        print(new_mm, ';')
    end
    return String(take!(new_mm)), new_ml
end

@inline uppercase_byte(b::UInt8) = UInt8('a') <= b <= UInt8('z') ? b - 0x20 : b

# Case-insensitive tag prefix, e.g. `MM:Z:` also matches `Mm:Z:`; returns the
# rest of the field or `nothing`.
function strip_tag(field::AbstractString, prefix::String)
    b = codeunits(field)
    p = codeunits(prefix)
    length(b) < length(p) && return nothing
    for i in eachindex(p)
        uppercase_byte(b[i]) == uppercase_byte(p[i]) || return nothing
    end
    return SubString(field, length(p) + 1)
end

"""
    rewrite_mod_tags(fields, seq, range) -> Vector{String}

Rewrite the base modification tags among the tab-separated SAM tag `fields`
of a read with sequence `seq` (bytes or string) so that they describe the
subsequence `seq[range]`:

- `MM:Z:` keeps only the calls inside `range`, with recomputed skip counts;
- `ML:B:C` keeps the probabilities of those calls;
- `MN:i:` becomes `length(range)`; it must equal the length of `seq`.

Other fields are returned unchanged, including tags that trimming also
invalidates (`qs`, `ns`, `ts`, `du`, ...). Tag names are matched without
regard to case (`Mm`/`Ml` of older basecallers). Throws [`ModTagError`](@ref)
if the tags are malformed or do not describe `seq`.
"""
function rewrite_mod_tags(fields::AbstractVector{<:AbstractString}, seq, range::AbstractUnitRange{<:Integer})
    s = seq isa AbstractString ? codeunits(seq) : seq
    (isempty(range) || (first(range) >= 1 && last(range) <= length(s))) ||
        throw(BoundsError(s, range))
    out, _ = rewrite_mod_fields(fields, s, first(range) - 1, last(range))
    return out
end

# Returns the rewritten fields and whether the read carried an MM tag.
function rewrite_mod_fields(fields::AbstractVector{<:AbstractString}, seq::AbstractVector{UInt8},
                            start::Int, stop::Int)
    out = String[String(f) for f in fields]
    i = findfirst(f -> strip_tag(f, "MN:i:") !== nothing, fields)
    if i !== nothing
        mn = parse_unsigned(strip_tag(fields[i], "MN:i:"))
        mn === nothing && throw(ModTagError(:malformed_mn, "malformed MN tag: '$(fields[i])' is not a length"))
        mn == length(seq) || throw(ModTagError(:stale_mn,
            "MN tag says the modification tags describe a $mn base sequence, but the read is $(length(seq)) bases"))
        out[i] = "MN:i:$(stop - start)"
    end
    mm_i = findfirst(f -> strip_tag(f, "MM:Z:") !== nothing, fields)
    if mm_i === nothing
        # ML without MM would otherwise be written through with stale
        # probabilities, which is the failure this rewriting prevents.
        any(f -> strip_tag(f, "ML:B:") !== nothing, fields) && throw(ModTagError(:ml_without_mm,
            "ML tag holds modification probabilities but there is no MM tag to place them"))
        return out, false
    end
    mm = strip_tag(fields[mm_i], "MM:Z:")
    ml_i = findfirst(f -> strip_tag(f, "ML:B:") !== nothing, fields)
    ml = nothing
    if ml_i !== nothing
        body = strip_tag(fields[ml_i], "ML:B:")
        startswith(body, 'C') ||
            throw(ModTagError(:malformed_ml, "malformed ML tag: '$(fields[ml_i])' is not a uint8 array"))
        values = SubString(body, 2)
        if startswith(values, ',')
            values = SubString(values, 2)
        elseif !isempty(values)
            throw(ModTagError(:malformed_ml,
                "malformed ML tag: '$(fields[ml_i])' is missing the comma after the array type"))
        end
        ml = parse_ml(values)
    end
    new_mm, new_ml = subset_mods(mm, ml, seq, start, stop)
    out[mm_i] = string(SubString(fields[mm_i], 1, 5), new_mm)
    if ml_i !== nothing
        io = IOBuffer()
        print(io, SubString(fields[ml_i], 1, 6))
        for v in new_ml
            print(io, ',', v)
        end
        out[ml_i] = String(take!(io))
    end
    return out, true
end
