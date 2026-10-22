# Contamination screen: does a read align to a contaminant sequence (e.g. the
# lambda phage DNA control strand of ONT kits) over most of its length?
#
# chopper maps every read with the minimap2 crate configured as
# `Aligner::builder().preset(LrHq).with_cigar().with_index(fasta)` and calls
# it a contaminant if any alignment has at least 90% of the read length as
# matching bases. The mapping is done here by Minimap2.jl, a port of the same
# minimap2 2.30 that gives the same alignments; only the decision is
# Chopper's.

"""
    ContaminantIndex(fasta; k = 19, w = 19)
    ContaminantIndex(names, sequences; k = 19, w = 19)

Minimap2 index of contaminant sequences (e.g. the lambda phage control DNA
of ONT kits, or a host genome) for [`is_contaminant`](@ref) and the
`contaminants` option of [`chopper`](@ref) and [`ReadFilter`](@ref). Reads
are mapped with the minimap2 preset `lr:hq` and base-level alignment, as
chopper does; `k` and `w` override the preset's minimizer length and window.

The index is only read while mapping and can be shared by any number of
tasks.
"""
struct ContaminantIndex
    names::Vector{String}
    index::Minimap2.MinimapIndex
    options::Minimap2.MapOptions
end

function ContaminantIndex(names::AbstractVector{<:AbstractString}, sequences::AbstractVector;
                          k::Integer = 19, w::Integer = 19)
    1 <= k <= 28 || throw(ArgumentError("k must be in 1:28, got $k"))
    1 <= w <= 255 || throw(ArgumentError("w must be in 1:255, got $w"))
    length(names) == length(sequences) || throw(ArgumentError("names and sequences differ in length"))
    isempty(sequences) && throw(ArgumentError("no contaminant sequences"))
    index_options, map_options = Minimap2.preset("lr:hq")
    index_options.k, index_options.w = k, w
    index = Minimap2.MinimapIndex(sequences, index_options; names)
    # As `with_index` of the minimap2 crate: options updated from the index.
    Minimap2.update_options!(map_options, index)
    map_options.flag |= Minimap2.F_CIGAR
    return ContaminantIndex(String.(names), index, map_options)
end

function ContaminantIndex(path::AbstractString; kwargs...)
    isfile(path) || throw(ArgumentError("contaminant FASTA file not found: $path"))
    names = String[]
    seqs = Vector{UInt8}[]
    open(FASTA.Reader, path) do reader
        for rec in reader
            push!(names, String(FASTX.identifier(rec)))
            push!(seqs, Vector{UInt8}(codeunits(sequence(rec))))
        end
    end
    isempty(seqs) && throw(ArgumentError("no sequences in contaminant FASTA file $path"))
    return ContaminantIndex(names, seqs; kwargs...)
end

Base.show(io::IO, idx::ContaminantIndex) = print(io, "ContaminantIndex(", length(idx.names),
    " sequences, ", sum(idx.index.lengths), " bases, k = ", idx.index.k, ", w = ", idx.index.w, ")")

"""
    best_match_length(index::ContaminantIndex, seq) -> Int

Largest number of matching bases ("match length", PAF column 10) over the
alignments of the read `seq` (bytes or string, either strand) to the
contaminant sequences of `index`; 0 if the read does not map.
"""
function best_match_length(idx::ContaminantIndex, seq)
    bytes = seq isa AbstractString ? codeunits(seq) : seq
    isempty(bytes) && return 0
    regions, _ = Minimap2.map_query(idx.index, bytes, idx.options)
    return isempty(regions) ? 0 : Int(maximum(r -> r.mlen, regions))
end

"""
    is_contaminant(index::ContaminantIndex, seq; min_fraction = 0.9) -> Bool

Whether the read `seq` aligns to a contaminant with at least
`min_fraction × length(seq)` matching bases, the criterion of chopper's
`--contam` filter.
"""
function is_contaminant(idx::ContaminantIndex, seq; min_fraction::Real = 0.9)
    bytes = seq isa AbstractString ? codeunits(seq) : seq
    isempty(bytes) && return false
    return best_match_length(idx, bytes) >= min_fraction * length(bytes)
end
