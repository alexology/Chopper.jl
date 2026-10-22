# FASTQ header of a written segment.
#
# The header (the line after '@') is: the read name, up to the first space or
# tab; an optional free-text description after a space; SAM tags as
# tab-separated fields after the first tab (`samtools fastq -T ...`, dorado).
# A split read gets the suffix `_segment_<n>` on the name, never at the end of
# the line, where it would corrupt the last tag. FASTX keeps the whole header
# line, so it never needs to be reassembled (chopper's FASTQ parser splits it
# at the first space, which may fall inside a tag value).

const ASCII_TRAILING_SPACE = (UInt8(' '), UInt8('\t'), UInt8('\v'), UInt8('\f'), UInt8('\r'))

# chopper's FASTQ parser trims trailing whitespace from the header line; the
# same is done here so that headers are written byte for byte as chopper does.
function header_end(h::AbstractVector{UInt8})
    n = length(h)
    @inbounds while n > 0 && h[n] in ASCII_TRAILING_SPACE
        n -= 1
    end
    return n
end

# Writes '@' and the header of segment `index` of `nsegments` for the range
# `start+1:stop` of `seq`. Returns whether the read carried an MM tag (only
# looked at when `update_mods`).
function write_header(io::IO, header::AbstractVector{UInt8}, seq::AbstractVector{UInt8},
                      start::Int, stop::Int, nsegments::Int, index::Int, update_mods::Bool)
    n = header_end(header)
    write(io, UInt8('@'))
    if nsegments == 1 && !update_mods
        write(io, view(header, 1:n))
        return false
    end
    tab = findfirst(==(UInt8('\t')), view(header, 1:n))
    name_end = tab === nothing ? n : tab - 1
    space = findfirst(==(UInt8(' ')), view(header, 1:name_end))
    name_stop = space === nothing ? name_end : space - 1
    write(io, view(header, 1:name_stop))
    nsegments > 1 && print(io, "_segment_", index)
    write(io, view(header, (name_stop + 1):name_end))       # " description", if any
    tab === nothing && return false
    if !update_mods
        write(io, view(header, tab:n))
        return false
    end
    tags = StringView(view(header, (tab + 1):n))
    fields, has_mm = rewrite_mod_fields(split(tags, '\t'), seq, start, stop)
    for f in fields
        write(io, UInt8('\t'))
        write(io, f)
    end
    return has_mm
end
