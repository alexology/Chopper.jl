"""
    Chopper

Filtering and trimming of long-read FASTQ files (Oxford Nanopore, PacBio): a
native Julia port of the Rust tool [chopper](https://github.com/wdecoster/chopper)
0.14.1 by Wouter De Coster.

Reads are filtered on mean quality, length, GC content and similarity to
contaminant sequences, then optionally trimmed with one of four strategies
([`FixedCrop`](@ref), [`TrimByQuality`](@ref), [`BestReadSegment`](@ref),
[`SplitByLowQuality`](@ref)). The main entry point is [`chopper`](@ref); the
same filters are available per read through [`read_segments`](@ref).
"""
module Chopper

import FASTX
import Minimap2
using FASTX: FASTQ, FASTA, description, sequence, quality, seqsize
using StringViews: StringView
using TranscodingStreams: TranscodingStreams, TranscodingStream, NoopStream
using CodecZlibNG: GzipDecompressorStream, GzipCompressorStream

export chopper, ChopperStats,
       ReadFilter, read_segments, read_status, ReadStatus,
       TrimStrategy, NoTrim, FixedCrop, TrimByQuality, BestReadSegment,
       SplitByLowQuality, trim_ranges, trim_ranges!,
       mean_quality, gc_content, phred_to_probability,
       ContaminantIndex, is_contaminant, best_match_length,
       rewrite_mod_tags, ModTagError,
       open_fastq

include("quality.jl")
include("trimming.jl")
include("modtags.jl")
include("header.jl")
include("contaminants.jl")
include("filtering.jl")
include("io.jl")
include("pipeline.jl")
include("cli.jl")
include("precompile.jl")

end # module
