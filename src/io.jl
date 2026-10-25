# Input and output streams, with transparent decompression.

# Decompressor stream constructors by format. gzip is always available;
# the package extensions for CodecBzip2, CodecXz and CodecZstd register theirs.
const DECOMPRESSORS = Dict{Symbol,Any}(:gzip => GzipDecompressorStream)

const CODEC_PACKAGES = Dict(
    :bzip2 => Base.PkgId(Base.UUID("523fee87-0ab8-5b00-afb7-3ecf72e48cfd"), "CodecBzip2"),
    :xz => Base.PkgId(Base.UUID("ba30903b-d9e8-5048-a5ec-d1f5b0d4b47b"), "CodecXz"),
    :zstd => Base.PkgId(Base.UUID("6b39b394-51ab-5f42-8807-6242bab2b4c2"), "CodecZstd"),
)

# Compression format from the first bytes of a stream (magic numbers).
function sniff_format(magic::AbstractVector{UInt8})
    starts(m) = length(magic) >= length(m) && all(magic[i] == m[i] for i in eachindex(m))
    starts((0x1f, 0x8b)) && return :gzip
    starts(codeunits("BZh")) && return :bzip2
    starts((0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00)) && return :xz
    starts((0x28, 0xb5, 0x2f, 0xfd)) && return :zstd
    return :none
end

# The decompressor of `format`, loading its codec package (if installed in
# the active environment) the first time it is needed.
function decompressor(format::Symbol)
    haskey(DECOMPRESSORS, format) && return DECOMPRESSORS[format]
    pkg = CODEC_PACKAGES[format]
    try
        Base.require(pkg)
    catch
        throw(ArgumentError("$format-compressed input needs the package $(pkg.name): " *
                            "install it with `import Pkg; Pkg.add(\"$(pkg.name)\")`"))
    end
    haskey(DECOMPRESSORS, format) ||
        error("the extension of Chopper for $(pkg.name) did not load")
    return DECOMPRESSORS[format]
end

"""
    open_fastq(path::AbstractString) -> TranscodingStream
    open_fastq(io::IO) -> TranscodingStream

Open a FASTQ file or stream for reading, decompressing it if needed. The
format is recognised from the first bytes, not from the file name: gzip
(always available, also multi-member files such as those of `pigz` or
`bgzip`), bzip2, xz and zstd (with the packages CodecBzip2, CodecXz and
CodecZstd installed). Closing the returned stream closes `path`'s file, or
`io`.
"""
function open_fastq(path::AbstractString)
    isfile(path) || throw(ArgumentError("input file $path is invalid"))
    file = open(path, "r")
    try
        return open_fastq(file)
    catch
        close(file)
        rethrow()
    end
end

function open_fastq(io::IO)
    stream = NoopStream(io)
    magic = read(stream, 6)
    TranscodingStreams.unread(stream, magic)
    format = sniff_format(magic)
    format === :none && return stream
    return Base.invokelatest(decompressor(format), stream)
end

# Ends the use of a stream from `open_fastq`. If it wraps a caller's IO, only
# the decompressor is freed (its native memory has no finalizer) and the IO
# stays open; otherwise the stream is closed with the file below it.
function release(stream::TranscodingStream; owned::Bool)
    owned ? close(stream) : TranscodingStreams.finalize(stream.codec)
    return nothing
end

# Decompression in a separate task, a few chunks ahead of the reader, so that
# it overlaps with filtering; gzip decompression is otherwise the slowest
# step. `PrefetchedInput` is a minimal read-only IO over the chunks.
#
# The producer ends in one of three ways, kept distinct: the end of the
# source (the task succeeds after the last chunk), a read error (the task
# fails with it and the reader rethrows it) or a cancellation by `close`
# from the reader (the task succeeds without reading further; errors of the
# read in progress are ignored). The source is read only by the producer,
# which calls `done(source)` when it ends: a read that never returns cannot
# be interrupted, and then the reader must not wait for it.
mutable struct PrefetchedInput <: IO
    chunks::Channel{Vector{UInt8}}
    current::Vector{UInt8}
    pos::Int
    @atomic cancelled::Bool
    task::Task
    PrefetchedInput(chunks) = new(chunks, UInt8[], 1, false)
end

function prefetch(source::IO, done = identity; chunk_size::Int = 1 << 20, depth::Int = 8)
    chunks = Channel{Vector{UInt8}}(depth)
    input = PrefetchedInput(chunks)
    input.task = Threads.@spawn begin
        try
            while !(@atomic input.cancelled) && !eof(source)
                put!(chunks, read(source, chunk_size))
            end
        catch
            # `put!` on the channel closed by the reader, or a read cut short
            # after the cancellation (a caller closing its stream)
            (@atomic input.cancelled) || rethrow()
        finally
            close(chunks)
            done(source)
        end
    end
    return input
end

# The channel is closed in the producer's `finally`, before the task is marked
# as failed: wait for the task before reporting the end of the input.
function _producer_done(io::PrefetchedInput)
    try
        wait(io.task)
    catch e
        throw(e isa TaskFailedException ? e.task.exception : e)
    end
    return true
end

function Base.eof(io::PrefetchedInput)
    while io.pos > length(io.current)
        !isready(io.chunks) && !isopen(io.chunks) && return _producer_done(io)
        chunk = try
            take!(io.chunks)
        catch e
            e isa InvalidStateException || rethrow()
            return _producer_done(io)
        end
        io.current, io.pos = chunk, 1
    end
    return false
end

Base.bytesavailable(io::PrefetchedInput) = length(io.current) - io.pos + 1

function Base.unsafe_read(io::PrefetchedInput, p::Ptr{UInt8}, n::UInt)
    done = 0
    while done < n
        eof(io) && throw(EOFError())
        k = min(Int(n) - done, bytesavailable(io))
        GC.@preserve io unsafe_copyto!(p + done, pointer(io.current, io.pos), k)
        io.pos += k
        done += k
    end
    return nothing
end

function Base.read(io::PrefetchedInput, ::Type{UInt8})
    eof(io) && throw(EOFError())
    b = io.current[io.pos]
    io.pos += 1
    return b
end

Base.isopen(io::PrefetchedInput) = !eof(io)
function Base.close(io::PrefetchedInput)
    @atomic io.cancelled = true
    close(io.chunks)
    return nothing
end

# Output stream for a path: gzip-compressed if the name ends in `.gz`.
function open_output(path::AbstractString; gzip::Bool = endswith(lowercase(path), ".gz"))
    io = open(path, "w")
    return gzip ? GzipCompressorStream(io) : io
end

# Output to a path is written to a temporary file in the same directory and
# renamed over `path` only when the run succeeded, so a failed run leaves
# any previous file untouched and no partial one. Paths that exist but are
# not regular files (`/dev/null`, `NUL`, named pipes) are written directly.
struct OutputTarget
    path::String
    temp::Union{Nothing,String}
end

function OutputTarget(path::AbstractString)
    path = String(path)
    dir = dirname(abspath(path))
    isdir(dir) || throw(ArgumentError("the output directory does not exist: $dir"))
    ispath(path) && !isfile(path) && return OutputTarget(path, nothing)
    temp = new_temp(dir)
    # Keep the permissions of a file being replaced, as writing it in place would
    isfile(path) && Sys.isunix() && chmod(temp, filemode(path) & 0o7777)
    return OutputTarget(path, temp)
end

# An empty new file in `dir`, created with the permissions `open(path, "w")`
# gives a new file (0o666 less the umask); `mktemp` would make it private.
function new_temp(dir::AbstractString)
    flags = Base.Filesystem.JL_O_WRONLY | Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL
    while true
        temp = tempname(dir; cleanup = false)
        try
            close(Base.Filesystem.open(temp, flags, 0o666))
            return temp
        catch e
            e isa Base.IOError && e.code == Base.UV_EEXIST || rethrow()
        end
    end
end

open_output(t::OutputTarget) =
    open_output(something(t.temp, t.path); gzip = endswith(lowercase(t.path), ".gz"))

# A single rename, which replaces `path` atomically or fails leaving it as it
# was. `mv(...; force = true)` is not used: after a failed rename it deletes
# `path` and falls back to copying, so an error could still replace the file.
publish(t::OutputTarget) = (t.temp === nothing || Base.Filesystem.rename(t.temp, t.path); nothing)

# Called after a failure: a temporary file that cannot be removed (still
# locked by another process) is left behind rather than hiding that failure.
function discard(t::OutputTarget)
    t.temp === nothing && return nothing
    try
        rm(t.temp; force = true)
    catch e
        e isa Base.IOError || rethrow()
    end
    return nothing
end
