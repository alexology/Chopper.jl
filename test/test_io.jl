# Compressed input and FASTQ edge cases.

@testset "Compressed input" begin
    raw = read(joinpath(DATA, "test.fastq"))
    expected = read(joinpath(REF, "q10_l500.fastq"))
    half = findall(==(UInt8('\n')), raw)[400]          # end of record 100
    mktempdir() do dir
        files = Dict(
            "plain.fq" => raw,
            "a.gz" => transcode(GzipCompressor, raw),
            # Multi-member gzip, as written by pigz/bgzip or `cat a.gz b.gz`.
            "multi.gz" => vcat(transcode(GzipCompressor, raw[1:half]), transcode(GzipCompressor, raw[(half + 1):end])),
            "noext" => transcode(GzipCompressor, raw),   # recognised by content, not name
            "a.bz2" => transcode(Bzip2Compressor, raw),
            "a.xz" => transcode(XzCompressor, raw),
            "a.zst" => transcode(ZstdCompressor, raw),
        )
        for (name, bytes) in files
            path = joinpath(dir, name)
            write(path, bytes)
            out = IOBuffer()
            s = chopper(path, out; min_quality = 10, min_length = 500)
            @test take!(out) == expected
            @test s.reads == 250
            # From a stream (stdin) as well.
            out = IOBuffer()
            open(path) do io
                chopper(io, out; min_quality = 10, min_length = 500)
            end
            @test take!(out) == expected
        end
    end
    @test Chopper.sniff_format(UInt8[]) == :none
    @test Chopper.sniff_format(codeunits("BZh91")) == :bzip2
end

@testset "FASTQ edge cases" begin
    # Empty input, empty reads, CRLF line ends, lower-case bases.
    @test chopper(IOBuffer(""), IOBuffer()).reads == 0
    out = IOBuffer()
    s = chopper(IOBuffer("@e\n\n+\n\n@r desc\r\nacgtNN\r\n+\r\nIIII##\r\n"), out)
    @test s.reads == 2 && s.failed_length == 1           # chopper itself panics on empty reads
    @test String(take!(out)) == "@r desc\nacgtNN\n+\nIIII##\n"
    out = IOBuffer()
    chopper(IOBuffer("@e\n\n+\n\n"), out; min_length = 0)
    @test String(take!(out)) == "@e\n\n+\n\n"
    @test_throws Exception chopper(IOBuffer("@r\nACGT\n+\nIII\n"), IOBuffer())
    @test_throws Exception chopper(IOBuffer("r\nACGT\n+\nIIII\n"), IOBuffer())
end

# A source that fails after its first byte.
mutable struct FailingIO <: IO
    reads::Int
end
Base.eof(io::FailingIO) = io.reads > 0 ? error("injected read failure") : false
Base.read(io::FailingIO, ::Integer) = (io.reads += 1; UInt8[0x41])

@testset "Prefetch errors and output safety" begin
    # The producer's error always reaches the reader (it used to look like
    # the end of the input when the channel closed before the task failed)
    for _ in 1:500
        p = Chopper.prefetch(FailingIO(0); chunk_size = 1, depth = 1)
        @test read(p, UInt8) == 0x41
        @test_throws ErrorException eof(p)
    end
    # The input file is never truncated by using it as output
    mktempdir() do dir
        p = joinpath(dir, "reads.fastq")
        write(p, read(joinpath(DATA, "test.fastq")))
        n = filesize(p)
        @test_throws ArgumentError chopper(p, p)
        @test_throws ArgumentError chopper(p, joinpath(dir, ".", "READS.fastq"))
        @test filesize(p) == n
        # A failed run leaves an existing output unchanged and no temporary file
        bad = joinpath(dir, "bad.fastq")
        write(bad, read(p, String) * "@r\nACGT\n+\nIII\n")
        out = joinpath(dir, "out.fastq")
        write(out, "previous")
        for threads in (1, 4)
            @test_throws Exception chopper(bad, out; threads)
            @test read(out, String) == "previous"
        end
        @test sort(readdir(dir)) == ["bad.fastq", "out.fastq", "reads.fastq"]
        # A successful run replaces it with the full output
        expected = IOBuffer()
        chopper(p, expected; min_length = 500)
        chopper(p, out; min_length = 500)
        @test read(out) == take!(expected)
        # A missing output directory is reported before anything is written
        @test_throws ArgumentError chopper(p, joinpath(dir, "missing", "out.fastq"))
        @test_throws ArgumentError chopper(joinpath(dir, "none.fastq"), joinpath(dir, "x.fastq"))
        @test sort(readdir(dir)) == ["bad.fastq", "out.fastq", "reads.fastq"]
    end
end

@static if Sys.iswindows()
    const GENERIC_READ = 0x80000000
    # Opens `path` for reading, sharing it with readers only (no rename or delete).
    lock_file(path) = ccall((:CreateFileW, "kernel32"), stdcall, Ptr{Cvoid},
        (Cwstring, UInt32, UInt32, Ptr{Cvoid}, UInt32, UInt32, Ptr{Cvoid}),
        path, GENERIC_READ, 0x00000001, C_NULL, 0x00000003, 0x00000000, C_NULL)
    unlock_file(handle) = ccall((:CloseHandle, "kernel32"), stdcall, Cint, (Ptr{Cvoid},), handle)
end

@testset "Publishing the output" begin
    mktempdir() do dir
        path = joinpath(dir, "out.fastq")
        # A failed rename never removes what is at the destination (`mv` with
        # `force` deleted it and fell back to copying)
        mkdir(path)
        target = Chopper.OutputTarget(joinpath(dir, "new.fastq"))
        write(target.temp, "replacement")
        target = Chopper.OutputTarget(path, target.temp)
        write(joinpath(path, "kept"), "kept")
        @test_throws Base.IOError Chopper.publish(target)
        @test read(joinpath(path, "kept"), String) == "kept"
        Chopper.discard(target)
        @test readdir(dir) == ["out.fastq"]
        rm(path; recursive = true)
        @static if Sys.iswindows()
            # The temporary file or the output locked by another process
            for locked in (:temp, :output)
                write(path, "previous")
                target = Chopper.OutputTarget(path)
                write(target.temp, "replacement")
                handle = lock_file(locked === :temp ? target.temp : path)
                @test handle != Ptr{Cvoid}(-1)
                try
                    @test_throws Base.IOError Chopper.publish(target)
                    @test read(path, String) == "previous"
                    # Cleanup after the failure does not hide it
                    @test Chopper.discard(target) === nothing
                finally
                    unlock_file(handle)
                end
                Chopper.discard(target)
                @test readdir(dir) == ["out.fastq"]
            end
        end
        @static if Sys.isunix()
            # A new output gets the default permissions, not those of a private
            # temporary file; a replaced one keeps its own
            rm(path; force = true)
            probe = joinpath(dir, "probe")
            touch(probe)
            target = Chopper.OutputTarget(path)
            Chopper.publish(target)
            @test filemode(path) & 0o777 == filemode(probe) & 0o777
            chmod(path, 0o640)
            target = Chopper.OutputTarget(path)
            write(target.temp, "replacement")
            Chopper.publish(target)
            @test filemode(path) & 0o777 == 0o640
            @test read(path, String) == "replacement"
        end
    end
end

# A source that never ends.
struct EndlessIO <: IO end
Base.eof(::EndlessIO) = false
Base.read(::EndlessIO, n::Integer) = fill(UInt8('A'), n)

@testset "Prefetch cancellation" begin
    # Closing the reader stops the producer without an error
    for _ in 1:100
        p = Chopper.prefetch(EndlessIO(); chunk_size = 4, depth = 1)
        @test read(p, UInt8) == UInt8('A')
        close(p)
        wait(p.task)
        @test istaskdone(p.task) && !istaskfailed(p.task)
    end
end

using Random

# Records the decompressor streams that `chopper` creates.
function with_gzip_streams(f)
    original = Chopper.DECOMPRESSORS[:gzip]
    streams = Any[]
    Chopper.DECOMPRESSORS[:gzip] = io -> push!(streams, original(io))[end]
    try
        f(streams)
    finally
        Chopper.DECOMPRESSORS[:gzip] = original
    end
end

# An output whose writes fail, after signalling that they were reached.
struct BrokenOutput <: IO
    reached::Channel{Nothing}
end
Base.write(io::BrokenOutput, ::Vector{UInt8}) = (put!(io.reached, nothing); error("injected output failure"))
Base.flush(::BrokenOutput) = nothing

@testset "Release of the input" begin
    released(stream) = stream.codec.zstream.state == C_NULL
    good = transcode(GzipCompressor, codeunits("@r\nACGT\n+\nIIII\n"))
    bad = transcode(GzipCompressor, codeunits("@r\nACGT\n+\nIII\n"))
    with_gzip_streams() do streams
        # The decompressor of a caller's stream is freed (its native memory has
        # no finalizer), and the stream stays open
        for threads in (1, 4), (data, ok) in ((good, true), (bad, false))
            input = IOBuffer(data)
            if ok
                @test chopper(input, IOBuffer(); threads).reads == 1
            else
                @test_throws Exception chopper(input, IOBuffer(); threads)
            end
            @test released(pop!(streams))
            @test isopen(input)
        end
        # A file opened by chopper is closed
        mktempdir() do dir
            path = joinpath(dir, "in.fq.gz")
            write(path, good)
            for threads in (1, 4)
                chopper(path, IOBuffer(); threads)
                stream = pop!(streams)
                @test released(stream) && !isopen(stream.stream)
            end
            rm(path)
        end
        # A failure does not wait for a read of a caller's stream that waits for
        # data: the read-ahead task stops, and frees the decompressor, when the
        # stream is closed
        if Threads.nthreads() > 1
            raw = join("@r$i\n" * String(rand(MersenneTwister(i), collect(codeunits("ACGT")), 10_000)) *
                       "\n+\n" * repeat("I", 10_000) * "\n" for i in 1:150)
            source = Base.BufferStream()
            write(source, transcode(GzipCompressor, codeunits(raw)))
            output = BrokenOutput(Channel{Nothing}(1))
            task = Threads.@spawn try
                chopper(source, output; threads = 2, batch_reads = 1)
            catch e
                e
            end
            try
                @test timedwait(() -> istaskdone(task), 30) == :ok
                @test fetch(task) == ErrorException("injected output failure")
                @test isopen(source)
            finally
                close(source)
            end
            stream = pop!(streams)
            @test timedwait(() -> released(stream), 10) == :ok
        end
    end
end
