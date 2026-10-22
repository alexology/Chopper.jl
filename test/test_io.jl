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
