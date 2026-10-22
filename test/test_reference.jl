# Byte-for-byte comparison with the output of the original chopper 0.14.1
# (test/reference/, made by tools/make_reference.sh with one thread).

const CASES = let lines = readlines(joinpath(REF, "cases.tsv"))
    map(filter(!isempty, lines[2:end])) do l
        parts = split(l, '\t')
        (parts[1], parts[2], length(parts) >= 3 ? parts[3] : "")
    end
end

@testset "Same output as chopper 0.14.1" begin
    @test strip(read(joinpath(REF, "VERSION"), String)) == "chopper 0.14.1"
    for (name, input, arguments) in CASES
        @testset "$name" begin
            expected = read(joinpath(REF, "$name.fastq"))
            expected_err = read(joinpath(REF, "$name.stderr"), String)
            expected_status = parse(Int, readchomp(joinpath(REF, "$name.status")))
            args = vcat(["-i", joinpath(DATA, input)], split(arguments))
            for t in unique([1, Threads.nthreads()])
                # From the data folder, as the reference: `-c` paths are relative.
                status, out, err = cd(() -> run_cli(vcat(args, ["-t", string(t)])), DATA)
                @test status == expected_status
                @test out == expected
                @test err == expected_err
            end
        end
    end
end

@testset "Output independent of threads and batch size" begin
    path = joinpath(DATA, "test.fastq")
    reference = IOBuffer()
    s1 = chopper(path, reference; trim = SplitByLowQuality(12), min_length = 20, threads = 1)
    ref = take!(reference)
    for threads in (2, 4), batch_reads in (1, 7, 1000)
        io = IOBuffer()
        s = chopper(path, io; trim = SplitByLowQuality(12), min_length = 20, threads, batch_reads)
        @test take!(io) == ref
        @test s == s1
    end
    @test s1.reads == 250
    @test s1.segments == count(==(UInt8('\n')), ref) ÷ 4
end

@testset "Statistics" begin
    path = joinpath(DATA, "test.fastq")
    s = chopper(path, IOBuffer(); min_quality = 10, min_length = 500)
    @test s.reads == 250
    @test s.segments == s.kept_reads == s.passed == 205
    @test s.failed_length + s.failed_quality == 45
    @test s.bases == sum(length, readlines(path)[2:4:end])
    s = chopper(path, IOBuffer(); min_quality = 10, inverse = true)
    @test s.segments == 43 && s.failed_quality == 43 && s.passed == 207
    s = chopper(path, IOBuffer(); trim = FixedCrop(10_000, 0))
    @test s.trimmed_away == s.passed - s.kept_reads
    io = IOBuffer()
    show(io, MIME("text/plain"), s)
    @test occursin("kept $(s.segments) reads out of 250 reads", String(take!(io)))
    @test Chopper.thousands(1234567) == "1,234,567" && Chopper.thousands(999) == "999"
    @test_throws ArgumentError chopper(path, IOBuffer(); filter = ReadFilter(), min_length = 3)
    # Invalid batch sizes are rejected before the output is created
    mktempdir() do dir
        outpath = joinpath(dir, "out.fastq")
        for kw in ((batch_reads = 0,), (batch_bases = 0,), (batch_reads = -1,), (batch_bases = -5,))
            @test_throws ArgumentError chopper(path, outpath; kw...)
        end
        @test !isfile(outpath)
    end
end

@testset "Command line" begin
    status, out, _ = run_cli(["--help"])
    @test status == 0 && occursin("--trim-approach", String(out))
    status, out, _ = run_cli(["-V"])
    @test status == 0 && startswith(String(out), "Chopper.jl")
    status, _, err = run_cli(["--bogus"])
    @test status == 2 && occursin("unexpected argument '--bogus'", err)
    status, _, err = run_cli(["-q", "abc"])
    @test status == 2 && occursin("invalid value 'abc'", err)
    status, _, err = run_cli(["--mingc", "1.5"])
    @test status == 2 && occursin("Expected a value from 0 to 1", err)
    status, _, err = run_cli(["-q", "1", "-q", "2"])
    @test status == 2 && occursin("cannot be used multiple times", err)
    status, _, err = run_cli(["--trim-approach", "nope"])
    @test status == 2 && occursin("possible values", err)
    status, _, err = run_cli(["-i", joinpath(DATA, "missing.fastq")])
    @test status == 1 && occursin("is invalid", err)
    status, _, err = run_cli(["-c", joinpath(DATA, "missing.fa")]; input = IOBuffer())
    @test status == 1 && occursin("contaminants is invalid", err)
    # Option syntaxes: --name=value, -qVALUE, -q=VALUE, stdin input.
    fq = read(joinpath(DATA, "test.fastq"))
    expected = read(joinpath(REF, "q10_l500.fastq"))
    for args in (["--quality=10", "--minlength=500"], ["-q10", "-l500"], ["-q=10", "-l", "500"])
        status, out, err = run_cli(args; input = IOBuffer(fq))
        @test status == 0 && out == expected && err == "Kept 205 reads out of 250 reads\n"
    end
    _, _, err = run_cli(["-q", "10", "--summary"]; input = IOBuffer(fq))
    @test occursin("failed quality", err)
    mktempdir() do dir
        path = joinpath(dir, "out.fastq.gz")
        status, out, _ = run_cli(["-q", "10", "-l", "500", "-o", path]; input = IOBuffer(fq))
        @test status == 0 && isempty(out)
        @test transcode(GzipDecompressor, read(path)) == expected
    end
end
