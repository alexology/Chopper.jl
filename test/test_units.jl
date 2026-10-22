# Unit tests, most of them ported from chopper's Rust tests.

# The six reads of chopper's trimmers.rs tests.
const TRIM_READS = [
    ("AAAAAAAAAAAAAAATTTAA", "&#3-G27C:(@G7B55+C4I"),
    ("TTTTTTTTTTTTTTTTTTTT", "77%'24)FAF9@=94'%054"),
    ("AAAAAAAAAAAAAAATTTTA", "'8\$-BF2!C;+59->H@91#"),
    ("AAAAAAAAAAAAAAAAAAAA", "%,42\$CH*#0+0C6=0,*6/"),
    ("AAAAAAAAAAAAAAAAAAAT", "-------------------J"),
    ("TAAAAAAAAAAAAAAAAAAA", "I-------------------"),
]

@testset "Phred probabilities and mean quality" begin
    @test phred_to_probability(20) == 0.01
    @test phred_to_probability(30) == 0.001
    @test phred_to_probability(15) == 0.03162277660168379
    @test phred_to_probability(25) == 0.0031622776601683794
    @test phred_to_probability(0) == 1.0
    q(x...) = UInt8[v + 33 for v in x]
    @test mean_quality(q(10)) == 10.0
    @test mean_quality(q(10, 11, 12)) ≈ 10.923583702678473 atol = 1e-5
    @test mean_quality(q(10, 11, 12, 20, 30, 40, 50)) ≈ 14.408827647036087 atol = 1e-5
    # First 100 qualities of the first read of test.fastq (chopper's test_quals).
    first_quals = [17, 19, 11, 5, 3, 19, 22, 24, 20, 22, 30, 31, 32, 20, 21, 30, 28, 10, 13, 12, 18,
                   18, 18, 19, 24, 25, 35, 33, 34, 35, 34, 27, 29, 25, 21, 18, 19, 12, 14, 15, 24, 26,
                   24, 7, 12, 17, 17, 19, 17, 8, 14, 15, 13, 15, 9, 3, 4, 23, 23, 29, 23, 10, 29, 30,
                   31, 27, 25, 14, 2, 13, 19, 14, 13, 13, 3, 2, 10, 17, 19, 25, 27, 20, 19, 11, 5, 7,
                   8, 8, 5, 2, 10, 12, 16, 18, 16, 14, 12, 15, 2, 3]
    qual1 = readlines(joinpath(DATA, "test.fastq"))[4]
    @test Int.(codeunits(qual1)[1:100]) .- 33 == first_quals
    # chopper's 200-value test of ave_qual (first 200 qualities of that read).
    @test mean_quality(codeunits(qual1)[1:200]) ≈ 10.017407548271677 atol = 1e-5
    @test isnan(mean_quality(UInt8[]))
    @test gc_content(codeunits("GCgcAT")) ≈ 4 / 6
    @test gc_content(codeunits("NNGC")) == 0.5           # N counts in the denominator
    @test isnan(gc_content(UInt8[]))
end

@testset "Trimming strategies" begin
    quals = [codeunits(q) for (_, q) in TRIM_READS]
    @testset "BestReadSegment" begin
        cases = [(20, [(10, 16)]), (7, [(0, 20)]), (15, [(11, 19)]), (40, []), (40, [(19, 20)]), (40, [(0, 1)])]
        for ((cutoff, expected), q) in zip(cases, quals)
            @test trim_ranges(BestReadSegment(cutoff), q) == ranges(expected)
        end
        @test BestReadSegment(7).probability == 0.19952623149688797
    end
    @testset "TrimByQuality" begin
        cases = [(20, [(4, 20)]), (7, [(0, 20)]), (15, [(1, 19)]), (40, []), (40, [(19, 20)]), (40, [(0, 1)])]
        for ((cutoff, expected), q) in zip(cases, quals)
            @test trim_ranges(TrimByQuality(cutoff), q) == ranges(expected)
        end
    end
    @testset "FixedCrop" begin
        cases = [(5, 3, [(5, 17)]), (30, 4, []), (1, 1, [(1, 19)]), (15, 30, []), (19, 0, [(19, 20)]), (0, 19, [(0, 1)])]
        for ((head, tail, expected), q) in zip(cases, quals)
            @test trim_ranges(FixedCrop(head, tail), q) == ranges(expected)
        end
        @test_throws ArgumentError FixedCrop(0, 0)
        @test_throws ArgumentError FixedCrop(-1, 3)
        @test FixedCrop(head = 2) == FixedCrop(2, 0)
    end
    @testset "SplitByLowQuality" begin
        # (cutoff, min_length, expected); the length filter is applied after
        # the strategy, as `read_segments` does.
        cases = [(20, 3, [(6, 9), (10, 16)]), (7, 3, [(4, 15), (17, 20)]), (15, 3, [(4, 7), (14, 19)]),
                 (40, 3, []), (40, 1, [(19, 20)]), (40, 1, [(0, 1)])]
        for ((cutoff, minlen, expected), q) in zip(cases, quals)
            @test Base.filter(r -> length(r) >= minlen, trim_ranges(SplitByLowQuality(cutoff), q)) == ranges(expected)
        end
        # chopper's extra cases on reads 1 and 2.
        @test Base.filter(r -> length(r) >= 5, trim_ranges(SplitByLowQuality(10), quals[1])) == ranges([(2, 9), (10, 20)])
        @test Base.filter(r -> length(r) >= 3, trim_ranges(SplitByLowQuality(10), quals[2])) == ranges([(7, 15), (17, 20)])
        @test trim_ranges(SplitByLowQuality(10), "III###III###III###III") == ranges([(0, 3), (6, 9), (12, 15), (18, 21)])
        @test trim_ranges(SplitByLowQuality(10), "") == UnitRange{Int}[]
        @test trim_ranges(SplitByLowQuality(10), "####") == UnitRange{Int}[]
        @test trim_ranges(SplitByLowQuality(10), "IIII") == [1:4]
        q = "III#IIII###III"
        @test trim_ranges(SplitByLowQuality(10; window = 1), q) == ranges([(0, 3), (4, 8), (11, 14)])
        @test trim_ranges(SplitByLowQuality(10; window = 3), q) == ranges([(0, 8), (11, 14)])
        @test trim_ranges(SplitByLowQuality(10; window = 4), q) == ranges([(0, 14)])
        @test SplitByLowQuality(10; window = 0).window == 1
    end
    @test trim_ranges(NoTrim(), "IIII") == [1:4]
    @test_throws ArgumentError TrimByQuality(256)
    out = UnitRange{Int}[1:2, 3:4]
    @test trim_ranges!(out, TrimByQuality(10), codeunits("##II#")) === out && out == [3:4]
end

# A user-defined strategy, like chopper's MockTrimmer.
struct MockTrim <: TrimStrategy
    segments::Vector{UnitRange{Int}}
end
Chopper.trim_ranges!(out::Vector{UnitRange{Int}}, m::MockTrim, qual::AbstractVector{UInt8}) =
    append!(empty!(out), m.segments)

@testset "Filters and segments (get_valid_segments)" begin
    f = ReadFilter(min_length = 5)
    @test read_segments("ATGC", "IIII"; filter = f) == UnitRange{Int}[]
    @test read_segments("", ""; filter = ReadFilter(min_length = 0)) == [1:0]
    @test read_segments("ATGCGTACGA", "IIIIIIIIII"; filter = f, trim = MockTrim([1:10])) == [1:10]
    @test read_segments("ATGCGTACGA", "IIIIIIIIII"; filter = ReadFilter(min_length = 4),
                        trim = MockTrim([1:3, 4:5, 6:10])) == [6:10]
    @test read_status(ReadFilter(min_quality = 41), "ACGT", "IIII") == Chopper.FAILED_QUALITY
    @test read_status(ReadFilter(max_quality = 39), "ACGT", "IIII") == Chopper.FAILED_QUALITY
    @test read_status(ReadFilter(min_gc = 0.6), "ACGT", "IIII") == Chopper.FAILED_GC
    @test read_status(ReadFilter(max_length = 3), "ACGT", "IIII") == Chopper.FAILED_LENGTH
    @test read_status(ReadFilter(), "ACGT", "!!!!") == Chopper.PASSED   # quality filter inactive
    # Inverse: failing reads are kept, and trimmed segments are length-checked.
    @test read_segments("ACGT", "IIII"; filter = ReadFilter(min_length = 5), inverse = true) == [1:4]
    @test read_segments("ACGTACGT", "IIIIIIII"; inverse = true) == UnitRange{Int}[]
    @test read_segments("ACGTAC", "II#III"; filter = ReadFilter(min_length = 7), inverse = true,
                        trim = SplitByLowQuality(10)) == UnitRange{Int}[]
    @test_throws ArgumentError ReadFilter(min_gc = 1.5)
    @test_throws DimensionMismatch read_segments("ACGT", "III")
end

# Formats one record as `chopper` writes it.
function format_record(header, seq, qual, start, stop, nseg, index; update_mods = false)
    io = IOBuffer()
    s = codeunits(seq)
    Chopper.write_header(io, codeunits(header), s, start, stop, nseg, index, update_mods)
    print(io, '\n', seq[(start + 1):stop], "\n+\n", qual[(start + 1):stop], '\n')
    return String(take!(io))
end

@testset "Record headers (records.rs)" begin
    @test format_record("10-bases", "AAAAAAAAAA", "IIIIIIIIII", 0, 10, 1, 1) == "@10-bases\nAAAAAAAAAA\n+\nIIIIIIIIII\n"
    @test format_record("10-bases", "TTAAAAAATT", "KKIIIIIIKK", 2, 8, 1, 1) == "@10-bases\nAAAAAA\n+\nIIIIII\n"
    @test format_record("10-bases", "TTAAAAAATT", "KKIIIIIIKK", 2, 8, 2, 2) == "@10-bases_segment_2\nAAAAAA\n+\nIIIIII\n"
    @test format_record("10-bases description", "TTAAAAAATT", "KKIIIIIIKK", 2, 8, 2, 2) ==
          "@10-bases_segment_2 description\nAAAAAA\n+\nIIIIII\n"
    @test format_record("10-bases\tMM:Z:C+m?,1;\tML:B:C,200", "TTAAAAAATT", "KKIIIIIIKK", 2, 8, 1, 1) ==
          "@10-bases\tMM:Z:C+m?,1;\tML:B:C,200\nAAAAAA\n+\nIIIIII\n"
    @test format_record("10-bases\tMM:Z:C+m?,1;\tML:B:C,200", "TTAAAAAATT", "KKIIIIIIKK", 2, 8, 2, 2) ==
          "@10-bases_segment_2\tMM:Z:C+m?,1;\tML:B:C,200\nAAAAAA\n+\nIIIIII\n"
    @test format_record("10-bases\tMM:Z:C+m?,1; description", "TTAAAAAATT", "KKIIIIIIKK", 2, 8, 2, 2) ==
          "@10-bases_segment_2\tMM:Z:C+m?,1; description\nAAAAAA\n+\nIIIIII\n"
    @test format_record("read1\tMM:Z:C+m?,1;\tML:B:C,200\tqs:f:20.1", "ACGTCAAAC", "IIIIIIIII", 2, 9, 1, 1;
                        update_mods = true) == "@read1\tMM:Z:C+m?,0;\tML:B:C,200\tqs:f:20.1\nGTCAAAC\n+\nIIIIIII\n"
    @test format_record("read1\tMM:Z:C+m?,1;\tML:B:C,200", "ACGTCAAAC", "IIIIIIIII", 2, 9, 1, 1) ==
          "@read1\tMM:Z:C+m?,1;\tML:B:C,200\nGTCAAAC\n+\nIIIIIII\n"
    @test_throws ModTagError format_record("read1\tMM:Z:C+m?,9;\tML:B:C,200", "ACGTCAAAC", "IIIIIIIII", 0, 9, 1, 1;
                                           update_mods = true)
    # A tag value with a space: chopper's parser splits the header there.
    @test format_record("read1\tCO:Z:some comment\tMM:Z:C+m?,1;\tML:B:C,200", "ACGTCAAAC", "IIIIIIIII", 2, 9, 1, 1;
                        update_mods = true) ==
          "@read1\tCO:Z:some comment\tMM:Z:C+m?,0;\tML:B:C,200\nGTCAAAC\n+\nIIIIIII\n"
    @test format_record("read1 a free text description", "ACGTCAAAC", "IIIIIIIII", 0, 9, 2, 2) ==
          "@read1_segment_2 a free text description\nACGTCAAAC\n+\nIIIIIIIII\n"
    @test format_record("trailing  \t ", "ACGT", "IIII", 0, 4, 1, 1) == "@trailing\nACGT\n+\nIIII\n"
end

@testset "Base modification tags (modtags.rs)" begin
    SEQ = codeunits("ACGTCAAACGGTCAAACGT")   # C at 1, 4, 8, 12, 16 (0-based); G at 2, 9, 10, 17
    n = length(SEQ)
    sub(mm, ml, s, e) = Chopper.subset_mods(mm, ml, SEQ, s, e)
    @test sub("C+m?,1,1;", UInt8[200, 10], 0, n) == ("C+m?,1,1;", UInt8[200, 10])
    @test sub("C+m?,1,1;", UInt8[200, 10], 2, n) == ("C+m?,0,1;", UInt8[200, 10])
    @test sub("C+m?,1,1;", UInt8[200, 10], 0, 10) == ("C+m?,1;", UInt8[200])
    @test sub("C+m?,1,1;", UInt8[200, 10], 10, n) == ("C+m?,0;", UInt8[10])
    @test sub("C+m?,1,1;", UInt8[200, 10], 5, 12) == ("C+m?;", UInt8[])
    @test sub("C+mh?,1,1;", UInt8[200, 5, 10, 3], 10, n) == ("C+mh?,0;", UInt8[10, 3])
    @test sub("C+m?,1,1;G-m,0;", UInt8[200, 10, 250], 10, n) == ("C+m?,0;G-m;", UInt8[10])
    @test sub("C+76792?,1,1;", UInt8[200, 10], 0, n) == ("C+76792?,1,1;", UInt8[200, 10])
    @test sub("C+m?,1,1;", nothing, 10, n) == ("C+m?,0;", UInt8[])
    @test sub("U+m?,1;", UInt8[200], 0, n)[1] == "U+m?,1;"
    err(f) = try f(); nothing catch e; e end
    @test err(() -> sub("C+m?,1,1;", UInt8[200], 0, n)).kind == :count_mismatch
    @test err(() -> sub("C+m?,1,4;", UInt8[200, 10], 0, n)).kind == :beyond_sequence
    for mm in ("C,1,1;", "X+m,1;", "C+,1;", "C+m,x;", "C;", "C+m,-1;", "C+m, 1;", "C+m,1,;")
        @test err(() -> sub(mm, nothing, 0, n)).kind == :malformed_mm
    end
    @test sub("C+m,+1;", nothing, 0, n) == ("C+m,1;", UInt8[])     # Rust accepts a leading '+'
    @test sub("N+n?,3;", UInt8[200], 0, n) == ("N+n?,3;", UInt8[200])
    @test sub("N+n?,3;", UInt8[200], 3, n) == ("N+n?,0;", UInt8[200])
    @test sub("N+n?,5,5;", UInt8[200, 10], 0, n) isa Tuple

    rw(fields, s, e) = Chopper.rewrite_mod_fields(fields, SEQ, s, e)[1]
    @test rw(["MM:Z:C+m?,1,1;", "ML:B:C,200,10", "MN:i:19", "qs:f:20.1"], 10, n) ==
          ["MM:Z:C+m?,0;", "ML:B:C,10", "MN:i:9", "qs:f:20.1"]
    @test rw(["Mm:Z:C+m?,1,1;", "Ml:B:C,200,10"], 10, n) == ["Mm:Z:C+m?,0;", "Ml:B:C,10"]
    @test rw(["MM:Z:C+m?,1,1;", "ML:B:C,200,10"], 5, 12) == ["MM:Z:C+m?;", "ML:B:C"]
    @test rw(["qs:f:20.1", "RG:Z:abc"], 2, 10) == ["qs:f:20.1", "RG:Z:abc"]
    @test err(() -> rw(["MM:Z:C+m?,1,1;", "ML:B:C,200,10", "MN:i:25"], 0, n)).kind == :stale_mn
    @test err(() -> rw(["MM:Z:C+m?,1;", "ML:B:C,300"], 0, n)).kind == :malformed_ml
    @test err(() -> rw(["ML:B:C,200", "MN:i:19"], 2, 10)).kind == :ml_without_mm
    @test err(() -> rw(["MM:Z:C+m?,1;", "MN:i:x"], 0, n)).kind == :malformed_mn
    @test err(() -> rw(["MM:Z:C+m?,1;", "ML:B:C200"], 0, n)).kind == :malformed_ml
    @test rewrite_mod_tags(["MM:Z:C+m?,1,1;", "ML:B:C,200,10"], "ACGTCAAACGGTCAAACGT", 11:19) ==
          ["MM:Z:C+m?,0;", "ML:B:C,10"]
    @test_throws BoundsError rewrite_mod_tags(["MM:Z:C+m?,1;"], "ACGT", 2:9)
    e = ModTagError(:stale_mn, "x", "r1")
    @test sprint(showerror, e) == "read r1: x"
end
