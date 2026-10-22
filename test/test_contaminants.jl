# Contamination screen (mapping by Minimap2.jl, decision by Chopper).

using Random

const COMP = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')
revcomp(s) = String(reverse([COMP[c] for c in s]))

# Read sampled from `genome` with substitution, insertion and deletion rate
# `err / 3` each.
function simulate_read(rng, genome, start, len, err)
    out = Char[]
    i = start
    while i < start + len && i <= length(genome)
        r = rand(rng)
        if r < err / 3
            push!(out, rand(rng, setdiff("ACGT", genome[i])))
            i += 1
        elseif r < 2err / 3
            push!(out, rand(rng, "ACGT"))
        elseif r < err
            i += 1
        else
            push!(out, genome[i])
            i += 1
        end
    end
    return String(out)
end

@testset "Contaminant screen" begin
    rng = MersenneTwister(2026)
    genome = String(rand(rng, collect("ACGT"), 30_000))
    other = String(rand(rng, collect("ACGT"), 30_000))
    idx = ContaminantIndex(["g"], [genome])
    @test occursin("1 sequences", sprint(show, idx))
    @testset "simulated reads, error $err" for err in (0.0, 0.03, 0.06)
        for trial in 1:20
            len = rand(rng, 500:4000)
            start = rand(rng, 1:(length(genome) - len))
            read = simulate_read(rng, genome, start, len, err)
            rev = isodd(trial)
            r = rev ? revcomp(read) : read
            @test is_contaminant(idx, r)
        end
    end
    # Unrelated reads, chimeras with half the read from the contaminant, and
    # reads with 25% errors are not contaminants.
    for trial in 1:10
        @test !is_contaminant(idx, other[(1000trial):(1000trial + 1500)])
        half = genome[(2000trial):(2000trial + 999)] * other[(2000trial):(2000trial + 999)]
        @test !is_contaminant(idx, half)
        @test best_match_length(idx, half) >= 900
        noisy = simulate_read(rng, genome, 1000trial, 2000, 0.25)
        @test !is_contaminant(idx, noisy)
    end
    @test !is_contaminant(idx, "")
    @test best_match_length(idx, "ACGT") == 0
    @test is_contaminant(idx, genome[101:300]; min_fraction = 1.0)
    # Reads that run past the end of the contaminant.
    @test best_match_length(idx, other[1:300] * genome[1:1000]) >= 990
    @test best_match_length(idx, genome[(end - 999):end] * other[1:300]) >= 990
    # Several contaminant sequences, ambiguous bases.
    idx2 = ContaminantIndex(["o", "g"], [other, genome[1:10_000] * "NNNNN" * genome[10_001:end]])
    @test is_contaminant(idx2, genome[9000:11_000])
    @test is_contaminant(idx2, revcomp(other[5000:6000]))
    # Minimizers above the occurrence threshold are not ignored when the read
    # has no others: minimap2 2.30 maps this read (match length 2,000).
    rep = repeat("ACGTTGCAAGGCTTACCGATAGGCTA", 600)
    idx_rep = ContaminantIndex(["r"], [rep])
    @test idx_rep.options.mid_occ == 500
    @test best_match_length(idx_rep, rep[1:2000]) == 2000
    # From a FASTA file, as the -c option.
    mktempdir() do dir
        path = joinpath(dir, "c.fa")
        write(path, ">g desc\n", genome[1:15_000], "\n", genome[15_001:end], "\n")
        idx3 = ContaminantIndex(path)
        @test idx3.names == ["g"] && idx3.index.lengths == [30_000]
        f = ReadFilter(contaminants = path)
        @test read_status(f, genome[5000:7000], "I"^2001) == Chopper.CONTAMINANT
        @test_throws ArgumentError ContaminantIndex(joinpath(dir, "none.fa"))
        write(joinpath(dir, "empty.fa"), "")
        @test_throws ArgumentError ContaminantIndex(joinpath(dir, "empty.fa"))
    end
    @test_throws ArgumentError ContaminantIndex(["a"], ["ACGT"]; k = 40)
    @test_throws ArgumentError ContaminantIndex(String[], String[])
end

@testset "chopper's contamination tests" begin
    idx = ContaminantIndex(joinpath(DATA, "random_contam.fa"))
    first_seq(path) = readlines(path)[2]
    @test is_contaminant(idx, first_seq(joinpath(DATA, "test.fastq")))
    @test !is_contaminant(idx, first_seq(joinpath(DATA, "other-test.fastq")))
end

@testset "Match lengths identical to minimap2" begin
    # minimap2 2.30 -x lr:hq -c (tools/make_contaminant_reference.jl)
    expected = Dict(first(l) => parse(Int, last(l)) for l in split.(readlines(joinpath(REF, "contaminant_match_lengths.tsv"))[2:end], '	'))
    idx = ContaminantIndex(joinpath(DATA, "random_contam.fa"))
    lines = readlines(joinpath(DATA, "test.fastq"))
    names = [split(lines[i][2:end])[1] for i in 1:4:length(lines)]
    ours = Vector{Int}(undef, length(names))
    Threads.@threads for i in eachindex(names)
        ours[i] = best_match_length(idx, lines[4i - 2])
    end
    @test length(expected) == 250
    @test ours == [expected[n] for n in names]
end
