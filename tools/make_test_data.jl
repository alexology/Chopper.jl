# Generates the synthetic test inputs in test/data/ (the other files there are
# copied from chopper's test-data/):
#
#     julia --project=. tools/make_test_data.jl
#
# - mods.fastq: reads of test.fastq, cut to 150-1,500 bases, with base
#   modification tags as written by `samtools fastq -T MM,ML,MN` (several
#   modification lists, two-code lists, N lists, lower-case tag names, a
#   CO:Z: tag containing spaces before MM, MM without ML);
# - headers.fastq: header edge cases (descriptions, tabs, trailing spaces);
# - bad_mods.fastq: a read whose MM tag refers to more C bases than it has.

using Random

const DATA = joinpath(@__DIR__, "..", "test", "data")

function read_fastq(path)
    lines = readlines(path)
    return [(lines[i][2:end], lines[i + 1], lines[i + 3]) for i in 1:4:length(lines)]
end

# An MM list for canonical base `base` with calls on a random subset of its
# occurrences; returns the list body (skip counts) and the number of calls.
function mm_list(rng, seq, base)
    occ = [i for (i, c) in enumerate(seq) if base == 'N' || uppercase(c) == base]
    calls = sort!(shuffle(rng, occ)[1:min(length(occ), rand(rng, 0:12))])
    skips = Int[]
    last = 0
    for c in calls
        k = findfirst(==(c), occ)
        push!(skips, k - last - 1)
        last = k
    end
    return skips
end

function mods_reads(rng, reads)
    out = String[]
    for (n, (header, seq, qual)) in enumerate(reads[1:60])
        len = min(length(seq), rand(rng, 150:1500))
        seq, qual = seq[1:len], qual[1:len]
        name = split(header)[1]
        lists = String[]
        nvalues = 0
        specs = n % 5 == 0 ? [("C", "+mh?", 2), ("A", "+a.", 1)] :
                n % 5 == 1 ? [("C", "+m?", 1)] :
                n % 5 == 2 ? [("C", "+m?", 1), ("G", "-m", 1)] :
                n % 5 == 3 ? [("N", "+n?", 1)] : [("C", "+76792?", 1)]
        for (base, codes, stride) in specs
            skips = mm_list(rng, seq, only(base))
            push!(lists, base * codes * join(("," * string(s) for s in skips)) * ";")
            nvalues += stride * length(skips)
        end
        mm = join(lists)
        ml = join(("," * string(rand(rng, 0:255)) for _ in 1:nvalues))
        mmtag = n % 7 == 0 ? "Mm:Z:" : "MM:Z:"
        mltag = n % 7 == 0 ? "Ml:B:C" : "ML:B:C"
        tags = ["qs:f:$(round(10 + 20rand(rng); digits = 3))"]
        n % 6 == 0 && push!(tags, "CO:Z:some comment with spaces")
        push!(tags, mmtag * mm)
        n % 9 == 0 || push!(tags, mltag * ml)        # MM without ML is valid
        push!(tags, "MN:i:$len")
        push!(out, "@" * name * "\t" * join(tags, '\t'), seq, "+", qual)
    end
    return out
end

function header_reads(rng, reads)
    headers = ["plain", "with_desc some free text", "tabs\tqs:f:12.5\tRG:Z:run1",
               "trailing_space   ", "trailing_tab\t", "desc_and_tags a b\tqs:f:3\tCO:Z:x y",
               "space_in_tag\tCO:Z:a b\tRG:Z:z", "multi  spaces  inside"]
    out = String[]
    for (i, h) in enumerate(headers)
        _, seq, qual = reads[100 + i]
        # Quality dips so that the read is split by split-by-low-quality.
        q = collect(qual[1:400])
        q[120:140] .= '#'
        q[260:262] .= '#'
        push!(out, "@" * h, seq[1:400], "+", String(q))
    end
    return out
end

function main()
    rng = MersenneTwister(20261021)
    reads = read_fastq(joinpath(DATA, "test.fastq"))
    write(joinpath(DATA, "mods.fastq"), join(mods_reads(rng, reads), '\n') * '\n')
    write(joinpath(DATA, "headers.fastq"), join(header_reads(rng, reads), '\n') * '\n')
    _, seq, qual = reads[3]
    write(joinpath(DATA, "bad_mods.fastq"),
          "@ok\tMM:Z:C+m?,0;\tML:B:C,7\n$(seq[1:300])\n+\n$(qual[1:300])\n" *
          "@broken\tMM:Z:C+m?,100000;\tML:B:C,200\n$(seq[1:300])\n+\n$(qual[1:300])\n")
    println("written mods.fastq, headers.fastq, bad_mods.fastq in ", normpath(DATA))
end

main()
