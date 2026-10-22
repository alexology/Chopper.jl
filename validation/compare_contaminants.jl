# Compares the contamination screen with minimap2 and chopper on real reads:
#
#     julia -t 4 --project=. validation/compare_contaminants.jl reads.fq contaminants.fa minimap2.paf chopper_output.fq
#
# minimap2.paf: `minimap2 -x lr:hq -c contaminants.fa reads.fq` (minimap2
# 2.30, the version inside chopper 0.14.1); chopper_output.fq: `chopper -t 1
# -c contaminants.fa -i reads.fq`. Reports the agreement of the match lengths
# (largest over the alignments of a read), of the classification at 90%, and
# whether the output is identical to chopper's.

using Chopper

function main(reads, contam, paf, chopper_out)
    idx = ContaminantIndex(contam)
    lines = readlines(reads)
    names = [split(lines[i][2:end])[1] for i in 1:4:length(lines)]
    seqs = lines[2:4:end]
    mm2 = Dict{String,Int}()
    for l in eachline(paf)
        f = split(l, '\t')
        mm2[f[1]] = max(get(mm2, f[1], 0), parse(Int, f[10]))
    end
    ours = similar(seqs, Int)
    t = @elapsed Threads.@threads for i in eachindex(seqs)
        ours[i] = best_match_length(idx, seqs[i])
    end
    theirs = [get(mm2, n, 0) for n in names]
    d = ours .- theirs
    pct(x) = string(round(100x / length(d); digits = 2), "%")
    println("reads: ", length(d), "; aligned by minimap2: ", count(>(0), theirs), ", by Chopper.jl: ", count(>(0), ours))
    println("match length identical: ", pct(count(==(0), d)), "; within 5 bases: ", pct(count(x -> abs(x) <= 5, d)),
            "; within 1% of the read: ", pct(count(i -> abs(d[i]) <= 0.01 * length(seqs[i]), eachindex(d))))
    near = [i for i in eachindex(d) if abs(theirs[i] / length(seqs[i]) - 0.9) < 0.02]
    println("reads with minimap2 fraction within 0.02 of 0.9: ", length(near), ", match length identical: ",
            count(i -> d[i] == 0, near))
    c_ours = [ours[i] >= 0.9 * length(seqs[i]) for i in eachindex(d)]
    c_theirs = [theirs[i] >= 0.9 * length(seqs[i]) for i in eachindex(d)]
    println("contaminants: minimap2 ", count(c_theirs), ", Chopper.jl ", count(c_ours),
            "; classified differently: ", count(c_ours .!= c_theirs))
    for i in findall(c_ours .!= c_theirs)
        println("   ", names[i], " length ", length(seqs[i]), ": minimap2 ", theirs[i], " (",
                round(theirs[i] / length(seqs[i]); digits = 4), "), Chopper.jl ", ours[i], " (",
                round(ours[i] / length(seqs[i]); digits = 4), ")")
    end
    println("time of best_match_length on all reads: ", round(t; digits = 1), " s (", Threads.nthreads(), " threads)")
    out = IOBuffer()
    t = @elapsed chopper(reads, out; filter = ReadFilter(contaminants = idx))
    println("chopper(): ", round(t; digits = 1), " s; output identical to chopper: ", take!(out) == read(chopper_out))
end

main(ARGS...)
