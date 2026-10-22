# Writes test/reference/contaminant_match_lengths.tsv: for every read of
# test/data/test.fastq, the largest match length (PAF column 10) of minimap2
# 2.30 with the preset used by chopper (`-x lr:hq -c`) against
# test/data/random_contam.fa; 0 for reads that do not map.
#
#     julia tools/make_contaminant_reference.jl path/to/minimap2
#
# minimap2 must be version 2.30 (the one inside chopper 0.14.1), e.g. the
# build in ../Minimap2.jl/reference/minimap2.

minimap2 = only(ARGS)
version = readchomp(`$minimap2 --version`)
startswith(version, "2.30") || error("minimap2 2.30 is needed, found $version")
data = joinpath(@__DIR__, "..", "test", "data")
best = Dict{String,Int}()
for line in eachline(pipeline(`$minimap2 -x lr:hq -c -t 1 $(joinpath(data, "random_contam.fa")) $(joinpath(data, "test.fastq"))`;
                              stderr = devnull))
    f = split(line, '\t')
    best[f[1]] = max(get(best, f[1], 0), parse(Int, f[10]))
end
names = [split(l[2:end])[1] for (i, l) in enumerate(readlines(joinpath(data, "test.fastq"))) if i % 4 == 1]
open(joinpath(@__DIR__, "..", "test", "reference", "contaminant_match_lengths.tsv"), "w") do io
    println(io, "read\tmatch_length")
    for n in names
        println(io, n, '\t', get(best, n, 0))
    end
end
println("reads: ", length(names), ", mapped: ", length(best))
