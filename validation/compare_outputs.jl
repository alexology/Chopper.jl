# Compares Chopper.jl with outputs of the original chopper on a FASTQ file:
#
#     julia -t 4 --project=. validation/compare_outputs.jl reads.fastq.gz reference_dir
#
# `reference_dir` holds the outputs of `chopper -t 1` for the settings below
# (q20_t1.fq, split_t1.fq, best_t1.fq, gcinv_t1.fq; see validation/README.md).
# For every setting the script runs Chopper.jl with 1 thread and with all
# threads, reports the time and whether the output is byte-identical.

using Chopper

const SETTINGS = [
    ("q20", (; min_quality = 20)),
    ("split", (; trim = SplitByLowQuality(15), min_length = 200)),
    ("best", (; min_quality = 15, trim = BestReadSegment(18), min_length = 500)),
    ("gcinv", (; min_gc = 0.35, max_gc = 0.45, inverse = true, min_quality = 12, max_length = 3000)),
]

function main(input, refdir)
    # Warm-up (compilation).
    chopper(IOBuffer("@r\nACGT\n+\nIIII\n"), IOBuffer(); threads = Threads.nthreads())
    for (name, kw) in SETTINGS
        expected = read(joinpath(refdir, "$(name)_t1.fq"))
        for threads in unique([1, Threads.nthreads()])
            out = IOBuffer()
            t = @elapsed stats = chopper(input, out; threads, kw...)
            same = take!(out) == expected
            println(rpad(name, 6), " threads = ", threads, ": ", round(t; digits = 2), " s, ",
                    stats.segments, " records, identical to chopper: ", same)
        end
    end
end

main(ARGS[1], ARGS[2])
