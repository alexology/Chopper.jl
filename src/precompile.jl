# Precompilation workload: shortens the start-up of `bin/chopper.jl`.

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    fastq = "@r1\tqs:f:20\nACGTACGTTTGACCAGTAGGACCATGACTAGGG\n+\n" * "I"^10 * "#"^3 * "I"^20 * "\n" *
            "@r2 desc\nACGTACGTACGTACGT\n+\n" * "5"^16 * "\n"
    contam_seq = "ACGTACGTTTGACCAGTAGGACCATGACTAGGGTTACAGGATTACGATCAGGACTTAGCA"
    @compile_workload begin
        index = ContaminantIndex(["c"], [contam_seq])
        for trim in (NoTrim(), FixedCrop(2, 2), TrimByQuality(10), BestReadSegment(10),
                     SplitByLowQuality(10))
            chopper(IOBuffer(fastq), IOBuffer(); min_quality = 5, trim, threads = 1)
        end
        chopper(IOBuffer(fastq), IOBuffer(); contaminants = index, min_gc = 0.1, threads = 1)
        main(["-q", "5"]; input = IOBuffer(fastq), output = IOBuffer(), err = IOBuffer())
    end
end
