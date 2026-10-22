# Filtering and trimming of a whole FASTQ stream.
#
# Records are read in batches on the calling task. Each batch is processed by
# a separate task (filters, trimming, formatting of the output records into a
# buffer) and the buffers are written in input order, so the output does not
# depend on the number of threads. (chopper writes records in the order in
# which its threads finish them.)

"""
    ChopperStats

Counts returned by [`chopper`](@ref):

- `reads`, `bases`: reads and bases read;
- `segments`, `kept_bases`: records and bases written; `segments` is the
  number that chopper reports as "Kept N reads";
- `kept_reads`: reads with at least one written record;
- `passed`: reads that passed all filters (with `inverse`, these are the
  discarded ones);
- `failed_length`, `failed_quality`, `failed_gc`, `contaminants`: reads by
  the first filter they failed;
- `trimmed_away`: reads selected by the filters but with no segment left
  after trimming (or none long enough);
- `reads_with_mods`: reads that carried an `MM` tag (only counted with
  `update_mods`).
"""
Base.@kwdef mutable struct ChopperStats
    reads::Int = 0
    bases::Int = 0
    kept_reads::Int = 0
    segments::Int = 0
    kept_bases::Int = 0
    passed::Int = 0
    failed_length::Int = 0
    failed_quality::Int = 0
    failed_gc::Int = 0
    contaminants::Int = 0
    trimmed_away::Int = 0
    reads_with_mods::Int = 0
end

function Base.merge!(a::ChopperStats, b::ChopperStats)
    for f in fieldnames(ChopperStats)
        setfield!(a, f, getfield(a, f) + getfield(b, f))
    end
    return a
end

function reset!(s::ChopperStats)
    for f in fieldnames(ChopperStats)
        setfield!(s, f, 0)
    end
    return s
end

Base.:(==)(a::ChopperStats, b::ChopperStats) =
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(ChopperStats))

# 1234567 -> "1,234,567"
thousands(n::Integer) = replace(string(n), r"(?<=\d)(?=(\d{3})+$)" => ",")

function Base.show(io::IO, ::MIME"text/plain", s::ChopperStats)
    println(io, "ChopperStats: kept ", thousands(s.segments), " reads out of ", thousands(s.reads), " reads")
    pct(x, n) = n == 0 ? "" : string(" (", round(100 * x / n; digits = 2), "%)")
    rows = ("bases read" => s.bases, "bases written" => s.kept_bases,
            "reads with a written record" => s.kept_reads, "passed all filters" => s.passed,
            "failed length" => s.failed_length, "failed quality" => s.failed_quality,
            "failed GC content" => s.failed_gc, "contaminants" => s.contaminants,
            "nothing left after trimming" => s.trimmed_away, "reads with MM tags" => s.reads_with_mods)
    for (i, (label, v)) in enumerate(rows)
        print(io, "  ", rpad(label * ":", 30), thousands(v))
        label in ("bases written",) && print(io, pct(v, s.bases))
        label in ("passed all filters", "failed length", "failed quality", "failed GC content",
                  "contaminants") && print(io, pct(v, s.reads))
        i < length(rows) && println(io)
    end
end

# The name of a read, as chopper reports it: the header up to the first
# space or tab.
function read_name(header::AbstractVector{UInt8})
    stop = findfirst(b -> b == UInt8(' ') || b == UInt8('\t'), header)
    return String(view(header, 1:(stop === nothing ? length(header) : stop - 1)))
end

# One batch of records and everything a task needs to process it.
mutable struct BatchSlot
    records::Vector{FASTQ.Record}
    n::Int
    buffer::IOBuffer
    stats::ChopperStats
    segments::Vector{UnitRange{Int}}
end

BatchSlot(cfg::ChopperConfig) = BatchSlot(FASTQ.Record[], 0, IOBuffer(), ChopperStats(), UnitRange{Int}[])

# Reads up to `max_records` records or about `max_bases` bases into `slot`;
# returns the number read (0 at the end of the input).
function fill_batch!(slot::BatchSlot, reader::FASTQ.Reader, max_records::Int, max_bases::Int)
    records = slot.records
    n = 0
    bases = 0
    while n < max_records && bases < max_bases
        n += 1
        n > length(records) && push!(records, FASTQ.Record())
        try
            read!(reader, records[n])
        catch e
            e isa EOFError || rethrow()
            n -= 1
            break
        end
        bases += seqsize(records[n])
    end
    slot.n = n
    return n
end

function process_batch!(slot::BatchSlot, cfg::ChopperConfig)
    reset!(slot.stats)
    for i in 1:slot.n
        process_record!(slot, cfg, slot.records[i])
    end
    return slot
end

function process_record!(slot::BatchSlot, cfg::ChopperConfig, rec::FASTQ.Record)
    stats = slot.stats
    buf = slot.buffer
    segs = slot.segments
    seq = codeunits(sequence(StringView, rec))
    qual = codeunits(quality(StringView, rec))
    header = codeunits(description(rec))
    stats.reads += 1
    stats.bases += length(seq)
    status = read_segments!(segs, cfg, seq, qual)
    if status == PASSED
        stats.passed += 1
    elseif status == FAILED_LENGTH
        stats.failed_length += 1
    elseif status == FAILED_QUALITY
        stats.failed_quality += 1
    elseif status == FAILED_GC
        stats.failed_gc += 1
    else
        stats.contaminants += 1
    end
    if isempty(segs)
        (status == PASSED) != cfg.inverse && (stats.trimmed_away += 1)
        return
    end
    stats.kept_reads += 1
    has_mm = false
    nseg = length(segs)
    for (k, r) in enumerate(segs)
        try
            has_mm |= write_header(buf, header, seq, first(r) - 1, last(r), nseg, k, cfg.update_mods)
        catch e
            e isa ModTagError || rethrow()
            throw(ModTagError(e.kind, e.msg, read_name(header)))
        end
        write(buf, UInt8('\n'))
        write(buf, view(seq, r))
        write(buf, UInt8('\n'), UInt8('+'), UInt8('\n'))
        write(buf, view(qual, r))
        write(buf, UInt8('\n'))
        stats.segments += 1
        stats.kept_bases += length(r)
    end
    has_mm && (stats.reads_with_mods += 1)
    return
end

function flush_slot!(out::IO, total::ChopperStats, slot::BatchSlot)
    write(out, take!(slot.buffer))
    merge!(total, slot.stats)
    return
end

"""
    chopper(input, output = stdout; kwargs...) -> ChopperStats

Filter and trim the reads of the FASTQ `input` and write the selected reads
to `output`, as the command-line tool chopper does. `input` is a path or an
`IO` (compressed input is recognised automatically, see [`open_fastq`](@ref));
`output` is a path (gzip-compressed if it ends in `.gz`) or an `IO` (not
closed). Returns the counts as a [`ChopperStats`](@ref).

A path output is written to a temporary file in the same directory and
renamed over `output` only when the run succeeds: after an error the file
is missing or unchanged, never partial. Bytes already written to an `IO`
(such as `stdout`) cannot be taken back: after an error that output is
incomplete.

Keywords:

- the filters of [`ReadFilter`](@ref): `min_quality`, `max_quality`,
  `min_length`, `max_length`, `min_gc`, `max_gc`, `contaminants`,
  `contaminant_fraction`; or `filter`, a `ReadFilter`, instead of them;
- `trim = NoTrim()`: the trimming strategy, [`FixedCrop`](@ref),
  [`TrimByQuality`](@ref), [`BestReadSegment`](@ref) or
  [`SplitByLowQuality`](@ref);
- `inverse = false`: write the reads that fail the filters instead;
- `update_mods = false`: rewrite the base modification tags `MM`, `ML` and
  `MN` for the trimmed reads (see [`rewrite_mod_tags`](@ref)); a read whose
  tags cannot be rewritten raises a [`ModTagError`](@ref);
- `threads = Threads.nthreads()`: tasks processing batches in parallel. The
  output is identical for any number of threads;
- `batch_bases = 4_000_000`, `batch_reads = 10_000`: size of a batch.

Filters judge the whole read; trimming is applied afterwards to the reads
they select, and only the minimum length is checked again on each segment.
In particular the mean quality is not recomputed for trimmed segments.
"""
function chopper(input, output = stdout; filter::Union{Nothing,ReadFilter} = nothing,
                 trim::TrimStrategy = NoTrim(), inverse::Bool = false, update_mods::Bool = false,
                 threads::Integer = Threads.nthreads(), batch_bases::Integer = 4_000_000,
                 batch_reads::Integer = 10_000, filter_options...)
    if filter === nothing
        filter = ReadFilter(; filter_options...)
    elseif !isempty(filter_options)
        throw(ArgumentError("pass either `filter` or the filter keywords, not both"))
    end
    threads >= 1 || throw(ArgumentError("threads must be at least 1"))
    batch_reads >= 1 || throw(ArgumentError("batch_reads must be at least 1"))
    batch_bases >= 1 || throw(ArgumentError("batch_bases must be at least 1"))
    # Opening the output truncates it: refuse the input file (also through a
    # different spelling, link or letter case)
    input isa IO || output isa IO || !ispath(String(output)) ||
        !samefile(String(input), String(output)) ||
        throw(ArgumentError("the output file is the input file: $(output)"))
    cfg = ChopperConfig(filter, trim, inverse, update_mods)
    input isa IO || isfile(String(input)) || throw(ArgumentError("input file $input is invalid"))
    target = output isa IO ? nothing : OutputTarget(String(output))
    instream = ahead = out = nothing
    success = false
    try
        instream = input isa IO ? open_fastq(input) : open_fastq(String(input))
        # Compressed input is decompressed by a separate task when there are
        # threads to spare.
        ahead = instream isa NoopStream || threads == 1 || Threads.nthreads() == 1 ? nothing : prefetch(instream)
        out = target === nothing ? output : open_output(target)
        reader = FASTQ.Reader(ahead === nothing ? instream : NoopStream(ahead); copy = false)
        stats = run_batches(reader, out, cfg, Int(threads), Int(batch_reads), Int(batch_bases))
        flush(out)
        if target !== nothing
            o, out = out, nothing
            close(o)
            publish(target)
        end
        success = true
        return stats
    finally
        if ahead !== nothing
            # After a complete run `eof` has already reported a producer error,
            # and a cancelled producer ends without one; after a failure the
            # first error is the one reported.
            close(ahead)
            try
                wait(ahead.task)
            catch
                success && rethrow()
            end
        end
        input isa IO || instream === nothing || close(instream)
        target === nothing || out === nothing || close(out)
        target === nothing || success || discard(target)
    end
end

function run_batches(reader::FASTQ.Reader, out::IO, cfg::ChopperConfig, threads::Int,
                     batch_reads::Int, batch_bases::Int)
    total = ChopperStats()
    if threads == 1 || Threads.nthreads() == 1
        slot = BatchSlot(cfg)
        while fill_batch!(slot, reader, batch_reads, batch_bases) > 0
            process_batch!(slot, cfg)
            flush_slot!(out, total, slot)
        end
        return total
    end
    # Up to 2 × threads batches in flight; the oldest is written first.
    free = [BatchSlot(cfg) for _ in 1:(2 * threads)]
    pending = Tuple{Task,BatchSlot}[]
    try
        while true
            if isempty(free)
                task, slot = popfirst!(pending)
                fetch(task)
                flush_slot!(out, total, slot)
                push!(free, slot)
            end
            slot = pop!(free)
            if fill_batch!(slot, reader, batch_reads, batch_bases) == 0
                push!(free, slot)
                break
            end
            push!(pending, (Threads.@spawn(process_batch!($slot, $cfg)), slot))
        end
        while !isempty(pending)
            task, slot = popfirst!(pending)
            fetch(task)
            flush_slot!(out, total, slot)
        end
    catch e
        # Let the other tasks finish before reporting the first error.
        for (task, _) in pending
            try
                wait(task)
            catch
            end
        end
        e isa TaskFailedException ? throw(e.task.exception) : rethrow()
    end
    return total
end
