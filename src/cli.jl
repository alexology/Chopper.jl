# Command-line interface with the options of chopper 0.14.1.

const UPSTREAM_VERSION = "0.14.1"

const HELP = """
Filtering and trimming of fastq files. Reads on stdin and writes to stdout.
Native Julia port of chopper $UPSTREAM_VERSION (https://github.com/wdecoster/chopper).

Usage: julia -t <N> --project=<Chopper.jl> bin/chopper.jl [OPTIONS]

Filtering Options:
  -q, --quality <MINQUAL>        Sets a minimum Phred average quality score [default: 0]
      --maxqual <MAXQUAL>        Sets a maximum Phred average quality score [default: 1000]
  -l, --minlength <MINLENGTH>    Sets a minimum read length [default: 1]
      --maxlength <MAXLENGTH>    Sets a maximum read length [default: INF]
      --mingc <MINGC>            Filter min GC content (0 to 1)
      --maxgc <MAXGC>            Filter max GC content (0 to 1)
  -c, --contam <CONTAM>          Filter contaminants against a fasta

Trimming Options:
      --trim-approach <TRIM_APPROACH>
          Select the trimming strategy to apply to the reads:
          - fixed-crop:           Remove a fixed number of bases from the start and/or end
                                  of the read. At least one of --headcrop and --tailcrop
                                  must be greater than 0
          - trim-by-quality:      Trim low-quality bases from the ends of the read until
                                  reaching a base with quality >= --cutoff
          - best-read-segment:    Extract the highest-quality read segment based on
                                  --cutoff, trimming low-quality bases from both ends
          - split-by-low-quality: Split reads by low-quality segments and output
                                  high-quality parts on the left and right, provided they
                                  pass the length filter
      --cutoff <CUTOFF>          Minimum quality score (Q-score) for trimming; required by
                                 trim-by-quality, best-read-segment and split-by-low-quality
      --split-window <SPLIT_WINDOW>
                                 Minimum number of consecutive bases below --cutoff
                                 required to split a read (split-by-low-quality) [default: 1]
      --headcrop <HEADCROP>      Trim N bases from the start of each read (fixed-crop) [default: 0]
      --tailcrop <TAILCROP>      Trim N bases from the end of each read (fixed-crop) [default: 0]
      --update-mods              Recompute the base modification tags (MM, ML and MN)
                                 written by `samtools fastq -T MM,ML,MN` so that they match
                                 the trimmed read. Other tags are passed through unchanged.
                                 Reads whose tags are malformed or do not describe their
                                 sequence are reported as an error

Setup Options:
  -t, --threads <THREADS>        Use N parallel threads, at most the number of Julia
                                 threads (julia -t N) [default: all Julia threads]
  -i, --input <INPUT>            Input filename, plain or compressed [default: read from stdin]
  -o, --output <OUTPUT>          Output filename, gzip-compressed if it ends in .gz
                                 [default: write to stdout]
      --inverse                  Output the opposite of the normal results
      --summary                  Print counts of reads per filter to stderr
  -h, --help                     Print help
  -V, --version                  Print version
"""

# name => (short, kind)
const CLI_OPTIONS = Dict(
    "quality" => ('q', :float), "maxqual" => (nothing, :float),
    "minlength" => ('l', :uint), "maxlength" => (nothing, :uint_inf),
    "mingc" => (nothing, :fraction), "maxgc" => (nothing, :fraction),
    "contam" => ('c', :string), "trim-approach" => (nothing, :trim),
    "cutoff" => (nothing, :u8), "split-window" => (nothing, :uint),
    "headcrop" => (nothing, :uint), "tailcrop" => (nothing, :uint),
    "update-mods" => (nothing, :flag), "threads" => ('t', :uint),
    "input" => ('i', :string), "output" => ('o', :string),
    "inverse" => (nothing, :flag), "summary" => (nothing, :flag),
    "help" => ('h', :flag), "version" => ('V', :flag),
)
const CLI_SHORT = Dict(s => name for (name, (s, _)) in CLI_OPTIONS if s !== nothing)
const TRIM_APPROACHES = ("fixed-crop", "trim-by-quality", "best-read-segment", "split-by-low-quality")

struct CLIError <: Exception
    msg::String
end

function parse_value(name::String, kind::Symbol, v::String)
    bad() = throw(CLIError("invalid value '$v' for '--$name'"))
    if kind === :float
        x = tryparse(Float64, v)
        (x === nothing || any(isspace, v)) && bad()
        return x
    elseif kind === :fraction
        x = tryparse(Float64, v)
        (x === nothing || any(isspace, v)) && bad()
        0 <= x <= 1 || throw(CLIError("invalid value '$v' for '--$name': " *
                                      "Value out of range. Expected a value from 0 to 1."))
        return x
    elseif kind === :uint || kind === :u8
        x = parse_unsigned(v, kind === :u8 ? 255 : typemax(Int))
        x === nothing && bad()
        return x
    elseif kind === :uint_inf
        uppercase(v) == "INF" && return typemax(Int)
        x = parse_unsigned(v)
        x === nothing && bad()
        return x
    elseif kind === :trim
        v in TRIM_APPROACHES || throw(CLIError("invalid value '$v' for '--$name': possible values: " *
                                               join(TRIM_APPROACHES, ", ")))
        return v
    end
    return v
end

# Parses `args` into a Dict of option name => value.
function parse_cli(args::AbstractVector{<:AbstractString})
    opts = Dict{String,Any}()
    i = 1
    set!(name, value) = (haskey(opts, name) &&
        throw(CLIError("the argument '--$name' cannot be used multiple times")); opts[name] = value)
    while i <= length(args)
        a = String(args[i])
        if startswith(a, "--") && length(a) > 2
            name, value = occursin('=', a) ? split(a[3:end], '='; limit = 2) : (a[3:end], nothing)
            name = String(name)
            haskey(CLI_OPTIONS, name) || throw(CLIError("unexpected argument '$a' found"))
        elseif startswith(a, "-") && length(a) >= 2 && a != "--"
            c = a[2]
            haskey(CLI_SHORT, c) || throw(CLIError("unexpected argument '$a' found"))
            name = CLI_SHORT[c]
            rest = a[nextind(a, 2):end]
            startswith(rest, '=') && (rest = rest[2:end])
            value = isempty(rest) ? nothing : rest
        else
            throw(CLIError("unexpected argument '$a' found"))
        end
        kind = CLI_OPTIONS[name][2]
        if kind === :flag
            value === nothing || throw(CLIError("unexpected value '$value' for '--$name' found; no more were expected"))
            set!(name, true)
        else
            if value === nothing
                i += 1
                i <= length(args) || throw(CLIError("a value is required for '--$name' but none was supplied"))
                value = args[i]
            end
            set!(name, parse_value(name, kind, String(value)))
        end
        i += 1
    end
    return opts
end

# The trimming strategy, with chopper's warnings for options that the chosen
# strategy ignores and errors for missing ones.
function cli_trim_strategy(opts::Dict{String,Any}, err::IO)
    approach = get(opts, "trim-approach", nothing)
    cutoff_used = approach in ("trim-by-quality", "best-read-segment", "split-by-low-quality")
    crop_used = approach == "fixed-crop"
    window_used = approach == "split-by-low-quality"
    warn(flag) = println(err, approach === nothing ?
        "Warning: $flag is set but is ignored because no --trim-approach was selected." :
        "Warning: $flag is set but is ignored by the selected trimming approach.")
    haskey(opts, "cutoff") && !cutoff_used && warn("--cutoff")
    get(opts, "headcrop", 0) != 0 && !crop_used && warn("--headcrop")
    get(opts, "tailcrop", 0) != 0 && !crop_used && warn("--tailcrop")
    get(opts, "split-window", 1) != 1 && !window_used && warn("--split-window")
    approach === nothing && return NoTrim()
    if approach == "fixed-crop"
        head, tail = get(opts, "headcrop", 0), get(opts, "tailcrop", 0)
        head == 0 && tail == 0 && throw(CLIError("When using the 'fixed-crop' trimming approach, " *
            "at least one of --headcrop or --tailcrop must be greater than 0."))
        return FixedCrop(head, tail)
    end
    haskey(opts, "cutoff") || throw(CLIError("When using the '$approach' trimming approach, " *
                                             "the --cutoff parameter must be set."))
    cutoff = opts["cutoff"]
    approach == "trim-by-quality" && return TrimByQuality(cutoff)
    approach == "best-read-segment" && return BestReadSegment(cutoff)
    return SplitByLowQuality(cutoff; window = get(opts, "split-window", 1))
end

is_broken_pipe(e) = e isa Base.IOError && e.code == Base.UV_EPIPE

"""
    Chopper.main(args = ARGS; input = stdin, output = stdout, err = stderr) -> Int

Command-line entry point with the options of chopper (`--help` lists them);
returns the exit status. `bin/chopper.jl` calls it.
"""
function main(args::AbstractVector{<:AbstractString} = ARGS; input::IO = stdin, output::IO = stdout,
              err::IO = stderr)
    local opts, trim
    try
        opts = parse_cli(args)
    catch e
        e isa CLIError || rethrow()
        println(err, "Error: ", e.msg, "\n\nFor more information, try '--help'.")
        return 2                     # usage error, as clap
    end
    if get(opts, "help", false)
        print(output, HELP)
        return 0
    elseif get(opts, "version", false)
        println(output, "Chopper.jl ", pkgversion(@__MODULE__), " (chopper ", UPSTREAM_VERSION, ")")
        return 0
    end
    try
        trim = cli_trim_strategy(opts, err)
    catch e
        e isa CLIError || rethrow()
        println(err, "Error: ", e.msg)
        return 1
    end
    threads = get(opts, "threads", Threads.nthreads())
    threads == 0 && (threads = Threads.nthreads())
    if threads > Threads.nthreads()
        println(err, "Warning: $threads threads requested but Julia was started with ",
                Threads.nthreads(), "; start it with `julia -t $threads` to use them.")
    end
    source = get(opts, "input", nothing)
    if source === nothing && input isa Base.TTY
        println(err, "Error: stdin not detected")
        return 1
    end
    if source !== nothing && !isfile(source)
        println(err, "Error: Input file $source is invalid")
        return 1
    end
    contam = get(opts, "contam", nothing)
    if contam !== nothing && !isfile(contam)
        println(err, "Error: Fasta file for filtering contaminants is invalid: $contam")
        return 1
    end
    update_mods = get(opts, "update-mods", false)
    local stats
    try
        filter = ReadFilter(; min_quality = get(opts, "quality", 0.0), max_quality = get(opts, "maxqual", 1000.0),
                            min_length = get(opts, "minlength", 1), max_length = get(opts, "maxlength", typemax(Int)),
                            min_gc = get(opts, "mingc", nothing), max_gc = get(opts, "maxgc", nothing),
                            contaminants = contam)
        stats = chopper(something(source, input), get(opts, "output", output); filter, trim,
                        inverse = get(opts, "inverse", false), update_mods, threads = min(threads, Threads.nthreads()))
    catch e
        if e isa ModTagError
            println(err, "Error: cannot update base modification tags of read ", e.read, ": ", e.msg)
            println(err, haskey(opts, "output") ? "Note: no output file was written." : "Note: the output is incomplete.",
                    " Rerun without --update-mods to write these tags through unchanged.")
            return 1
        elseif is_broken_pipe(e)
            # A closed downstream pipe (`... | head`) is a normal way to stop.
            return 0
        end
        print(err, "Error: ")
        showerror(err, e)
        println(err)
        return 1
    end
    println(err, "Kept ", stats.segments, " reads out of ", stats.reads, " reads")
    get(opts, "summary", false) && (show(err, MIME("text/plain"), stats); println(err))
    if update_mods && stats.reads_with_mods == 0
        println(err, "Warning: --update-mods is set but no read carried an MM tag. ",
                "Did you keep the tags, e.g. `samtools fastq -T MM,ML,MN`?")
    end
    return 0
end
