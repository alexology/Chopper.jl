# Command-line chopper:
#
#     julia -t 4 --project=path/to/Chopper.jl path/to/Chopper.jl/bin/chopper.jl -q 10 -l 500 -i reads.fastq.gz > filtered.fastq
#
# `--help` lists the options (those of chopper 0.14.1, plus -o/--output and
# --summary).

using Chopper

exit(Chopper.main(ARGS))
