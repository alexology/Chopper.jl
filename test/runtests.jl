using Chopper
using Test
using CodecZlib, CodecBzip2, CodecXz, CodecZstd

const DATA = joinpath(@__DIR__, "data")
const REF = joinpath(@__DIR__, "reference")

# 0-based half-open (start, end) pairs of chopper's tests -> 1-based ranges.
ranges(pairs) = UnitRange{Int}[(s + 1):e for (s, e) in pairs]

# Runs the command line in-process; returns (exit status, stdout, stderr).
function run_cli(args; input = IOBuffer())
    out, err = IOBuffer(), IOBuffer()
    status = Chopper.main(args; input, output = out, err)
    return status, take!(out), String(take!(err))
end

include("test_units.jl")
include("test_reference.jl")
include("test_io.jl")
include("test_contaminants.jl")
