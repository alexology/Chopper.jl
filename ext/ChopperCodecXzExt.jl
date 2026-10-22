module ChopperCodecXzExt

using Chopper: Chopper
using CodecXz: XzDecompressorStream

function __init__()
    Chopper.DECOMPRESSORS[:xz] = XzDecompressorStream
    return
end

end # module
