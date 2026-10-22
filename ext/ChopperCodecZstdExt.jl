module ChopperCodecZstdExt

using Chopper: Chopper
using CodecZstd: ZstdDecompressorStream

function __init__()
    Chopper.DECOMPRESSORS[:zstd] = ZstdDecompressorStream
    return
end

end # module
