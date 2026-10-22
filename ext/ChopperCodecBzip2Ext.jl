module ChopperCodecBzip2Ext

using Chopper: Chopper
using CodecBzip2: Bzip2DecompressorStream

function __init__()
    Chopper.DECOMPRESSORS[:bzip2] = Bzip2DecompressorStream
    return
end

end # module
