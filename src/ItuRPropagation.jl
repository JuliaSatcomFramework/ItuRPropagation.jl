module ItuRPropagation

using Artifacts

export attenuations

#region package export

export ItuRP840

export ItuRP453
export ItuRP1144
export ItuRP1511
export ItuRP2145
export ItuRP835
export ItuRP838
export ItuRP839
export ItuRP837
export ItuRP676

export ItuRP618

#endregion package export

#region include

include("iturcommon.jl")

include("iturP1144.jl") # Interpolations

include("iturP840.jl")

include("iturP453.jl")
include("iturP1511.jl")
include("iturP2145.jl")
include("iturP835.jl")
include("iturP838.jl")
include("iturP839.jl")
include("iturP837.jl")
include("iturP676.jl")

include("iturP618.jl")

using .ItuRP618: attenuations

#endregion include


end # module ItuRPropagation
