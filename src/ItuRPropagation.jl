module ItuRPropagation

using Artifacts

export attenuations
export downlinkparameters
export uplinkparameters
export linkparameters

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
export ItuRP372

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
include("iturP372.jl")

include("iturP618.jl")

#endregion include

"""
    attenuations(latlon::LatLon, f::Real, p::Real, θ::Real, D::Real; η::Real=60.0, hs::Union{Real,Missing}=missing, polarization::IturEnum=EnumCircularPolarization)

Computes scintillation attenuation based on Section 2.4.1 of ITU-R P618-13.
    
# Arguments
- `latlon::LatLon`: latitude and longitude (degrees)
- `f::Real`: frequency (GHz)
- `p::Real`: exceedance probability (%)
- `θ::Real`: elevation angle (degrees)
- `D::Real`: antenna diameter (m)
- `η::Real=60`: antenna efficiency (% from 0 to 100, typically 60)
- `hs::Union{Real,Missing}=missing`: altitude of ground station (km)
- `polarization::IturEnum=EnumCircularPolarization`: polarization (EnumHorizontalPolarization, EnumVerticalPolarization, or EnumCircularPolarization)

# Return
- `(Ac=Ac, Ag=Ag, Ar=Ar, As=As, At=At)`: cloud, gas, rain, scintillation, total attenuations (dB)
"""
function attenuations(
    latlon::LatLon,
    f::Real,
    el::Real,
    p::Real,
    ;
    D,
    efficiency = 60, η=efficiency,
    alt = nothing,
    polarization::IturEnum=EnumCircularPolarization
)
    Ac = ItuRP840.cloudattenuation(latlon, f, el, max(5, p))
    Ag = ItuRP676.gaseousattenuation(latlon, f, el, max(5, p); alt)
    Ar = ItuRP618.rainattenuation(latlon, f, el, p, polarization)
    As = ItuRP618.scintillationattenuation(latlon, f, el, p; D, η)

    At = Ag + sqrt((Ar + Ac)^2 + As * As)
    return (Ac=Ac, Ag=Ag, Ar=Ar, As=As, At=At)
end


end # module ItuRPropagation
