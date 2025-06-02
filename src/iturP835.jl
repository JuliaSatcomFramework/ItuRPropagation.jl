module ItuRP835

#=
Recommendation ITU-R P.835 provides expressions and data for reference standard atmospheres required for
the calculation of gaseous attenuation on Earth-space paths.
=#

using ..ItuRPropagation
using Artifacts

const version = ItuRVersion("ITU-R", "P.835", 7, "(08/2024)")

# Geometric height to geopotential height
geopotentialheight(z::Real) = z * 6356.766 / (z + 6356.766)

"""
    standardtemperature(Z::Float64)

Standard temperature for geometric heights <= 100 km, based on equations 2a-2g and 4a-4b of Section 1.1

# Arguments
- `Z::Float64`: geometric height (km)

# Return
- temperature (°K)
"""
function standardtemperature(Z::Real)
    0 ≤ Z ≤ 100 || throw(ArgumentError("Z=$Z\n (geometric height) in IturP835.standardtemperature must be between 0 and 100"))

    hp = geopotentialheight(Z)
    if 0 <= hp <= 11
        return (288.15 - 6.5 * hp)
    elseif 11 < hp <= 20
        return 216.65
    elseif 20 < hp <= 32
        return (216.65 + (hp - 20))
    elseif 32 < hp <= 47
        return (228.65 + 2.8(hp - 32))
    elseif 47 < hp <= 51
        return 270.65
    elseif 51 < hp <= 71
        return (270.65 - 2.8 * (hp - 51))
    elseif 71 < hp <= 84.852
        return (214.65 - 2(hp - 71))
    elseif 86 <= Z <= 91
        return 186.8673
    elseif 91 < Z <= 100
        return (263.1905 - 76.3232 * sqrt(1 - ((Z - 91) / 19.9429)^2))
    end
end


"""
    standardpressure(Z::Float64)

Standard pressure for geometric heights <= 100 km, based on equations 3a-3g and 5 of Section 1.1.

# Arguments
- `Z::Float64`: geometric height (km)

# Return
- dry pressure (hPa)
"""
function standardpressure(Z::Real)
    0 ≤ Z ≤ 100 || throw(ArgumentError("Z=$Z\n (geometric height) in IturP835.standardpressure must be between 0 and 100"))
    hp = geopotentialheight(Z)
    if 0 <= hp <= 11
        return (1013.25 * (288.15 / (288.15 - 6.5 * hp))^(-34.1632 / 6.5))
    elseif 11 < hp <= 20
        return (226.3226 * exp(-34.1632 * (hp - 11) / 216.65))
    elseif 20 < hp <= 32
        return (54.74980 * (216.65 / (216.65 + (hp - 20)))^34.1632)
    elseif 32 < hp <= 47
        return (8.680422 * (228.65 / (228.65 + 2.8 * (hp - 32)))^(34.1632 / 2.8))
    elseif 47 < hp <= 51
        return (1.109106 * exp(-34.1632 * (hp - 47) / 270.65))
    elseif 51 < hp <= 71
        return (0.6694167 * (270.65 / (270.65 - 2.8 * (hp - 51)))^(-34.1632 / 2.8))
    elseif 71 < hp <= 84.852
        return (0.03956649 * (214.65 / (214.65 - 2.0 * (hp - 71)))^(-34.1632 / 2.0))
    elseif 86 <= Z <= 100
        a₀ = 95.571899
        a₁ = -4.011801
        a₂ = 6.424731e-2
        a₃ = -4.789660e-4
        a₄ = 1.340543e-6
        return exp(a₀ + a₁ * Z + a₂ * Z^2 + a₃ * Z^3 + a₄ * Z^4)
    end
end


"""
    standardwatervapourdensity(Z; kwargs...)

Standard water vapour density for geometric heights <= 100 km, based on equations 6-8 of Section 1.2.

# Arguments
- `Z::Real`: geometric height (km)

# Keyword arguments
- `rho_0::Real=7.5`: standard ground level water vapour density (g/m^3). Can also be provided as `ρ₀`.
- `T::Real`: Temperature (°K) at given height. Computed with `standardtemperature(Z)` if not provided.
- `P::Real`: Total barometric pressure (hPa) at given height. Computed with `standardpressure(Z)` if not provided.
"""
function standardwatervapourdensity(
        Z::Real; 
        rho_0 = nothing,
        ρ₀ = nothing,
        T = standardtemperature(Z),
        P = standardpressure(Z),
)
    ρ₀ = @something ρ₀ rho_0 7.5
    # equation 6 where h₀ = 2; equation 7 where ρ₀=7.5
    ρ = ρ₀ * exp(-Z / 2)

    # see paragraph below equation 8 regarding mixing ratio
    e = ρ * T / 216.7
    mixingratio = e / P
    if mixingratio < 2e-6
        # see paragraph below equation 8 regarding mixing ratio
        # and recalculate ρ
        enew = P * 2e-6
        ρ = enew * 216.7 / T
    end
    return ρ
end

end # module ItuRP835
