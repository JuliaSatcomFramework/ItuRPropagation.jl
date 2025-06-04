module ItuRP676

#=
This Recommendation provides:
a) methods in Annex 1 to calculate the slant path gaseous attenuation, phase nonlinearity, atmospheric
    bending, excess atmospheric path length and downwelling and upwelling noise temperatures due to oxygen
    and water vapour for the frequency range from 1 to 1 000 GHz for arbitrary known pressure, temperature and
    water vapour height profiles;
b) an approximate method in Annex 2 to estimate the instantaneous slant path gaseous attenuation due
    to oxygen and water vapour for the frequency range from 1 to 350 GHz when the instantaneous surface
    pressure, surface temperature and surface water vapour density or integrated water vapour content1 are known
    from local data, a reference profile, or referenced digital maps;
c) an approximate method in Annex 2 to estimate the statistics of slant path gaseous attenuation due to
    oxygen and water vapour for the frequency range from 1 to 350 GHz when the surface pressure, surface
    temperature and surface water vapour density or integrated water vapour content statistics are known either
    from local data, a reference profile, or referenced digital maps;
d) a Weibull approximation to the slant path water vapour attenuation for use in
    Recommendation ITU-R P.1853.
=#

using ..ItuRPropagation
using Artifacts
const version = ItuRVersion("ITU-R", "P.676", 13, "(08/2022)")

#region coefficients

# from Table 1 of Annex 1, Section 1 of ITU-R P.676-13
# spectroscopic data for oxygen
struct Table1Row
    f₀::Float64
    a₁::Float64
    a₂::Float64
    a₃::Float64
    a₄::Float64
    a₅::Float64
    a₆::Float64
end
const numberoxygencoefficients = 44
const TABLE1_DATA = [
    Table1Row(50.474214, 0.975, 9.651, 6.69, 0, 2.566, 6.850),
    Table1Row(50.987745, 2.529, 8.653, 7.17, 0, 2.246, 6.800),
    Table1Row(51.50336, 6.193, 7.709, 7.64, 0, 1.947, 6.729),
    Table1Row(52.021429, 14.32, 6.819, 8.11, 0, 1.667, 6.640),
    Table1Row(52.542418, 31.24, 5.983, 8.58, 0, 1.388, 6.526),
    Table1Row(53.066934, 64.29, 5.201, 9.06, 0, 1.349, 6.206),
    Table1Row(53.595775, 124.6, 4.474, 9.55, 0, 2.227, 5.085),
    Table1Row(54.130025, 227.3, 3.8, 9.96, 0, 3.170, 3.750),
    Table1Row(54.67118, 389.7, 3.182, 10.37, 0, 3.558, 2.654),
    Table1Row(55.221384, 627.1, 2.618, 10.89, 0, 2.560, 2.952),
    Table1Row(55.783815, 945.3, 2.109, 11.34, 0, -1.172, 6.135),
    Table1Row(56.264774, 543.4, 0.014, 17.03, 0, 3.525, -0.978),
    Table1Row(56.363399, 1331.8, 1.654, 11.89, 0, -2.378, 6.547),
    Table1Row(56.968211, 1746.6, 1.255, 12.23, 0, -3.545, 6.451),
    Table1Row(57.612486, 2120.1, 0.91, 12.62, 0, -5.416, 6.056),
    Table1Row(58.323877, 2363.7, 0.621, 12.95, 0, -1.932, 0.436),
    Table1Row(58.446588, 1442.1, 0.083, 14.91, 0, 6.768, -1.273),
    Table1Row(59.164204, 2379.9, 0.387, 13.53, 0, -6.561, 2.309),
    Table1Row(59.590983, 2090.7, 0.207, 14.08, 0, 6.957, -0.776),
    Table1Row(60.306056, 2103.4, 0.207, 14.15, 0, -6.395, 0.699),
    Table1Row(60.434778, 2438, 0.386, 13.39, 0, 6.342, -2.825),
    Table1Row(61.150562, 2479.5, 0.621, 12.92, 0, 1.014, -0.584),
    Table1Row(61.800158, 2275.9, 0.91, 12.63, 0, 5.014, -6.619),
    Table1Row(62.41122, 1915.4, 1.255, 12.17, 0, 3.029, -6.759),
    Table1Row(62.486253, 1503, 0.083, 15.13, 0, -4.499, 0.844),
    Table1Row(62.997984, 1490.2, 1.654, 11.74, 0, 1.856, -6.675),
    Table1Row(63.568526, 1078, 2.108, 11.34, 0, 0.658, -6.139),
    Table1Row(64.127775, 728.7, 2.617, 10.88, 0, -3.036, -2.895),
    Table1Row(64.67891, 461.3, 3.181, 10.38, 0, -3.968, -2.590),
    Table1Row(65.224078, 274, 3.8, 9.96, 0, -3.528, -3.680),
    Table1Row(65.764779, 153, 4.473, 9.55, 0, -2.548, -5.002),
    Table1Row(66.302096, 80.4, 5.2, 9.06, 0, -1.660, -6.091),
    Table1Row(66.836834, 39.8, 5.982, 8.58, 0, -1.680, -6.393),
    Table1Row(67.369601, 18.56, 6.818, 8.11, 0, -1.956, -6.475),
    Table1Row(67.900868, 8.172, 7.708, 7.64, 0, -2.216, -6.545),
    Table1Row(68.431006, 3.397, 8.652, 7.17, 0, -2.492, -6.600),
    Table1Row(68.960312, 1.334, 9.65, 6.69, 0, -2.773, -6.650),
    Table1Row(118.750334, 940.3, 0.01, 16.64, 0, -0.439, 0.079),
    Table1Row(368.498246, 67.4, 0.048, 16.4, 0, 0.000, 0.000),
    Table1Row(424.76302, 637.7, 0.044, 16.4, 0, 0.000, 0.000),
    Table1Row(487.249273, 237.4, 0.049, 16, 0, 0.000, 0.000),
    Table1Row(715.392902, 98.1, 0.145, 16, 0, 0.000, 0.000),
    Table1Row(773.83949, 572.3, 0.141, 16.2, 0, 0.000, 0.000),
    Table1Row(834.145546, 183.1, 0.145, 14.7, 0, 0.000, 0.000)
]

# from Table 2 of Annex 1, Section 1 of ITU-R P.676-13
# spectroscopic data for water vapour
struct Table2Row
    f₀::Float64 # water vapour line frequency [GHz]
    b₁::Float64
    b₂::Float64
    b₃::Float64
    b₄::Float64
    b₅::Float64
    b₆::Float64
end
const TABLE2_DATA = [
    Table2Row(22.23508, 0.1079, 2.144, 26.38, 0.76, 5.087, 1),
    Table2Row(67.80396, 0.0011, 8.732, 28.58, 0.69, 4.93, 0.82),
    Table2Row(119.99594, 0.0007, 8.353, 29.48, 0.7, 4.78, 0.79),
    Table2Row(183.310087, 2.273, 0.668, 29.06, 0.77, 5.022, 0.85),
    Table2Row(321.22563, 0.047, 6.179, 24.04, 0.67, 4.398, 0.54),
    Table2Row(325.152888, 1.514, 1.541, 28.23, 0.64, 4.893, 0.74),
    Table2Row(336.227764, 0.001, 9.825, 26.93, 0.69, 4.74, 0.61),
    Table2Row(380.197353, 11.67, 1.048, 28.11, 0.54, 5.063, 0.89),
    Table2Row(390.134508, 0.0045, 7.347, 21.52, 0.63, 4.81, 0.55),
    Table2Row(437.346667, 0.0632, 5.048, 18.45, 0.6, 4.23, 0.48),
    Table2Row(439.150807, 0.9098, 3.595, 20.07, 0.63, 4.483, 0.52),
    Table2Row(443.018343, 0.192, 5.048, 15.55, 0.6, 5.083, 0.5),
    Table2Row(448.001085, 10.41, 1.405, 25.64, 0.66, 5.028, 0.67),
    Table2Row(470.888999, 0.3254, 3.597, 21.34, 0.66, 4.506, 0.65),
    Table2Row(474.689092, 1.26, 2.379, 23.2, 0.65, 4.804, 0.64),
    Table2Row(488.490108, 0.2529, 2.852, 25.86, 0.69, 5.201, 0.72),
    Table2Row(503.568532, 0.0372, 6.731, 16.12, 0.61, 3.98, 0.43),
    Table2Row(504.482692, 0.0124, 6.731, 16.12, 0.61, 4.01, 0.45),
    Table2Row(547.67644, 0.9785, 0.158, 26, 0.7, 4.5, 1),
    Table2Row(552.02096, 0.184, 0.158, 26, 0.7, 4.5, 1),
    Table2Row(556.935985, 497, 0.159, 30.86, 0.69, 4.552, 1),
    Table2Row(620.700807, 5.015, 2.391, 24.38, 0.71, 4.856, 0.68),
    Table2Row(645.766085, 0.0067, 8.633, 18, 0.6, 4, 0.5),
    Table2Row(658.00528, 0.2732, 7.816, 32.1, 0.69, 4.14, 1),
    Table2Row(752.033113, 243.4, 0.396, 30.86, 0.68, 4.352, 0.84),
    Table2Row(841.051732, 0.0134, 8.177, 15.9, 0.33, 5.76, 0.45),
    Table2Row(859.965698, 0.1325, 8.055, 30.6, 0.68, 4.09, 0.84),
    Table2Row(899.303175, 0.0547, 7.914, 29.85, 0.68, 4.53, 0.9),
    Table2Row(902.611085, 0.0386, 8.429, 28.65, 0.7, 5.1, 0.95),
    Table2Row(906.205957, 0.1836, 5.11, 24.08, 0.7, 4.7, 0.53),
    Table2Row(916.171582, 8.4, 1.441, 26.73, 0.7, 5.15, 0.78),
    Table2Row(923.112692, 0.0079, 10.293, 29, 0.7, 5, 0.8),
    Table2Row(970.315022, 9.009, 1.919, 25.5, 0.64, 4.94, 0.67),
    Table2Row(987.926764, 134.6, 0.257, 29.85, 0.68, 4.55, 0.9),
    Table2Row(1780, 17506, 0.952, 196.3, 2, 24.15, 5)
]

# from Part 1 data file of Annex 2, Section 1.2 of ITU-R P.676-13
# coefficients for oxygen attenuation prediction
struct Part1Row
    f::Float64
    aₒ::Float64
    bₒ::Float64
    cₒ::Float64
    dₒ::Float64
end
const PART1_DATA = map(eachrow(read!(artifact"p676/Part1.bin", zeros(700, 5)))) do row
    Part1Row(row...)
end

# from Part 2 data file of Annex 2, Section 2.2 of ITU-R P.676-13
# coefficients for water vapour attenuation prediction
struct Part2Row
    f::Float64
    aᵥ::Float64
    bᵥ::Float64
    cᵥ::Float64
    dᵥ::Float64
end
const PART2_DATA = map(eachrow(read!(artifact"p676/Part2.bin", zeros(699, 5)))) do row
    Part2Row(row...)
end

#endregion coefficients

#region internal use functions

"""
    _Sₒ(a₁, a₂, θ, P)

Oxygen line strength - equation 3 of ITU-R P.676-13.

# Arguments
- `a₁`: spectroscopic data for oxygen attenuation from Table 1
- `a₂`: spectroscopic data for oxygen attenuation from Table 1
- `θ`: 300/T (temperature in K)
- `P`: dry atmospheric pressure (hPa)

# Return
- `Sᵢ`: line strength of oxygen line
"""
function _Sₒ(a₁, a₂, θ, P)
    Sᵢ = a₁ * (1.0e-7 * P * θ^3) * exp(a₂ * (1.0 - θ))
    return Sᵢ
end
_Sₒ(r::Table1Row, θ, P) = _Sₒ(r.a₁, r.a₂, θ, P)


"""
    _Fₒ(f₀, a₃, a₄, a₅, a₆, f, θ, P, e)

Oxygen line shape factor - equation 5 of ITU-R P.676-13.

# Arguments
- `f₀`: spectroscopic data for oxygen attenuation from Table 1, oxygen line frequency
- `a₃`: spectroscopic data for oxygen attenuation from Table 1
- `a₄`: spectroscopic data for oxygen attenuation from Table 1
- `a₅`: spectroscopic data for oxygen attenuation from Table 1
- `a₆`: spectroscopic data for oxygen attenuation from Table 1
- `f`: frequency (GHz)
- `θ`: 300/T (temperature in K)
- `P`: dry atmospheric pressure (hPa)
- `e`: water vapour partial pressure (K * g/m^3)

# Return
- `Fᵢ`: oxygen line shape factor
"""
function _Fₒ(f₀, a₃, a₄, a₅, a₆, f, θ, P, e)
    df = a₃ * (1.0e-4) * (P * θ^(0.8 - a₄) + (1.1 * e * θ))    # equation 6a - width of the line
    Δf = sqrt(df * df + (2.25e-6))     # equation 6b - accounts for Zeeman splitting of oxygen lines
    δ = (a₅ + (θ * a₆)) * (1.0e-4 * (P + e) * θ^0.8)     # equation 7 - correction factor due to interference
    # effects in oxygen lines
    Fᵢ = (f / f₀) * (
        (
            (Δf - (δ * (f₀ - f))) / ((f₀ - f)^2 + Δf^2)
        )
        +
        (
            (Δf - (δ * (f₀ + f))) / ((f₀ + f)^2 + Δf^2)
        )
    )
    return Fᵢ
end
_Fₒ(r::Table1Row, f, θ, P, e) = _Fₒ(r.f₀, r.a₃, r.a₄, r.a₅, r.a₆, f, θ, P, e)


"""
    _Sᵥ(b₁, b₂, θ, e)

Water vapour line strength - equation 3 of ITU-R P.676-13.

# Arguments
- `b₁`: spectroscopic data for water vapour attenuation from Table 2
- `b₂`: spectroscopic data for water vapour attenuation from Table 2
- `θ`: 300/T (temperature in K)
- `P`: dry atmospheric pressure (hPa)

# Return
- `Sᵢ`: line strength of water vapour line
"""
function _Sᵥ(b₁, b₂, θ, e)
    Sᵢ = b₁ * 0.1 * e * θ^3.5 * exp(b₂ * (1.0 - θ))
    return Sᵢ
end
_Sᵥ(r::Table2Row, θ, e) = _Sᵥ(r.b₁, r.b₂, θ, e)


"""
    _Fᵥ(fₒ, b₃, b₄, b₅, b₆, f, θ, P, e)

Water vapour line shape factor - equation 5 of ITU-R P.676-13.

# Arguments
- `fₒ`: spectroscopic data for water vapour attenuation from Table 2, water vapour line frequency
- `b₃`: spectroscopic data for water vapour attenuation from Table 2
- `b₄`: spectroscopic data for water vapour attenuation from Table 2
- `b₅`: spectroscopic data for water vapour attenuation from Table 2
- `b₆`: spectroscopic data for water vapour attenuation from Table 2
- `f`: frequency (GHz)
- `θ`: 300/T (temperature in K)
- `P`: dry atmospheric pressure (hPa)
- `e`: water vapour partial pressure (K * g/m^3)

# Return
- `Fᵢ`: water vapour line shape factor
"""
function _Fᵥ(f₀, b₃, b₄, b₅, b₆, f, θ, P, e)
    Δf = b₃ * 1.0e-4 * (P * θ^b₄ + b₅ * e * θ^b₆)     # equation 6a - width of the line
    Δf = 0.535 * Δf + sqrt(0.217 * Δf * Δf + (2.1316e-12 * f₀^2 / θ))     # equation 6b - accounts for Doppler
    # broadening of water vapour lines
    # δ = 0     # equation 7
    Fᵢ = f / f₀ * (
        (
            (Δf) / ((f₀ - f)^2 + Δf^2)
        )
        +
        (
            (Δf) / ((f₀ + f)^2 + Δf^2)
        )
    )
    return Fᵢ
end
_Fᵥ(r::Table2Row, f, θ, P, e) = _Fᵥ(r.f₀, r.b₃, r.b₄, r.b₅, r.b₆, f, θ, P, e)


"""
    _gammaoxygen(f::Real, T::Real, P::Real, ρ::Real)

Oxygen helper for gamma. Part of equation 1 of ITU-R P.676-13.

# Arguments
- `f::Real`: frequency (GHz)
- `T::Real`: absolute temperature (K)
- `P::Real`: dry atmospheric pressure (hPa)
- `ρ::Real`: water vapour density (g/m^3)

# Return
- `γₒ`: specific attenuation due to dry air (oxygen, pressure-induced nitrogen, and non-resonant Debye attenuation)
"""
function _gammaoxygen(
    f::Real,
    T::Real,
    P::Real,
    ρ::Real
)
    theta = 300 / T     # where portion of equation 3
    e = ρ * T / 216.7     # equation 4

    d = 5.6e-4 * (P + e) * theta^0.8      # equation 9

    # equation 8
    NppD = f * P * theta^2 * (
               (6.14e-5) / (d * (1 + (f / d)^2))
               +
               (1.4e-12 * P * theta^1.5) / (1 + 1.9e-5 * f^1.5)
           )

    Nppₒ = NppD
    for row in TABLE1_DATA
        S = _Sₒ(row, theta, P) # Equation 3
        F = _Fₒ(row, f, theta, P, e) # Equation 5
        Nppₒ += S * F # Equation 2a
    end

    γₒ =  0.1820 * f * Nppₒ
    return (; γₒ, Nppₒ, d)
end


"""
    _gammawater(f::Real, P::Real, T::Real, ρ::Real)

Water vapour helper for gamma. Part of equation 1 of ITU-R P.676-13.

# Arguments
- `f::Real`: frequency (GHz)
- `T::Real`: absolute temperature (K)
- `P::Real`: dry atmospheric pressure (hPa)
- `ρ::Real`: water vapour density (g/m^3)

# Return
- `γᵥ`: specific attenuation due to water vapour
"""
function _gammawater(
    f::Real,
    T::Real,
    P::Real,
    ρ::Real
)
    theta = 300 / T     # where portion of equation 3    
    e = ρ * T / 216.7     # equation 4

    Nppᵥ = 0.0
    for row in TABLE2_DATA
        S = _Sᵥ(row, theta, e) # Equation 3
        F = _Fᵥ(row, f, theta, P, e) # Equation 5
        Nppᵥ += S * F # Equation 2b
    end
    γᵥ = 0.1820 * f * Nppᵥ
    return (; γᵥ, Nppᵥ)
end


"""
    _gamma(f::Real, P::Real, T::Real, ρ::Real)

Primarily to check against Itu-Rpy values. Equation 1 of ITU-R P.676-13.

# Arguments
- `f::Real`: frequency (GHz)
- `T::Real`: absolute temperature (K)
- `P::Real`: dry atmospheric pressure (hPa)
- `ρ::Real`: water vapour density (g/m^3)

# Return
- `gamma`: specific gaseous attenuation
"""
function _gamma(
    f::Real,
    T::Real,
    P::Real,
    ρ::Real
)
    outₒ = _gammaoxygen(f, T, P, ρ)
    outᵥ = _gammawater(f, T, P, ρ)
    γ = outₒ.γₒ + outᵥ.γᵥ
    return (; γ, outₒ..., outᵥ...)
end

# #endregion internal use functions

# #region initialization

# const δᵢ = 0.0001 .* exp.(([1:922;] .- 1) ./ 100)     # equation 14
# const δᵢ2 = δᵢ .* δᵢ
# const hᵢ = (δᵢ .- 0.0001) ./ (exp(0.01) - 1) .+ (δᵢ ./ 2)     # equation 15
# const rᵢ = hᵢ .+ (6371) .- (δᵢ ./ 2)

# const Tᵢ = ItuRP835.standardtemperature.(hᵢ)
# const Pᵢ = ItuRP835.standardpressure.(hᵢ)

# #endregion initialization


# """
#     function gaseousattenuation(latlon::LatLon, f::Real, p::Real, θ::Real, hs::Real=missing)

# Computes exact slant path gaseous attenuation based on Annex 1. δᵢ and hᵢ are computed
# using Equations 14 and 15.

# Per equation 62 of ITU-R P.618-13, a large part of the cloud attenuation and gaseous attenuation is already 
# included in the rain attenuation prediction for time percentages below 1%.

# # Arguments
# - `latlon::LatLon`: latitude and longitude (degrees)
# - `f::Real`: frequency (GHz)
# - `p::Real`: exceedance probability 1-100 (%)
# - `θ::Real`: elevation angle (degrees)
# - `hs::Real=missing`: altitude of ground station (km)

# # Return
# - `Agas`: gaseous attenuation (dB)
# """
# function gaseousattenuation(
#     latlon::LatLon,
#     f::Real,
#     p::Real,
#     θ::Real,
#     hs::Union{Real,Missing}=missing
# )
#     ρ₀ = ItuRP2145.surfacewatervapourdensityannual(latlon, p, hs)     # equation 62 of ItuRP618
#     ρᵢ = map(zip(hᵢ, Tᵢ, Pᵢ)) do (Z, T, P)
#         ItuRP835.standardwatervapourdensity(Z; T, P, ρ₀)
#     end
#     eᵢ = (Tᵢ .* ρᵢ) ./ (216.7)     # equation 4
#     dryPᵢ = Pᵢ - eᵢ     # must be changed to dry air pressure

#     nᵢ = ItuRP453.radiorefractiveindex.(Tᵢ, dryPᵢ, eᵢ)

#     β₁ = π / 2 - θ * π / 180
#     n₁ = nᵢ[1]
#     r₁ = rᵢ[1]

#     βᵢ = asin.((n₁ * r₁ * sin(β₁)) ./ (nᵢ .* rᵢ))     # equation 19b
#     cosβ = cos.(βᵢ)

#     # equation 17
#     aᵢ = -1 .* rᵢ .* cos.(βᵢ) .+ sqrt.(rᵢ .* rᵢ .* cosβ .* cosβ + (2) .* rᵢ .* δᵢ + δᵢ2)
#     γᵢ = _gammaoxygen.(f, Tᵢ, dryPᵢ, ρᵢ) .+ _gammawater.(f, Tᵢ, dryPᵢ, ρᵢ)

#     Agas = sum(aᵢ .* γᵢ)     # equation 13
#     return Agas
# end

end # module ItuRP676
