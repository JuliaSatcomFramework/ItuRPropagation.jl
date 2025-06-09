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

The current implementation only covers points a) and c) in the above list
=#

using ..ItuRPropagation: _todeg, ItuRP835, ItuRP453, ItuRVersion, ItuRP1511, LatLon, ItuRP2145, _toghz
using Artifacts
const version = ItuRVersion("ITU-R", "P.676", 13, "(08/2022)")

# Exports and constructor with separate latitude and longitude arguments
for name in (:gaseousattenuation,)
    @eval $name(lat::Number, lon::Number, args...; kwargs...) = $name(LatLon(lat, lon), args...; kwargs...)
    @eval export $name
end

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
    _Sₒ(a₁, a₂, θ, Pd)

Oxygen line strength - equation 3 of ITU-R P.676-13.

# Arguments
- `a₁`: spectroscopic data for oxygen attenuation from Table 1
- `a₂`: spectroscopic data for oxygen attenuation from Table 1
- `θ`: 300/T (temperature in K)
- `Pd`: dry atmospheric pressure (hPa)

# Return
- `Sᵢ`: line strength of oxygen line
"""
function _Sₒ(a₁, a₂, θ, Pd)
    Sᵢ = a₁ * (1.0e-7 * Pd * θ^3) * exp(a₂ * (1.0 - θ))
    return Sᵢ
end
_Sₒ(r::Table1Row, θ, Pd) = _Sₒ(r.a₁, r.a₂, θ, Pd)


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
- `Pd`: dry atmospheric pressure (hPa)
- `e`: water vapour partial pressure (K * g/m^3)

# Return
- `Fᵢ`: oxygen line shape factor
"""
function _Fₒ(f₀, a₃, a₄, a₅, a₆, f, θ, Pd, e)
    df = a₃ * (1.0e-4) * (Pd * θ^(0.8 - a₄) + (1.1 * e * θ))    # equation 6a - width of the line
    Δf = sqrt(df * df + (2.25e-6))     # equation 6b - accounts for Zeeman splitting of oxygen lines
    δ = (a₅ + (θ * a₆)) * (1.0e-4 * (Pd + e) * θ^0.8)     # equation 7 - correction factor due to interference
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
_Fₒ(r::Table1Row, f, θ, Pd, e) = _Fₒ(r.f₀, r.a₃, r.a₄, r.a₅, r.a₆, f, θ, Pd, e)


"""
    _Sᵥ(b₁, b₂, θ, e)

Water vapour line strength - equation 3 of ITU-R P.676-13.

# Arguments
- `b₁`: spectroscopic data for water vapour attenuation from Table 2
- `b₂`: spectroscopic data for water vapour attenuation from Table 2
- `θ`: 300/T (temperature in K)
- `e`: water vapour partial pressure (K * g/m^3)

# Return
- `Sᵢ`: line strength of water vapour line
"""
function _Sᵥ(b₁, b₂, θ, e)
    Sᵢ = b₁ * 0.1 * e * θ^3.5 * exp(b₂ * (1.0 - θ))
    return Sᵢ
end
_Sᵥ(r::Table2Row, θ, e) = _Sᵥ(r.b₁, r.b₂, θ, e)


"""
    _Fᵥ(f₀, b₃, b₄, b₅, b₆, f, θ, Pd, e)

Water vapour line shape factor - equation 5 of ITU-R P.676-13.

# Arguments
- `f₀`: spectroscopic data for water vapour attenuation from Table 2, water vapour line frequency
- `b₃`: spectroscopic data for water vapour attenuation from Table 2
- `b₄`: spectroscopic data for water vapour attenuation from Table 2
- `b₅`: spectroscopic data for water vapour attenuation from Table 2
- `b₆`: spectroscopic data for water vapour attenuation from Table 2
- `f`: frequency (GHz)
- `θ`: 300/T (temperature in K)
- `Pd`: dry atmospheric pressure (hPa)
- `e`: water vapour partial pressure (K * g/m^3)

# Return
- `Fᵢ`: water vapour line shape factor
"""
function _Fᵥ(f₀, b₃, b₄, b₅, b₆, f, θ, Pd, e)
    Δf = b₃ * 1.0e-4 * (Pd * θ^b₄ + b₅ * e * θ^b₆)     # equation 6a - width of the line
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
_Fᵥ(r::Table2Row, f, θ, Pd, e) = _Fᵥ(r.f₀, r.b₃, r.b₄, r.b₅, r.b₆, f, θ, Pd, e)


"""
    _gammaoxygen(f::Real, T::Real, Pd::Real, ρ::Real)

Oxygen helper for gamma. Part of equation 1 of ITU-R P.676-13.

# Arguments
- `f::Real`: frequency (GHz)
- `T::Real`: absolute temperature (K)
- `Pd::Real`: dry atmospheric pressure (hPa)
- `ρ::Real`: water vapour density (g/m^3)

# Return
- `γₒ`: specific attenuation due to dry air (oxygen, pressure-induced nitrogen, and non-resonant Debye attenuation)
"""
function _gammaoxygen(
    f::Real,
    T::Real,
    Pd::Real,
    ρ::Real
)
    theta = 300 / T     # where portion of equation 3
    e = ρ * T / 216.7     # equation 4

    d = 5.6e-4 * (Pd + e) * theta^0.8      # equation 9

    # equation 8
    NppD = f * Pd * theta^2 * (
               (6.14e-5) / (d * (1 + (f / d)^2))
               +
               (1.4e-12 * Pd * theta^1.5) / (1 + 1.9e-5 * f^1.5)
           )

    Nppₒ = NppD
    for row in TABLE1_DATA
        S = _Sₒ(row, theta, Pd) # Equation 3
        F = _Fₒ(row, f, theta, Pd, e) # Equation 5
        Nppₒ += S * F # Equation 2a
    end

    γₒ = 0.1820 * f * Nppₒ
    return (; γₒ, Nppₒ, d)
end


"""
    _gammawater(f::Real, Pd::Real, T::Real, ρ::Real)

Water vapour helper for gamma. Part of equation 1 of ITU-R P.676-13.

# Arguments
- `f::Real`: frequency (GHz)
- `T::Real`: absolute temperature (K)
- `Pd::Real`: dry atmospheric pressure (hPa)
- `ρ::Real`: water vapour density (g/m^3)

# Return
- `γᵥ`: specific attenuation due to water vapour
"""
function _gammawater(
    f::Real,
    T::Real,
    Pd::Real,
    ρ::Real
)
    theta = 300 / T     # where portion of equation 3    
    e = ρ * T / 216.7     # equation 4

    Nppᵥ = 0.0
    for row in TABLE2_DATA
        S = _Sᵥ(row, theta, e) # Equation 3
        F = _Fᵥ(row, f, theta, Pd, e) # Equation 5
        Nppᵥ += S * F # Equation 2b
    end
    γᵥ = 0.1820 * f * Nppᵥ
    return (; γᵥ, Nppᵥ)
end


"""
    _gamma(f::Real, Pd::Real, T::Real, ρ::Real)

Primarily to check against Itu-Rpy values. Equation 1 of ITU-R P.676-13.

# Arguments
- `f::Real`: frequency (GHz)
- `T::Real`: absolute temperature (K)
- `Pd::Real`: dry atmospheric pressure (hPa)
- `ρ::Real`: water vapour density (g/m^3)

# Return
- `gamma`: specific gaseous attenuation
"""
function _gamma(
    f::Real,
    T::Real,
    Pd::Real,
    ρ::Real
)
    outₒ = _gammaoxygen(f, T, Pd, ρ)
    outᵥ = _gammawater(f, T, Pd, ρ)
    γ = outₒ.γₒ + outᵥ.γᵥ
    return (; γ, outₒ..., outᵥ...)
end

# This is just a function to extract the idx of the Part1 row whose frequency is just below the given `f`

@inline function _part1_idx(f)
    # Doing searchsorted on a steprange of 0.5 is much slower than doing it on an integer steprange. So we just index on 1:700 with 2 * f. We also have to add 1 to the result if f ≥ 118.75 as the oxygen data in Part 1 has an additional line at 118.75 GHz.
    return searchsortedlast(1:700, 2f - 1) + (f ≥ 118.75 ? 1 : 0)
end

@inline function _part2_idx(f)
    # Doing searchsorted on a steprange of 0.5 is much slower than doing it on an integer steprange. So we just index on 1:700 with 2 * f.
    return searchsortedlast(1:700, 2f - 1)
end

#= 
This is the function hₒ which can be used for equation 31 and 34 of Annex 2 of ITU-R P.676-13.
The inputs are:
# Arguments
- `f`: frequency (GHz)

# Keyword arguments
- `P`: total pressure (hPa)
- `T`: temperature (K)
- `ρ`: water vapour density (g/m^3)
=#
function _hₒ(f; P, T, ρ)
    idx = _part1_idx(f)
    below = PART1_DATA[idx]
    above = PART1_DATA[min(idx + 1, end)]
    # Now we do linear interpolation
    δf = (f - below.f) / (above.f - below.f)
    aₒ = δf * (above.aₒ - below.aₒ) + below.aₒ
    bₒ = δf * (above.bₒ - below.bₒ) + below.bₒ
    cₒ = δf * (above.cₒ - below.cₒ) + below.cₒ
    dₒ = δf * (above.dₒ - below.dₒ) + below.dₒ
    hₒ = aₒ + bₒ * T + cₒ * P + dₒ * ρ
    return (; hₒ, aₒ, bₒ, cₒ, dₒ, idx, below, above)
end

#= 
This is the function hₒ which can be used for equation 31 and 34 of Annex 2 of ITU-R P.676-13.
The inputs are:
# Arguments
- `f`: frequency (GHz)

# Keyword arguments
- `P`: total pressure (hPa)
- `T`: temperature (K)
- `ρ`: water vapour density (g/m^3)
=#
function _Kᵥ(f; P, T, ρ)
    idx = _part2_idx(f)
    below = PART2_DATA[idx]
    above = PART2_DATA[min(idx + 1, end)]
    # Now we do linear interpolation
    δf = (f - below.f) / (above.f - below.f)
    aᵥ = δf * (above.aᵥ - below.aᵥ) + below.aᵥ
    bᵥ = δf * (above.bᵥ - below.bᵥ) + below.bᵥ
    cᵥ = δf * (above.cᵥ - below.cᵥ) + below.cᵥ
    dᵥ = δf * (above.dᵥ - below.dᵥ) + below.dᵥ
    Kᵥ = aᵥ + bᵥ * ρ + cᵥ * T + dᵥ * P
    return (; Kᵥ, aᵥ, bᵥ, cᵥ, dᵥ, idx, below, above)
end

# #endregion internal use functions

# #region initialization

const R_E = 6371.0     # radius of the Earth (km) for use in P676 calculation

# Struct used to store the data used to compute the gas attenuation for a slant path as defined in Section 2.2 of ITU-R P.676-13.
struct SlantPathLayer
    "Layer thickness (km)"
    δ::Float64
    "Height of the bottom part of layer (km), i.e. the altitude at which the layer starts"
    h::Float64
    "Distance from the center of the earth at layer bottom (km)"
    r::Float64
    "Temperature at layer midpoint (K)"
    T::Float64
    "Total pressure at layer midpoint (hPa)"
    P::Float64
    "Water vapour density at layer midpoint (g/m^3)"
    ρ::Float64
    "Water vapour partial pressure at layer midpoint (hPa)"
    e::Float64
    "Dry atmospheric pressure at layer midpoint (hPa)"
    Pd::Float64
    "Refractive index at layer midpoint"
    n::Float64
end
function SlantPathLayer(; δ, h, r=nothing, T=nothing, P=nothing, ρ=nothing, n=nothing)
    r = @something(r, h + R_E)
    h′ = h + δ / 2
    T = @something(T, ItuRP835.standardtemperature(h′))
    P = @something(P, ItuRP835.standardpressure(h′))
    ρ = @something(ρ, ItuRP835.standardwatervapourdensity(h′; T, P))
    e = ρ * T / 216.7     # equation 4
    Pd = P - e     # must be changed to dry air pressure
    n = @something(n, ItuRP453.radiorefractiveindex(T, Pd, e))
    return SlantPathLayer(δ, h, r, T, P, ρ, e, Pd, n)
end
function SlantPathLayer(i::Integer)
    coeff = exp((i - 1) / 100)
    δ = 0.0001 * coeff # Equation 14
    h = if i === 1
        0.0
    else
        0.0001 * (coeff - 1) / (exp(0.01) - 1) # Equation 15
    end
    return SlantPathLayer(; δ, h)
end

# These are the standard layers to be used for slant path gaseous attenuation prediction as per Section 2.2 of ITU-R P.676-13.
const STANDARD_LAYERS = map(SlantPathLayer, 1:922)

# This computes a single term in the sum of equation 13 of ITU-R P.676-13.
function layerattenuation(layer::SlantPathLayer, f, el; r₁=first(STANDARD_LAYERS).r, n₁=first(STANDARD_LAYERS).n, sinβ₁ = nothing) 
    # β = 90 - el so sin(β) = cos(el)
    sinβ₁ = @something(sinβ₁, cos(_todeg(el) |> deg2rad))
    (; δ, r, n, T, Pd, ρ) = layer
    sinβ = sinβ₁ * n₁ * r₁ / (n * r)
    cos²β = 1 - sinβ^2
    cosβ = sqrt(cos²β)
    β = asin(sinβ)
    a = -r * cosβ + sqrt(r^2 * cos²β + 2r * δ + δ^2)
    (; γ, γₒ, γᵥ) = _gamma(f, T, Pd, ρ) # Gamma expects the try air pressure, not the total so we have to pass Pd
    att = a * γ
    return (; att, a, γ, γₒ, γᵥ, β)
end
# #endregion initialization

# Computes the gaseous attenuation for a slant bath as per Section 2.2.1 of ITU-R P.676-13.
function _gasattenuation_layers(layers::Vector{SlantPathLayer}, f, el)
    r₁ = first(layers).r
    n₁ = first(layers).n
    sinβ₁ = cos(el |> deg2rad)
    att = 0.0
    for layer in layers
        outs = layerattenuation(layer, f, el; r₁, n₁, sinβ₁)
        att += outs.att
    end
    return att
end

#=
    _slantoxygenattenuation(latlon, f, el, p; γₒ = nothing, alt = nothing)
    _slantoxygenattenuation(latlon, f, el; γₒ = nothing, alt = nothing, P, T, ρ)

Computes the slant path gaseous attenuation for oxygen as per equation 32 in Section 1.2, Annex 2 of ITU-R P.676-13.

# Arguments
- `latlon::LatLon`: latitude and longitude (degrees)
- `f`: frequency (GHz)
- `el`: elevation angle (degrees)
- `p`: exceedance probability 1-100 (%). If provided (first method), the pressure, temperature and gas vapour density are computed at the target probability based on the provided location using the ITU-R P.2145 algorithms.

# Keyword arguments
- `γₒ`: specific attenuation due to oxygen (dB/km), computed based on mean atmospheric conditions at the provided location if not provided
- `alt`: Altitude of the provided location above sea level (km). If not provided, it is computed from the provided location using the ITU-R P.1511 model
- `P`: Surface total pressure at the desired exceedance probability (hPa), at the desired location.
- `T`: Surface temperature at the desired exceedance probability (K), at the desired location.
- `ρ`: Surface water vapour density at the desired exceedance probability (g/m^3), at the desired location.
=#
function _slantoxygenattenuation(latlon::LatLon, f, el; γₒ, P, T, ρ)
    (; hₒ, aₒ, bₒ, cₒ, dₒ) = _hₒ(f; P, T, ρ)
    sinθ = sin(el |> deg2rad)
    Aₒ_zenith = γₒ * hₒ
    Aₒ = Aₒ_zenith / sinθ
    return (; Aₒ, Aₒ_zenith, γₒ, hₒ, aₒ, bₒ, cₒ, dₒ, P, T, ρ)
end
function _slantoxygenattenuation(latlon, f, el, p; γₒ=nothing, alt=nothing)
    (; P, T, ρ, alt) = ItuRP2145.annual_surface_values(latlon, p; alt)
    γₒ = @something γₒ let
        params = ItuRP2145.annual_surface_values(latlon; alt)
        P̄ = params.P
        T̄ = params.T
        ρ̄ = params.ρ
        ē = ρ̄ * T̄ / 216.7
        P̄d = P̄ - ē
        _gammaoxygen(f, T̄, P̄d, ρ̄).γₒ
    end
    _slantoxygenattenuation(latlon, f, el; γₒ, P, T, ρ)
end

# This function compute the slant water vapour attenuation for a given frequency, elevation angle, and atmospheric conditions as per equation 40 of Section 2.3 in Annex 2 of ITU-R P.676-13. The first method expects all relevant surface statistical parameters to be provided as keyword arguments. The second method instead computes them based on the given exceedance probability `p` and location/altitude.
function _slantwaterattenuation(latlon, f, el; P, T, ρ, V)
    (; Kᵥ, aᵥ, bᵥ, cᵥ, dᵥ) = _Kᵥ(f; P, T, ρ)
    sinθ = sin(el |> deg2rad)
    Aᵥ_zenith = Kᵥ * V
    Aᵥ = Aᵥ_zenith / sinθ
    return (; Aᵥ, Aᵥ_zenith, Kᵥ, aᵥ, bᵥ, cᵥ, dᵥ, P, T, ρ, V)
end
function _slantwaterattenuation(latlon, f, el, p; alt=nothing)
    (; P, T, ρ, alt) = ItuRP2145.annual_surface_values(latlon; alt)
    V = ItuRP2145.surfacewatervapourcontentannual(latlon, p; alt)
    _slantwaterattenuation(latlon, f, el; P, T, ρ, V)
end


"""
    Ag = gaseousattenuation(latlon, f, el, p; alt = nothing, gamma_oxygen = nothing, γₒ = nothing)

Computes the statistical gaseous attenuation for a slant path following the approximate computation specified in Annex 2 of ITU-R P.676-13.

More specifically this computes ``Ag = Ao + Aw`` implementing:
- Equation 32 in Section 1.2 for Oxygen attenuation `Ao`
- Equation 40 in Section 2.3 for Water vapour attenuation `Aw`

# Arguments
- `latlon`: Object specifying the latitude and longitude of the location of interest, must be an object that can be converted to an instance of `ItuRPropagation.LatLon`
- `f`: frequency (GHz)
- `el`: elevation angle (degrees)

# Keyword arguments
- `alt`: Altitude at the provided location, to be used for computing the various intermediate varaibles. If not provided, default to the altitude computed with `ItuRP1511.topographicheight`
- `gamma_oxygen` (or `γₒ`): Specific attenuation due to oxygen (dB/km) computed from the average surface conditions, at the desired location. **This is computed automatically based on other inputs if not explicitly provided, but it is only location dependent and >70% of the time is spent in computing this, so consider precomputing and passing it directly for maximum speed**

See the extended help for the signature of the function with precomputed intermediate variables.

# Extended help

## Alternative method

When the values (both average and ccdf) of Pressure, Temeprature and water vapour density are already available, the computation can be made faster by using this alternative method that takes all the useful parameters as keyword arguments.

This method can also be used to simulated oxygen/gas attenuations based on instantaneous values of the surface parameters as per Sections 1.1 and 2.2 of Annex 2.
In this case, the instantaneous values of P, T and ρ should be explicitly provided for the ccdf values with the same name, but also to the kwargs related to mean values (i.e. P̄, T̄ and ρ̄ respectively).

    gaseousattenuation(latlon, f, el; kwargs...)

An additional method is available, for cases where the statistical values of the surface parameters are already available and do not need to be computed internally.

This method do not accept the outage probability `p` as last argument but expects he following keyword arguments specifying the various surface paramters (both in mean value and in ccdf statistical value):
- `P_mean` (or `P̄`): Average surface total pressure (hPa) at the desired location. In case this method is used for simulating attenuation based on instantaneous values, the instantaneous total pressure `P` should be provided also as value to this kwarg.
- `T_mean` (or `T̄`): Average surface temperature (K) at the desired location. In case this method is used for simulating attenuation based on instantaneous values, the instantaneous temperature `T` should be provided also as value to this kwarg.
- `rho_mean` (or `ρ̄`): Average surface water vapour density (g/m^3) at the desired location. In case this method is used for simulating attenuation based on instantaneous values, the instantaneous water vapour density `ρ` should be provided also as value to this kwarg.
- `P`: Surface total pressure (hPa) at the desired exceedance probability, at the desired location.
- `T`: Surface temperature (K) at the desired exceedance probability, at the desired location.
- `rho`: Surface water vapour density (g/m^3) at the desired exceedance probability, at the desired location.
- `V`: Surface water vapour content (kg/m^2) at the desired exceedance probability, at the desired location.
- `gamma_oxygen` (or `γₒ`): Specific attenuation due to oxygen (dB/km) computed from the average surface conditions, at the desired location. **This is computed automatically based on other inputs if not explicitly provided, but it is only location dependent and >70% of the time is spent in computing this, so consider precomputing and passing it directly for maximum speed**
"""
function gaseousattenuation(latlon, f, el; 
    P_mean = nothing, P̄ = nothing, # Average surface total pressure at the desired location
    T_mean = nothing, T̄ = nothing, # Average surface temperature at the desired location
    rho_mean = nothing, ρ̄ = nothing, # Average surface water vapour density at the desired location
    P, # Surface total pressure (hPa) at the desired exceedance probability, at the desired location.
    T, # Surface temperature (K) at the desired exceedance probability, at the desired location.
    rho = nothing, ρ = nothing, # Surface water vapour density (g/m^3) at the desired exceedance probability, at the desired location.
    V, # Surface water vapour content (kg/m^2) at the desired exceedance probability, at the desired location.
    gamma_oxygen = nothing, γₒ = nothing # Specific attenuation due to oxygen (dB/km) computed from the average surface conditions, at the desired location.
    )
    f = _toghz(f)
    1 ≤ f ≤ 350 || throw(ArgumentError("Frequency must be between 1 and 350 GHz for the computation of gaseous attenuation"))
    P̄ = @something P̄ P_mean throw(ArgumentError("The average total surafce pressure has to be provided using either the `P_mean` or `P̄` keyword argument"))
    T̄ = @something T̄ T_mean throw(ArgumentError("The average surface temperature has to be provided using either the `T_mean` or `T̄` keyword argument"))
    ρ̄ = @something ρ̄ rho_mean throw(ArgumentError("The average surface water vapour density has to be provided using either the `rho_mean` or `ρ̄` keyword argument"))
    ρ = @something ρ rho throw(ArgumentError("The ccdf value of the surface water vapour density has to be provided using either the `rho` or `ρ` keyword argument"))
    γₒ = @something γₒ gamma_oxygen let
        ē = ρ̄ * T̄ / 216.7
        P̄d = P̄ - ē
        _gammaoxygen(f, T̄, P̄d, ρ̄).γₒ
    end
    (; Aₒ) = _slantoxygenattenuation(latlon, f, el; γₒ, P, T, ρ)
    (; Aᵥ) = _slantwaterattenuation(latlon, f, el; P = P̄, T = T̄, ρ = ρ̄, V)
    Agas = Aₒ + Aᵥ
    return Agas
end

function gaseousattenuation(latlon, f, el, p; alt=nothing, gamma_oxygen=nothing, γₒ=nothing)
    mean_vals = ItuRP2145.annual_surface_values(latlon; alt)
    P̄ = mean_vals.P
    T̄ = mean_vals.T
    ρ̄ = mean_vals.ρ
    alt = mean_vals.alt
    (; P, T, ρ) = ItuRP2145.annual_surface_values(latlon, p; alt)
    V = ItuRP2145.surfacewatervapourcontentannual(latlon, p; alt)
    gaseousattenuation(latlon, f, el; P̄, T̄, ρ̄, V, P, T, ρ, gamma_oxygen, γₒ)
end

end # module ItuRP676