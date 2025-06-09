module ItuRP838

#=
Recommendation ITU-R P.838-3 recommends the procedure for obtaining the 
 specfic attenuation (gamma sub R in dB/km) from the rain rate R (mm/h).
=#

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion, _todeg, SUPPRESS_WARNINGS, tilt_from_polarization, IturEnum, _toghz, EnumCircularPolarization, _validel
using Artifacts

const version = ItuRVersion("ITU-R", "P.838", 8, "03/2005")

#region coefficients

# Table 1
const kₕaⱼ = [-5.33980, -0.35351, -0.23789, -0.94158]
const kₕbⱼ = [-0.10008, 1.2697, 0.86036, 0.64552]
const kₕcⱼ = [1.13098, 0.454, 0.15354, 0.16817]
const kₕm = -0.18961
const kₕc = 0.71147

# Table 2
const kᵥaⱼ = [-3.80595, -3.44965, -0.39902, 0.50167]
const kᵥbⱼ = [0.56934, -0.22911, 0.73042, 1.07319]
const kᵥcⱼ = [0.81061, 0.51059, 0.11899, 0.27195]
const kᵥm = -0.16398
const kᵥc = 0.63297

# Table 3
const αₕaⱼ = [-0.14318, 0.29591, 0.32177, -5.37610, 16.1721]
const αₕbⱼ = [1.82442, 0.77564, 0.63773, -0.96230, -3.29980]
const αₕcⱼ = [-0.55187, 0.19822, 0.13164, 1.47828, 3.43990]
const αₕm = 0.67849
const αₕc = -1.95537

# Table 4
const αᵥaⱼ = [-0.07771, 0.56727, -0.20238, -48.2991, 48.5833]
const αᵥbⱼ = [2.3384, 0.95545, 1.1452, 0.791669, 0.791459]
const αᵥcⱼ = [-0.76284, 0.54039, 0.26809, 0.116226, 0.116479]
const αᵥm = -0.053739
const αᵥc = 0.83433

#endregion coefficients

#region internal functions

"""
    _kₕkᵥαₕαᵥ(f)

Computes rain specific attenuation coefficients based on Section 1.

# Arguments
- `f::Float64`: frequency (GHz)

# Return
- `(kₕ::Real, kᵥ::Real, αₕ::Real, αᵥ::Real)`: rain specific attenuation coefficients
"""
function _kₕkᵥαₕαᵥ(f)
    f = _toghz(f)
    logf = log10(f)
    compute_sum(av,bv,cv) = sum(zip(av, bv, cv)) do (a,b,c)
        a * exp(-((logf - b) / c)^2)
    end

    # coefficient k_h based on equation 2
    kₕ = 10^(kₕm * logf + kₕc + compute_sum(kₕaⱼ, kₕbⱼ, kₕcⱼ))

    # coefficient k_v based on equation 2
    kᵥ = 10^(kᵥm * logf + kᵥc + compute_sum(kᵥaⱼ, kᵥbⱼ, kᵥcⱼ))

    # coefficient α_h based on equation 3
    αₕ = αₕm * logf + αₕc + compute_sum(αₕaⱼ, αₕbⱼ, αₕcⱼ)

    # coefficient α_v based on equation 3
    αᵥ = αᵥm * logf + αᵥc + compute_sum(αᵥaⱼ, αᵥbⱼ, αᵥcⱼ)
    return (;kₕ, kᵥ, αₕ, αᵥ)
end

#endregion internal functions

"""
    rainspecificattenuation(f, el; R, polarization::IturEnum, polarization_angle = nothing)

Computes rain specific attenuation for horizontal polarization based on equation 1 of Section 1.

# Arguments
- `f::Float64`: frequency (GHz)
- `el:Float64`: path elevation angle (degrees)

# Keyword Arguments
- `R`: rain rate (mm/hr)
- `polarization::IturEnum`: polarization (EnumHorizontalPolarization, EnumVerticalPolarization, or EnumCircularPolarization). Defaults to circular polarization.
- `polarization_angle`: Polarization tilt angle (degrees) relative to horizontal polarization (τ = 45° for circular polarization)
  - Note: If this is provided, the `polarization` argument is ignored

# Return
- specific attenuation at given rain rate (dB)
"""
function rainspecificattenuation(f, el; R, polarization::IturEnum = EnumCircularPolarization, polarization_angle = nothing)
    f = _toghz(f)
    el = _todeg(el)
    (; γᵣ) = _rainspecificattenuation(f, el; R, polarization, polarization_angle)
    return γᵣ
end

function _rainspecificattenuation(f, el; R, polarization::IturEnum = EnumCircularPolarization, polarization_angle = nothing)
    # Input processing
    polarization_angle = @something(polarization_angle, tilt_from_polarization(polarization)) |> _todeg
    τ = polarization_angle

    cosθ² = cos(el |> deg2rad)^2
    cos2τ = cos(2τ |> deg2rad)
    kₕ, kᵥ, αₕ, αᵥ = _kₕkᵥαₕαᵥ(f)
    k = (kₕ + kᵥ + (kₕ - kᵥ) * cosθ² * cos2τ) / 2
    α = (kₕ * αₕ + kᵥ * αᵥ + (kₕ * αₕ - kᵥ * αᵥ) * cosθ² * cos2τ) / (2k)
    γᵣ = k * R^α
    return (; γᵣ, α, k, kₕ, kᵥ, αₕ, αᵥ)
end

end # module ItuRP838
