module ItuRP618

#=
This Recommendation predicts the various propagation parameters needed in planning Earth-space
systems operating in either the Earth-to-space or space-to-Earth direction.
=#

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion, _tolatlon, _tokm, _todeg, _toghz, SUPPRESS_WARNINGS, IturEnum, tilt_from_polarization
using Artifacts

const version = ItuRVersion("ITU-R", "P.618", 14, "(08/2023)")

const Rₑ = 8500 # effective radius of the Earth (km)

"""
    scintillationattenuation(latlon::LatLon, f::Real, p::Real, θ::Real, D::Real, η::Real=60.0, hL::Real=1000.0)

Computes scintillation attenuation based on Section 2.4.1.
    
# Arguments
- `latlon::LatLon`: latitude and longitude (degrees)
- `f::Real`: frequency (GHz)
- `θ::Real`: elevation angle (degrees)
- `p::Real`: exceedance probability (%)

# Keyword Arguments
- `D::Real=1.0`: antenna diameter (m)
- `η::Real=60`: antenna efficiency (% from 0 to 100, typically 60). Can also be provided with the name `efficiency`.

# Return
- `Ascintillation::Real`: scintillation attenuation (dB)
"""
function scintillationattenuation(
    latlon,
    f,
    el,
    p;
    D=1.0,
    efficiency = 60, η=efficiency,
    warn=!SUPPRESS_WARNINGS[],
)
    (; As) = _scintillationattenuation(latlon, f, el, p; D, η, warn)
    return As
end

function _scintillationattenuation(latlon, f, el, p; D, η, Nwet = nothing, warn=!SUPPRESS_WARNINGS[])
    # We don't process latlon as that is only used in another function that does that
    
    # Inputs checks
    (4 ≤ f ≤ 55) || !warn || @noinline(@warn("ItuR618.scintillationattenuation only supports frequencies between 4 and 55 GHz.\nThe given frequency $f GHz is outside this range. The result may be inaccurate."))
    (0.01 ≤ p ≤ 50) || !warn || @noinline(@warn("ItuR618.scintillationattenuation only supports exceedance probabilities between 0.01% and 50%.\nThe given exceedance probability $p% is outside this range. The result may be inaccurate."))
    # list of parameters in Section 2.4.1
    (5 ≤ el ≤ 90) || !warn || @noinline(@warn("ItuR618.scintillationattenuation only supports elevation angles between 5 and 90 degrees.\nThe given elevation angle $el degrees is outside this range. The result may be inaccurate."))
    
    # step 2
    Nwet = @something(Nwet, ItuRP453.wettermsurfacerefractivityannual_50(latlon))

    # step 3
    σref = 3.6e-3 + 1.0e-4 * Nwet     # equation 40

    # step 4
    sinθ = sin(el |> deg2rad)
    hL = 1000 # Height of turbulent layer (m)
    L = 2 * hL / (sqrt(sinθ^2 + 2.35e-4) + sinθ)     # equation 41

    # step 5
    DeffSquared = (η / 100) * D^2     # based on equation 42

    # step 6
    x = 1.22 * DeffSquared * (f / L)     # equation 43a
    if x >= 7.0
        # from last paragraph of step 6
        # If the argument of the square root is negative (i.e. when x >= 7.0),
        # the predicted scintillation fade depth for any time percentage is zero
        # and the following steps are not required
        # That is: we return 0.0 for the predicted scintillation fade depth
        return 0.0
    end
    # for x < 7.0
    g = sqrt(3.86 * (x^2 + 1)^(11 / 12) * sin(11 / 6 * atan(1, x)) - 7.08 * x^(5 / 6))     # equation 43

    # step 7
    σ = σref * f^(7 / 12) * (g / sinθ^(1.2))     # equation 44

    # step 8
    logp = log10(p)
    aₚ = -0.061 * logp^3 + 0.072 * logp^2 - 1.71 * logp + 3     # equation 45

    # step 9
    As = aₚ * σ     # equation 46
    return (; As, Nwet, σref, L, x, g, σ, aₚ)
end


"""
    rainattenuation(latlon, f, el, p; polarization::IturEnum=EnumCircularPolarization, kwargs...)

Computes rain attenuation based on Section 2.2.1.1.
    
# Arguments
- `latlon`: object representing latitude and longitude, must be convertible to `ItuRPropagation.LatLon`
- `f::Real`: frequency (GHz),
- `el::Real`: elevation angle (degrees)
- `p::Real`: exceedance probability (%)

# Keyword Arguments
- `h_r` or `hᵣ`: Rain height [Km] to be used for the computation. Defaults to `ItuRP1511.topographicheight(lat, lon)`
- `polarization::IturEnum=EnumCircularPolarization`: polarization (EnumHorizontalPolarization, EnumVerticalPolarization, or EnumCircularPolarization) 
  - Note that this last argument is overridden by the keyword argument `polarization_angle` if provided
- `polarization_angle`: Tilt angle [degrees] of the electric field polarization w.r.t. horizontal polarization. Defaults to 45, corresponding to a circularly polarized field.

# Return
- `Aₚ::Float64`: rain attenuation (dB)
"""
function rainattenuation(
    latlon,
    f,
    el,
    p;
    polarization::IturEnum=EnumCircularPolarization,
    polarization_angle = tilt_from_polarization(polarization),
    alt = nothing,
    warn=!SUPPRESS_WARNINGS[],
    hr = nothing, hᵣ = hr,
    R001 = nothing,
)
    # Inputs Processing
    el = _todeg(el)
    f = _toghz(f)

    # Inputs checks
    # first paragraph of 2.2.1.1
    (1 ≤ f ≤ 55) || !warn || @noinline(@warn("ItuR618.rainattenuation only supports frequencies between 1 and 55 GHz.\nThe given frequency $f GHz is outside this range. The result may be inaccurate."))
    
    # per step 10
    (0.001 ≤ p ≤ 5) || !warn || @noinline(@warn("ItuR618.rainattenuation only supports exceedance probabilities between 0.001% and 5%.\nThe given exceedance probability $p% is outside this range."))
    (; Aₚ) = _rainattenuation(latlon, f, el, p; polarization_angle, alt, hᵣ, R001)
    return Aₚ
end

function _rainattenuation(latlon, f, el, p; polarization_angle = nothing, alt = nothing, hᵣ = nothing, R001 = nothing)
    # Inputs Processing
    hᵣ = @something(hᵣ, ItuRP839.rainheightannual(latlon))
    alt = @something(alt, ItuRP1511.topographicheight(latlon))
    polarization_angle = @something(polarization_angle, 45)
    R001 = @something(R001, ItuRP837.rainfallrate001(latlon))

    # We define the early exit output
    early_out = (; Aₚ = 0.0, Lₛ = NaN, Lg = NaN, Lᵣ = NaN, v001 = NaN, A001 = NaN, β = NaN, hᵣ, Lₑ = NaN, γᵣ = NaN, r001 = NaN, R001)
    
    # Early exit
    R001 ≈ 0 && return early_out # This is at the end of Step 4 in 618-14
    altdiff = hᵣ - alt
    altdiff > 0 || return early_out # This is at the end  of Step 2 in 618-14
    
    # Step 2
    sinθ, cosθ = sincos(el |> deg2rad)
    Lₛ = el >= 5 ? altdiff / sinθ : 2 * altdiff / (sqrt(sinθ^2 + (2 * altdiff) / Rₑ) + sinθ)

    # Step 3
    Lg = Lₛ * cosθ

    # step 5
    γᵣ = ItuRP838.rainspecificattenuation(f, el; R=R001, polarization_angle)

    # step 6 - horizontal reduction factor
    r001 = 1 / (1 + 0.78 * sqrt((Lg * γᵣ) / f) - 0.38 * (1 - exp(-2 * Lg)))

    # step 7
    abslat = abs(latlon.lat)
    ζ = rad2deg(atan((altdiff) / (Lg * r001)))
    Lᵣ = ζ > el ? (Lg * r001) / cosθ : altdiff / sinθ
    χ = abslat < 36 ? 36 - abslat : 0.0

    v001 = 1 / (1 + sqrt(sinθ) * (31 * (1 - exp(-el / (1 + χ))) * (sqrt(Lᵣ * γᵣ) / (f * f)) - 0.45))

    # step 8
    Lₑ = Lᵣ * v001

    # step 9
    A001 = γᵣ * Lₑ

    # step 10
    β = if p >= 1.0 || abslat >= 36
        0.0
    elseif p < 1.0 && abslat < 36 && el >= 25
        -0.005 * (abslat - 36)
    else
        -0.005 * (abslat - 36) + 1.8 - 4.25 * sinθ
    end

    Aₚ = A001 * (p / 0.01)^(-(0.655 + 0.033 * log(p) - 0.045 * log(A001) - β * (1 - p) * sinθ))

    out = (; Aₚ, Lₛ, Lg, Lᵣ, v001, A001, β, hᵣ, Lₑ, γᵣ, r001, R001)

    return out
end


"""
    raindiversitygain(f::Real, θ::Real, d::Real, A::Real, Ψ::Real)

Computes rain diversity gain based on Section 2.2.4.2.
    
# Arguments
- `f::Real`: frequency (GHz)
- `θ::Real`: path elevation angle (degrees)
- `d::Real`: separation between the two sites (km)
- `A::Real`: path rain attenuation for a single site (dB)
- `Ψ::Real`: angle made by azimuth of propagation path with respect to the baseline between sites,
             chosen such that Ψ <= 90 (degrees)

             # Return
- G::Real`: net diversity gain (dB)
"""
function raindiversitygain(
   f::Real,
   θ::Real,
   d::Real,
   A::Real,
   Ψ::Real
)
    # first paragraph of 2.2.4.2
    (d < 0 || d > 20) && @warn("ItuR618.raindiversitygain only supports site separations between 1 and 20 km.\nThe given site separation $d km is outside this range.")
    
    # step 1
    a = 0.78*A - 1.94*(1-exp(-0.11*A))
    b = 0.59*(1-exp(-0.1*A))
    Gd = a*(1 - exp(-b*d))     # equation 35

    # step 2
    Gf = exp(-0.025*f)     # equation 36

    # step 3
    Gθ = 1 + 0.006 * θ     # equation 37

    # step 4
    GΨ = 1 + 0.002 * Ψ     # equation 38

    # step 5
    G = Gd * Gf * Gθ * GΨ     # equation 39
    return G
end

"""
    crosspolarizationdiscrimination(f::Real, p::Real, θ::Real, Aᵣ::Real, polarization::IturEnum=EnumCircularPolarization)

Computes cross-polarization discrimination based on Section 4.1.
    
# Arguments
- `f::Real`: frequency (GHz)
- `p::Real`: exceedance probability (%)
- `θ::Real`: path elevation angle (degrees)
- `Aᵣ::Real`: rain attenuation (dB)
- `polarization::IturEnum=EnumCircularPolarization`: polarization (EnumHorizontalPolarization, EnumVerticalPolarization, or EnumCircularPolarization)

# Return
- `XPD::Real`: cross polarization discrimination from rain attenuation statistics (dB)
"""
function crosspolarizationdiscrimination(
    f::Real,
    p::Real,
    θ::Real,
    Aᵣ::Real,
    polarization::IturEnum=EnumCircularPolarization
)
    # first paragraph of 4.1
    (f < 4 || f > 55) && @warn("ItuR618.crosspolarizationdiscrimination only supports frequencies between 4 and 55 GHz.\nThe given frequency $f GHz is outside this range.")
    (θ < 0 || θ > 60) && @warn("ItuR618.crosspolarizationdiscrimination only supports elevation angles between 0 and 60 degrees.\nThe given elevation angle $θ degrees is outside this range.")

    if f > 55
        return 100  # large discrimination
    end

    shouldscale = false
    if 4 <= f < 6
        forig = f
        f = 6.0
        shouldscale = true
    end

    # step 1, equation 65
    logf = log10(f)
    if 6 <= f < 9
        Cf = 60 * logf - 28.3
    elseif 9 <= f < 36
        Cf = 26*logf + 4.1
    elseif 36 <= f <= 55
        Cf = 35.9*logf - 11.3
    end

    # step 2
    if 6 <= f < 9
        Vf = 30.8 * f^(-0.21)        
    elseif 9 <= f < 20
        Vf = 12.8 * f^0.19
    elseif 20 <= f < 40
        Vf = 22.6
    elseif 40 <= f <= 55
        Vf = 13.0 
    end

    Cₐ = Vf * log10(Aᵣ)     # equation 66

    # step 3, equation 67
    if polarization == EnumCircularPolarization
        Cτ = 0
    else
        Cτ = 14.948500216800937  # -10 * log10(1 - 0.484 * (1 + cos(4 * τ))), τ = 0 or 90
    end

    # step 4
    Cθ = -40 * log10(cos(deg2rad(θ)))     # equation 68

    # step 5
    if p <= 0.001
        σ = 15
    elseif p <= 0.01
        σ = 10
    elseif p <= 0.1
        σ = 5
    else
        σ = 0
    end

    Cσ = 0.0053 * σ * σ     # equation 69

    # step 6
    XPDrain = Cf - Cₐ + Cτ + Cθ + Cσ     # equation 70

    # step 7
    Cice = XPDrain * (0.3 + 0.1*log10(p))/2     # equation 71

    # step 8
    XPDp = XPDrain - Cice     # equation 72

    if shouldscale
        XPDp = XPDp - 20*log10(forig/f)
    end
    
    return XPDp
end

end # module ItuRP618
