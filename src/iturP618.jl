module ItuRP618

#=
This Recommendation predicts the various propagation parameters needed in planning Earth-space
systems operating in either the Earth-to-space or space-to-Earth direction.
=#

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion, tolatlon, _tokm, _todeg, _toghz, SUPPRESS_WARNINGS, IturEnum, tilt_from_polarization, _validel, ItuRP840, ItuRP676, ItuRP453, ItuRP838, ItuRP837, ItuRP1511, ItuRP839, EnumCircularPolarization
using Artifacts

const version = ItuRVersion("ITU-R", "P.618", 14, "(08/2023)")

# Exports and constructor with separate latitude and longitude arguments
for name in (:scintillationattenuation, :rainattenuation, :attenuations)
    @eval $name(lat::Number, lon::Number, args...; kwargs...) = $name(LatLon(lat, lon), args...; kwargs...)
    @eval export $name
end

const Rₑ = 8500 # effective radius of the Earth (km)

"""
    scintillationattenuation(latlon::LatLon, f::Real, el::Real, p::Real; D::Real=1.0, η::Real=60.0)

Computes scintillation attenuation based on Section 2.4.1.
    
# Arguments
- `latlon`: object representing latitude and longitude, must be convertible to `ItuRPropagation.LatLon`
- `f::Real`: frequency (GHz)
- `el::Real`: elevation angle (degrees)
- `p::Real`: exceedance probability (%)

# Keyword Arguments
- `D::Real=1.0`: antenna diameter (m)
- `η::Real=60`: antenna efficiency (% from 0 to 100, typically 60). Can also be provided with the name `efficiency`.
- `warn`: Whether to warn if the inputs are outside the supported range. Defaults to `!SUPPRESS_WARNINGS[]`

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
    # Input preprocessing
    latlon = tolatlon(latlon)
    f = _toghz(f)
    el = _todeg(el)
    
    (; As) = _scintillationattenuation(latlon, f, el, p; D, η, warn)
    return As
end

function _scintillationattenuation(latlon, f, el, p; D, η, Nwet = nothing, warn=!SUPPRESS_WARNINGS[], ignore_pwarn = false)
    # We don't process latlon as that is only used in another function that does that
    
    # Inputs checks
    (0 ≤ η ≤ 100) || throw(ArgumentError("ItuR618.scintillationattenuation only supports antenna efficiencies between 0 and 100%.\nThe given antenna efficiency $η% is outside this range."))
    (η > 1) || !warn || @noinline(@warn("ItuR618.scintillationattenuation expects antenna efficiencies in percentage (between 0 and 100).\nThe given antenna efficiency $η% is below 1% and seems to indicate that efficiency was not provided as percentage."))
    (4 ≤ f ≤ 55) || !warn || @noinline(@warn("ItuR618.scintillationattenuation only supports frequencies between 4 and 55 GHz.\nThe given frequency $f GHz is outside this range. The result may be inaccurate."))
    (0.01 ≤ p ≤ 50) || ignore_pwarn || !warn || @noinline(@warn("ItuR618.scintillationattenuation only supports exceedance probabilities between 0.01% and 50%.\nThe given exceedance probability $p% is outside this range. The result may be inaccurate."))
    _validel(el) # Check that elevation is between 0 and 90 degrees
    (5 ≤ el) || !warn || @noinline(@warn("ItuR618.scintillationattenuation only supports elevation angles between 5 and 90 degrees.\nThe given elevation angle $el degrees is outside this range. The result may be inaccurate."))
    
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
- `alt`: Altitude [km] of the receiver. Defaults to `ItuRP1511.topographicheight(latlon)`
- `R001`: Annual rain rate [mm/h] exceeded 0.01% of the time. Defaults to `ItuRP837.rainfallrate001(latlon)`
- `warn`: Whether to warn if the inputs are outside the supported range. Defaults to `!SUPPRESS_WARNINGS[]`

# Return
- `Ap::Float64`: rain attenuation (dB)
"""
function rainattenuation(
    latlon,
    f,
    el,
    p;
    polarization::IturEnum=EnumCircularPolarization,
    polarization_angle = tilt_from_polarization(polarization),
    alt = nothing,
    hr = nothing, hᵣ = hr,
    R001 = nothing,
    warn=!SUPPRESS_WARNINGS[],
)
    # Inputs Processing
    latlon = tolatlon(latlon)
    el = _todeg(el)
    f = _toghz(f)

    (; Aₚ) = _rainattenuation(latlon, f, el, p; polarization_angle, alt, hᵣ, R001, warn)
    return Aₚ
end

function _rainattenuation(latlon, f, el, p; polarization_angle = nothing, alt = nothing, hᵣ = nothing, R001 = nothing, warn=!SUPPRESS_WARNINGS[])
    # Inputs Processing
    hᵣ = @something(hᵣ, ItuRP839.rainheightannual(latlon)) |> _tokm
    alt = @something(alt, ItuRP1511.topographicheight(latlon)) |> _tokm
    polarization_angle = @something(polarization_angle, 45) |> _todeg
    R001 = @something(R001, ItuRP837.rainfallrate001(latlon))

    # Inputs checks
    (1 ≤ f ≤ 55) || !warn || @noinline(@warn("ItuR618.rainattenuation only supports frequencies between 1 and 55 GHz.\nThe given frequency $f GHz is outside this range. The result may be inaccurate."))
    (0.001 ≤ p ≤ 5) || !warn || @noinline(@warn("ItuR618.rainattenuation only supports exceedance probabilities between 0.001% and 5%.\nThe given exceedance probability $p% is outside this range. The result may be inaccurate."))
    _validel(el) # Check that elevation is between 0 and 90 degrees

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
    combine_attenuations(; Ac::Real, Ag::Real, Ar::Real, As::Real)

Combines the various attenuations into a total attenuation following the formula of ITU-R P618-14 (Eq 65-66):

```math
At = Ag + \\sqrt{(Ar + Ac)^2 + As^2}
```

This function does not take the outage probability as input, it just assumes that `Ar` will be 0 for `p > 5%`.

# Arguments
- `Ac::Real`: Cloud attenuation (dB). Defaults to 0 if not provided
- `Ag::Real`: Gaseous attenuation (dB). Defaults to 0 if not provided
- `Ar::Real`: Rain attenuation (dB). Defaults to 0 if not provided
- `As::Real`: Scintillation attenuation (dB). Defaults to 0 if not provided

# Return
- `At::Real`: Total attenuation (dB)
"""
function combine_attenuations(; Ac::Real = 0, Ag::Real = 0, Ar::Real = 0, As::Real = 0)
    At = Ag + sqrt((Ar + Ac)^2 + As^2)
    return (; At, Ac, Ag, Ar, As)
end

"""
    attenuations(latlon, f, el, p; D, η, alt, polarization, polarization_angle, warn, extra_kwargs...)

Computes the various attenuations for earth/space links based on Section 2.4.1 of ITU-R P618-13.
    
# Arguments
- `latlon`: object representing latitude and longitude, must be convertible to `ItuRPropagation.LatLon`
- `f`: frequency (GHz)
- `el`: elevation angle (degrees)
- `p`: exceedance probability (%)

# Keyword Arguments
- `D`: antenna diameter (m)
- `η` (also accepted as `efficiency`): antenna efficiency (% from 0 to 100, defaults to 60).
- `alt`: altitude of ground station (km). Defaults to `ItuRP1511.topographicheight(latlon)`
- `polarization`: polarization (EnumHorizontalPolarization, EnumVerticalPolarization, or EnumCircularPolarization)
- `polarization_angle`: tilt angle [degrees] of the electric field polarization w.r.t. horizontal polarization. **Note: This field is computed from the `polarization` argument if not provided, but this function disregards the `polarization` argument if `polarization_angle` is explicitly provided.**
- `warn`: Whether to warn if the inputs are outside the supported range. Defaults to `!SUPPRESS_WARNINGS[]`

See extended help for more details on extra keyword arguments and for maximizing the speed of the computation when some intermediate inputs are known.

# Return
- `(Ac=Ac, Ag=Ag, Ar=Ar, As=As, At=At)`: cloud, gas, rain, scintillation, total attenuations (dB)

# Extended help

This function accepts some additional keyword arguments that are intermediate inputs to the computation, and can be provided to maximize the speed of the computation. These `extra_kwargs` are:
- `Ag_zenith`: Zenith gaseous attenuation (dB/km) at zenith (90° elevation)
- `Ac_zenith`: Zenith cloud attenuation (dB/km) at zenith (90° elevation)
- `Nwet`: Wet term surface refractivity (N-units)
- `R001`: Annual rain rate [mm/h] exceeded 0.01% of the time.
- `hᵣ`: Rain height [Km] to be used for the computation.
- `gamma_oxygen` (also accepted as `γₒ`): Specific attenuation due to oxygen (dB/km) computed from the average surface conditions, at the desired location. **Note: This is used to speed up the computation of `Ag_zenith` but is ignored if the `Ag_zenith` argument is provided**

Additionally, the function can be called with `Val(true)` as last positional argument (after `p`) to also return the various intermediate inputs used in the computation for speeding up subsequent calls where only the elevation angle changes. Here is an example on how to use this when needing to compute attenuations for a given ground location towards multiple satellites at different elevation angles:

## Fast computation example

```julia
els = [10, 20, 30, 40, 50, 60, 70, 80, 90]
ll = LatLon(0,0)
p = .5 # .5% outage
D = 1 # 1m antenna diameter
# We compute the attenuations for the first elevation angle and save also the various intermediate values stored in the `kwargs` field of the output.
nt = ItuRP618.attenuations(ll, 30, first(els), p, Val(true); D)
# We create a vector of outputs for different elevations exploiting the precomputed kwargs
outs = map(eachindex(els)) do i
    i == 1 && return nt.attenuations
    ItuRP618.attenuations(ll, 30, els[i], p; nt.kwargs...)
end
```
"""
function attenuations(
    latlon, f, el, p, ::Val{full_output} = Val(false);
    D,
    efficiency = 60, η=efficiency,
    alt = nothing,
    polarization::IturEnum=EnumCircularPolarization,
    polarization_angle = tilt_from_polarization(polarization),
    warn=!SUPPRESS_WARNINGS[],
    Ag_zenith = nothing,
    Ac_zenith = nothing,
    Nwet = nothing,
    R001 = nothing,
    hᵣ = nothing,
    gamma_oxygen = nothing, γₒ = nothing,
) where {full_output}
    latlon = tolatlon(latlon)
    el = _todeg(el)
    f = _toghz(f)
    _validel(el) # Validate input elevation
    full = _attenuations(latlon, f, el, p; D, η, polarization_angle, alt, warn, Ag_zenith, Ac_zenith, Nwet, R001, hᵣ, gamma_oxygen, γₒ)
    if full_output
        return full
    else
        return full.attenuations
    end
end

function _attenuations(latlon, f, el, p;
    D,
    Ag_zenith = nothing,
    Ac_zenith = nothing,
    η = 60,
    polarization_angle = nothing,
    alt = nothing,
    Nwet = nothing,
    R001 = nothing,
    hᵣ = nothing,
    gamma_oxygen = nothing, γₒ = nothing,
    warn=!SUPPRESS_WARNINGS[],
)
    # Check the positional arguments ranges
    5 ≤ el ≤ 90 || !warn || @noinline(@warn("ItuR840.cloudattenuation only supports elevation angles between 5 and 90 degrees.\nThe given elevation angle $el degrees is outside this range so results may be inaccurate."))

    # Extract altitude
    alt = @something(alt, ItuRP1511.topographicheight(latlon)) |> _tokm
    # Default polarization angle to 45 (circular)
    polarization_angle = @something(polarization_angle, 45) |> _todeg
    # Extract the zenith gaseous attenuation
    Ag_zenith = @something(Ag_zenith, ItuRP676.gaseousattenuation(latlon, f, 90, max(5, p); alt, gamma_oxygen, γₒ))
    # Extract the zenith cloud attenuation
    Ac_zenith = @something(Ac_zenith, ItuRP840.cloudattenuation(latlon, f, 90, max(5, p); warn))

    # Get the scintillation attenuation
    (; As, Nwet) = _scintillationattenuation(latlon, f, el, p; D, η, Nwet, warn, ignore_pwarn = p < 0.01) # We ignore pwarning because we might have values below 0.01%

    (; Aₚ, hᵣ, R001) = if p ≤ 5 
        _rainattenuation(latlon, f, el, p; polarization_angle, alt, hᵣ, R001, warn)
    else
        hᵣ = @something(hᵣ, ItuRP839.rainheightannual(latlon)) |> _tokm
        R001 = @something(R001, ItuRP837.rainfallrate001(latlon))
        (; Aₚ = 0.0, hᵣ, R001)
    end

    sinθ = sin(el |> deg2rad) # For clouds and gas we start from zenith values
    Ac = Ac_zenith / sinθ
    Ag = Ag_zenith / sinθ

    attenuations = combine_attenuations(; Ac, Ag, Ar = Aₚ, As)
    kwargs = (; Ac_zenith, Ag_zenith, polarization_angle, alt, D, η, Nwet, R001, hᵣ, warn)

    return (; attenuations, kwargs, inps = (; latlon, f, el, p))
end

end # module ItuRP618
