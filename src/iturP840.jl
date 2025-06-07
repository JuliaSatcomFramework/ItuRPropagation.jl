module ItuRP840

#=
This Recommendation provides methods to predict the attenuation due to clouds and fog on Earth-space paths.
=#

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion, _tolatlon, _tokm, _torad, _toghz, SUPPRESS_WARNINGS
using ..ItuRP1144: ItuRP1144, AbstractSquareGridITP, SquareGridData, SquareGridStatisticalData
using Artifacts: Artifacts, @artifact_str

const version = ItuRVersion("ITU-R", "P.840", 9, "(08/2023)")

# Exports and constructor with separate latitude and longitude arguments
for name in (:liquidwatercontent, :cloudattenuation)
    @eval $name(lat::Number, lon::Number, args...; kwargs...) = $name(LatLon(lat, lon), args...; kwargs...)
    @eval export $name
end

#region initialization

const δlat = 0.25
const δlon = 0.25
const latrange = range(-90, 90, step=δlat)
const lonrange = range(-180, 180, step=δlon)
const datasize = (length(latrange), length(lonrange))

# exceedance probabilities
const psannual = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 60, 70, 80, 90, 95, 99]

# exceedance probability values for reading files
const filespsannual = ["001", "002", "003", "005", "01", "02", "03", "05", "1", "2", "3", "5", "10", "20", "30", "50", "60", "70", "80", "90", "95", "99"]

const SGD_TYPE = let
    T = Float64
    R = typeof(latrange)
    SquareGridData{T, R, String}
end

# This will hold the mean and annual ccdf data of integrated liquid water content for computation of the cloud attenuation. The underlying data is already a interpolator that will use bilinear interpolation to find the value at the given location
@kwdef mutable struct LiquidContentData
    ccdf::Union{Nothing, SquareGridStatisticalData{SGD_TYPE}} = nothing
    mean::Union{Nothing, SGD_TYPE} = nothing
end

const ANNUAL_DATA = LiquidContentData()

function initialize!()
    isnothing(ANNUAL_DATA.mean) || return nothing
    @info "P840: Loading annual data of Integrated Liquid Water Content"
    Lmean = read!(artifact"p840_annual/L_mean.bin", zeros(datasize))
    ANNUAL_DATA.mean = SquareGridData(latrange, lonrange, Lmean, "Mean Integrated Liquid Water Content")
    items = map(zip(psannual, filespsannual)) do (p, suffix)
        data = read!(joinpath(artifact"p840_annual", "L_$(suffix).bin"), zeros(datasize))
        name = "Integrated Liquid Water Content at $p% exceedance probability"
        SquareGridData(latrange, lonrange, data, name)
    end
    ANNUAL_DATA.ccdf = SquareGridStatisticalData(psannual, items)
    return nothing
end


#endregion initialization

#region internal functions

"""
    _Kₗ(f::Real, T::Real)

Computes cloud specific attenuation coefficient based on Section 2 (Equation 2). 
    
# Arguments
- `f::Real`: frequency (GHz)
- `T::Real`: temperature (°C)

# Return
- `K::Real`: attenuation coefficient ((dB/km)/(g/m^3))
"""
function _Kₗ( 
    f::Real,
    T::Real=273.75,
)
    coeff = 300 / T - 1 # Term repeated in many places
    ϵ₀ = 77.66 + 103.3 * coeff     # equation 6
    ϵ₁ = 0.0671 * ϵ₀     # equation 7
    ϵ₂ = 3.52     # equation 8
    fₚ = 20.2 - 146 * coeff + 316 * coeff^2     # equation 9
    fₛ = 39.8 * fₚ     # equation 11

    rfₚ = f / fₚ
    rfₛ = f / fₛ
    # equation 4
    ϵ′′ = (
        ((f * (ϵ₀ - ϵ₁)) / (fₚ * (1 + rfₚ^2)))
        +
        ((f * (ϵ₁ - ϵ₂)) / (fₛ * (1 + rfₛ^2)))
    )

    # equation 5
    ϵ′ = (
        (ϵ₀ - ϵ₁) / (1 + rfₚ^2)
        +
        (ϵ₁ - ϵ₂) / (1 + rfₛ^2)
        + ϵ₂
    )


    η = (2 + ϵ′) / ϵ′′     # equation 3
    Kₗ = 0.819 * f / (ϵ′′ * (1 + η^2))     # equation 2
    return (; Kₗ, η, ϵ′′, ϵ′)
end

# This is implementing Equation 12/14
function _K_L(f)
    nt = _Kₗ(f, 273.75)
    A₁ = 0.1522
    A₂ = 11.51
    A₃ = -10.4912
    f₁ = -23.9589
    f₂ = 219.2096
    σ₁ = 3.2991e3
    σ₂ = 2.7595e6
    K_L = nt.Kₗ * (
        A₁ * exp(-(f - f₁)^2 / σ₁) +
        A₂ * exp(-(f - f₂)^2 / σ₂) +
        A₃
    )
    return (; K_L, nt...)
end

#endregion internal functions

"""
    liquidwatercontent(latlon, p)

Computes the integrated liquid water content at a given location and exceedance probability based on the digital annual maps in Part 1 of the Recommendation P.840-8.

# Arguments
- `latlon`: object representing latitude and longitude, must be convertible to `ItuRPropagation.LatLon`
- `p`: exceedance probability (%)   
"""
function liquidwatercontent(latlon, p; warn=!SUPPRESS_WARNINGS[])
    itp = @something(ANNUAL_DATA.ccdf,let
        initialize!()
        ANNUAL_DATA.ccdf
    end)::SquareGridStatisticalData{SGD_TYPE}
    return itp(latlon, p; warn, kind = "the integrated liquid water content")
end

"""
    cloudattenuation(latlon, f, elevation, p)

Computes annual cloud attenuation along a slant path based on Section 3. 
    
# Arguments
- `latlon`: object representing latitude and longitude, must be convertible to `ItuRPropagation.LatLon`
- `f`: frequency (GHz)
- `el`: elevation angle (degrees)
- `p`: exceedance probability (%)

# Return
- `Acloud::Real`: slant path cloud attenuation (dB)
"""
function cloudattenuation(latlon, f, el, p; warn=!SUPPRESS_WARNINGS[])
    L = liquidwatercontent(latlon, p)
    return cloudattenuation(latlon, f, el; L, warn)
end
function cloudattenuation(
    latlon,
    f,
    el;
    L,
    warn=!SUPPRESS_WARNINGS[]
)
    el = _torad(el)
    f = _toghz(f)
    deg2rad(5) ≤ el ≤ deg2rad(90) || !warn || @noinline(@warn("ItuR840.cloudattenuation only supports elevation angles between 5 and 90 degrees.\nThe given elevation angle $el degrees is outside this range so results may be inaccurate."))
    1 ≤ f ≤ 200 || !warn || @noinline(@warn("ItuR840.cloudattenuation only supports frequencies between 1 and 200 GHz.\nThe given frequency $f GHz is outside this range so results may be inaccurate."))
    
    # Section 3.2, equation 12
    (; K_L) = _K_L(f)

    # equation 13
    Acloud = L * K_L / sin(el)
    return Acloud
end

end # module ItuRP840
