module ItuRP453

#=
Recommendation ITU-R P.453 provides methods to estimate the radio refractive index 
 and its behaviour for locations worldwide; describes both surface and vertical profile
 characteristics; and provides global maps for the distribution of refractivity
 parameters and their statistical variation.
=#

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion, _tolatlon, _tokm, IturEnum, EnumWater, EnumIce, SUPPRESS_WARNINGS
using ..ItuRP1144: ItuRP1144, AbstractSquareGridITP, SquareGridData, SquareGridStatisticalData
using Artifacts
const version = ItuRVersion("ITU-R", "P.453", 14, "(08/2019)")

# Exports and constructor with separate latitude and longitude arguments
for name in (:wettermsurfacerefractivityannual, :wettermsurfacerefractivityannual_50)
    @eval $name(lat::Number, lon::Number, args...; kwargs...) = $name(LatLon(lat, lon), args...; kwargs...)
    @eval export $name
end

#region initialization

const δlat = 0.75
const δlon = 0.75
const latrange = range(-90, 90, step=δlat)
const lonrange = range(-180, 180, step=δlon)
const datasize = (length(latrange), length(lonrange))

# exceedance probabilities
const psannual = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 60, 70, 80, 90, 95, 99]

# exceedance probability values for reading files
const filespsannual = ["01", "02", "03", "05", "1", "2", "3", "5", "10", "20", "30", "50", "60", "70", "80", "90", "95", "99"]

const SGD_TYPE = let
    T = Float64
    R = typeof(latrange)
    SquareGridData{T, R, String}
end

# This will hold the annual ccdf data of wet term of surface refractivity for computation of the cloud attenuation. The underlying data is already a interpolator that will use bilinear interpolation to find the value at the given location
@kwdef mutable struct WetTermData
    ccdf::Union{Nothing, SquareGridStatisticalData{SGD_TYPE}} = nothing
end

const NWET_ANNUAL = WetTermData()

function initialize!()
    isnothing(NWET_ANNUAL.ccdf) || return nothing
    @info "P453: Loading annual data of wet term of surface refractivity"
    items = map(zip(psannual, filespsannual)) do (p, suffix)
        data = read!(joinpath(artifact"p453_nwet_annual", "NWET_Annual_$(suffix).bin"), zeros(datasize))
        name = "Wet term of surface refractivity at $p% exceedance probability"
        SquareGridData(latrange, lonrange, data, name)
    end
    NWET_ANNUAL.ccdf = SquareGridStatisticalData(psannual, items)
    return nothing
end


#endregion initialization

#region internal functions

"""
    _saturationvapourpressure(t::Real, P::Real, ::Val{EnumWater})

Compute the saturation water vapour pressure for water-type hydrometer based on Section 1.

# Arguments
- `T::Real`: temperature in (°C)
- `P::Real`: total atmospheric pressure (hPa)
- `::Val{EnumWater}`: water-type hydrometer

# Return
- `eₛ::Real`: saturation water vapour pressure (hPa)
"""
function _saturationwatervapourpressure(T::Real, P::Real, ::Val{EnumWater})
    (T < -40 || T > 50) && throw(ArgumentError("T=$T\nT (temperature in deg C) in ItuRP453.saturationvapourpressure must be between -80 and 0."))
    EF = 1 + 1e-4 * (7.2 + P * (0.0320 + 5.9e-6 * T * T))
    a, b, c, d = 6.1121, 18.678, 257.14, 234.5
    eₛ = EF * a * exp((b - T / d) * T / (T + c))     # equation 9
    return eₛ
end


"""
    _saturationvapourpressure(t::Real, P::Real, ::Val{EnumIce})

Compute the saturation water vapour pressure for ice-type hydrometer based on Section 1.

# Arguments
- `T::Real`: temperature (°C)
- `P::Real`: total atmospheric pressure (hPa)
- `::Val{EnumIce}`: ice-type hydrometer

# Return
- `eₛ::Real`: saturation water vapour pressure (hPa)
"""
function _saturationwatervapourpressure(T::Real, P::Real, ::Val{EnumIce})
    (T < -80 || T > 0) && throw(ArgumentError("T=$T\nT (temperature in deg C) in ItuRP453.saturationvapourpressure must be between -80 and 0."))
    EF = 1 + 1e-4 * (2.2 + P * (0.0383 + 6.4e-6 * T * T))
    a, b, c, d = 6.1115, 23.036, 279.82, 333.7
    eₛ = EF * a * exp((b - T / d) * T / (T + c))     # equation 9
    return eₛ
end

#endregion internal functions

"""
    radiorefractiveindex(Pd::Real, T::Real, e::Real)

Compute the atmospheric radio refractive index \$\\sqrt[n]{1 + x + x^2 + \\ldots}\$ based on Section 1.

# Arguments
- `T::Real`: absolute temperature (°K)
- `Pd::Real`: dry atmospheric pressure (hPa)
- `e::Real`: saturation water vapour pressure (hPa)

# Return
- `n::Real`: atmospheric radio refractive index
"""
function radiorefractiveindex(
    T, 
    Pd, 
    e
)
    Ndry = drytermradiorefractivity(T, Pd)
    Nwet = wettermradiorefractivity(T, e)
    N = Ndry + Nwet      # equation 2
    n = 1 + N * (1e-6)     # equation 1
    return n
end


"""
    drytermradiorefractivity(Pd, T)

Compute the dry term of the radio refractivity based on Section 1.

# Arguments
- `T`: absolute temperature (°K)
- `Pd`: dry atmospheric pressure (hPa)

# Return
- `Ndry`: dry term of the radio refractivity
"""
function drytermradiorefractivity(
    T, 
    Pd
)
    Ndry = 77.6 * (Pd / T)     # equation 3
    return Ndry
end


"""
    wettermradiorefractivity(Pd::Real, T::Real)

Compute the wet term of the radio refractivity based on Section 1.

# Arguments
- `T`: absolute temperature (°K)
- `e`: water vapour pressure (hPa)

# Return
- `Nwet`: wet term of the radio refractivity
"""
function wettermradiorefractivity(
    T, 
    e
)
    Nwet = 72 * (e / T) + 3.75e5 * (e / (T * T))     # equation 4
    return Nwet
end


"""
    vapourpressure(T::Real, P::Real, H::Real, hydrometer::IturEnum)

Computes the vapour pressure dependent on hydrometer type based on Section 1.

# Arguments
- `T::Real`: temperature in (°C)
- `P::Real`: total atmospheric pressure (hPa)
- `H::Real`: relative humidity (%)
- `hydrometer::IturEnum"`: type of hydrometer, either water or ice

# Return
- `e::Real`: water vapour pressure (hPa)
"""
function vapourpressure(
    T::Real, 
    P::Real, 
    H::Real,
    hydrometer::IturEnum
)
    if hydrometer == EnumIce
        eₛ = _saturationwatervapourpressure(T, P, EnumIce)
    elseif hydrometer == EnumWater
        eₛ = _saturationwatervapourpressure(T, P, EnumWater)
    else
        throw(ArgumentError("Invalid hydrometer value in ItuRP453.vapourpressure."))
    end
    e = H * eₛ / 100.0     # equation 8
    return e
end

"""
    wettermsurfacerefractivityannual(latlon, p)

Interpolates wet term of the surface refractivity at an exceedance probability `p` provided as input in % based on Section 2.2.

# Arguments
- `latlon`: object representing latitude and longitude, must be convertible to `ItuRPropagation.LatLon`
- `p`: exceedance probability (%)   

# Return
- `Nwet::Real`: wet term of the surface refractivity (ppm)
"""
function wettermsurfacerefractivityannual(latlon, p; warn=!SUPPRESS_WARNINGS[])
    itp = @something(NWET_ANNUAL.ccdf,let
        initialize!()
        NWET_ANNUAL.ccdf
    end)::SquareGridStatisticalData{SGD_TYPE}
    return itp(latlon, p; warn, kind="the wet term of surface refractivity")
end

"""
    wettermsurfacerefractivityannual_50(latlon)

Returns the wet term of surface refractivity at 50% exceedance probability (Based on annual data) at the desired location
"""
function wettermsurfacerefractivityannual_50(latlon)
    ccdf = @something(NWET_ANNUAL.ccdf,let
        initialize!()
        NWET_ANNUAL.ccdf
    end)::SquareGridStatisticalData{SGD_TYPE}
    itp =  ccdf.items[12] # index 12 is the one corresponding to 50% exceedance probability
    return itp(latlon)
end

end # module ItuRP453
