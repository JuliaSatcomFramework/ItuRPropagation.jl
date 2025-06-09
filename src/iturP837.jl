module ItuRP837

#=
Rainfall rate statistics with a 1-min integration time are required for the prediction of rain attenuation
 in terrestrial links (e.g. Recommendation ITU-R P.530) and Earth-space links (e.g. Recommendation ITU-R P.618).

When reliable long-term local rainfall rate data is not available, Annex 1 of this Recommendation provides
 a rainfall rate prediction method for the prediction of rainfall rate statistics with a 1-min integration time
 This prediction method is based on: a) total monthly rainfall data generated from the GPCC Climatology (V 2015) database
 over land and from the European Centre for Medium-Range Weather Forecast (ECMWF) ERA Interim re-analysis database
 over water, and b) monthly mean surface temperature data in Recommendation ITU-R P.1510.

When reliable long-term local rainfall rate data is available with integration times greater than 1-min,
 Annex 2 of this Recommendation provides a method for converting rainfall rate statistics with integration
 times that exceed 1-min to rainfall rate statistics with a 1-min integration time.
=#

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion, tolatlon
using ..ItuRP1144: ItuRP1144, SquareGridData
using Artifacts

const version = ItuRVersion("ITU-R", "P.837", 7, "(06/2017)")

# Exports and constructor with separate latitude and longitude arguments
for name in (:rainfallrate001,)
    @eval $name(lat::Number, lon::Number, args...; kwargs...) = $name(LatLon(lat, lon), args...; kwargs...)
    @eval export $name
end

#region initialization

const δlat = 0.125
const δlon = 0.125
const latrange = range(-90, 90, step=δlat)
const lonrange = range(-180, 180, step=δlon)
const datasize = (length(latrange), length(lonrange))

const SGD_TYPE = let
    T = Float64
    R = typeof(latrange)
    SquareGridData{T, R, String}
end

@kwdef mutable struct RainfallRate001
    itp::Union{SGD_TYPE, Nothing} = nothing
end

const R001_DATA = RainfallRate001()

function initialize!()
    data = read!(joinpath(artifact"p837_R001", "R001.bin"), zeros(datasize))
    R001_DATA.itp = SquareGridData(latrange, lonrange, data, "Rainfall rate exceeded 0.01% of the year")
    return nothing
end

#endregion initialization

"""
    rainfallrate001(latlon)
    rainfallrate001(lat, lon)

Computes rainfall rate exceeded 0.01% via bi-linear interpolation as described in Annex 1.

# Arguments
- `latlon`: object representing latitude and longitude, must be convertible to `ItuRPropagation.LatLon`

# Return
- `R::Float64`: annual rainfall rate exceeded 0.01%
"""
function rainfallrate001(latlon)
    latlon = tolatlon(latlon)
    itp = @something(R001_DATA.itp, let
        initialize!()
        R001_DATA.itp
    end)::SGD_TYPE
    return itp(latlon)
end

end # module ItuRP837
