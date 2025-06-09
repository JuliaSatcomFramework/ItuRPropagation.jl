module ItuRP839

#=
This Recommendation provides a method to predict the rain height for propagation prediction.
=#

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion, tolatlon, _tokm, _toghz, SUPPRESS_WARNINGS
using ..ItuRP1144: ItuRP1144, AbstractSquareGridITP, SquareGridData, SquareGridStatisticalData
using Artifacts
const version = ItuRVersion("ITU-R", "P.839", 4, "(09/2013)")

# Exports and constructor with separate latitude and longitude arguments
for name in (:isothermheight, :rainheightannual)
    @eval $name(lat::Number, lon::Number, args...; kwargs...) = $name(LatLon(lat, lon), args...; kwargs...)
    @eval export $name
end

#region initialization

const δlat = 1.5
const δlon = 1.5
const latrange = range(-90, 90, step=δlat)
const lonrange = range(-180, 180, step=δlon)
const datasize = (length(latrange), length(lonrange))

const H0_DATA = SquareGridData(latrange, lonrange, read!(joinpath(artifact"p839", "h0.bin"), zeros(datasize)), "Isotherm height")

#endregion initialization

"""
    rainheightannual(latlon::LatLon)

Computes rain height based on the equation in Section 2.
h0 will be interpolated for the given latitude and longitude.

# Arguments
- `latlon`: object representing latitude and longitude, must be convertible to `ItuRPropagation.LatLon`

# Return
- `hR`: mean annual rain height (km) above mean sea level
"""
function rainheightannual(latlon)
    h0 = isothermheight(latlon)

    hR = h0 + 0.36  # equation in section 2
    return hR
end


"""
    isothermheight(latlon)

Calculates isothermic height based on bilinear interpolation.
h0 will be interpolated for the given latitude and longitude.

# Arguments
- `latlon`: object representing latitude and longitude, must be convertible to `ItuRPropagation.LatLon`

# Return
- `h0`: mean annual 0°C isotherm height (km) above mean sea level
"""
function isothermheight(latlon)
    latlon = tolatlon(latlon)
    return H0_DATA(latlon)
end

end # module ItuRP839
