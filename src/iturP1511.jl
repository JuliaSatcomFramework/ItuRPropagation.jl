module ItuRP1511

#=
This Recommendation provides global topographical data, information on geographic coordinates, and
height data for the prediction of propagation effects for Earth-space paths in ITU-R recommendations.
=#

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion, tolatlon, _tokm, ItuRP1144
using Artifacts
const version = ItuRVersion("ITU-R", "P.1511", 3, "(08/2023)")

# Exports and constructor with separate latitude and longitude arguments
for name in (:topographicheight,)
    @eval $name(lat::Number, lon::Number, args...; kwargs...) = $name(LatLon(lat, lon), args...; kwargs...)
    @eval export $name
end

#region initialization

const GRID_DATA = (;
    topo = let 
        latrange = range(-90.125, 90.125, step=1/12)
        lonrange = range(-180.125, 180.125, step=1/12)
        matsize = (length(latrange), length(lonrange))
        data = read!(artifact"p1511/TOPO.bin", zeros(matsize))
        (; latrange, lonrange, data)
    end,
    # For the moment we skip the geoid undulation as we don't use it
    # egm = let
    #     latrange = range(-90 - 2/12, 90 + 2/12, step=1/12)
    #     lonrange = range(-180 - 2/12, 180 + 2/12, step=1/12)
    #     matsize = (length(latrange), length(lonrange))
    #     data = read!(artifact"p1511/EGM2008.bin", zeros(matsize))
    #     (; latrange, lonrange, data)
    # end,
)
#endregion initialization

"""
    topographicheight(latlon::LatLon)
    topographicheight(lat, lon)

Calculates topographic height as per Section 1.1 of ITU-R P.1511-3.

# Arguments
- `latlon::LatLon`: latitude and longitude (degrees)

# Return
- `I::Real`: height (km)
"""
function topographicheight(latlon)
    latlon = tolatlon(latlon)
    grid_data = GRID_DATA.topo
    alt = ItuRP1144.bicubic_interpolation(grid_data.data, latlon, grid_data.latrange, grid_data.lonrange) / 1e3
    return alt
end

end # module ItuRP1511
