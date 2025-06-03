module ItuRP1144

export bilinear_interpolation

using ..ItuRPropagation: ItuRPropagation, LatLon

#= 
Bare computation of bilinear interpolation following section 1b of ITU-R P.1144-12.
This already expects the 4 surrounding values used for interpolation and normalized fractional row/column indices: δr = r - R, δc = c - C
This barebone implementation is used because at least in P2145 you need to scale the original surrounding values by using data from multiple matrices.
=#
function bilinear_interpolation(vals::NTuple{4, Real}, δr::Real, δc::Real)
    vals[1] * (1 - δr) * (1 - δc) + # R,C
    vals[2] * δr * (1 - δc) + # R+1,C
    vals[3] * (1 - δr) * δc + # R,C+1
    vals[4] * δr * δc # R+1,C+1
end

# This will perform bilinear interpolation of the points stored in `data` assuming the 4 neighboring indices are stored in `idxs` and expecting as input δr = r - R and δc = c - C, where r, R, c, C are the variables used in ITU-R P.1144-12
# function bilinear_interpolation(data::Matrix, idxs::NTuple{4, CartesianIndex}, δr::Real, δc::Real)
#     vals = ntuple(4) do i
#         data[idxs[i]]
#     end
#     bilinear_interpolation(vals, δr, δc)
# end

# Function that takes a LatLon position and the lat and lon ranges of the grid on which to perform bilinear interpolation and returns the 4 indices of the points on the grid surrounding the target position as well as the normalized fractional row/column indices δr and δc
@inline function bilinear_itp_inputs(latlon::LatLon, latrange::AbstractRange, lonrange::AbstractRange)
    R = searchsortedlast(latrange, latlon.lat)
    C = searchsortedlast(lonrange, latlon.lon)
    δlat = step(latrange)
    δlon = step(lonrange)
    R₊₁ = min(R + 1, length(latrange))
    C₊₁ = min(C + 1, length(lonrange))
    δr = (latlon.lat - latrange[R]) / δlat
    δc = (latlon.lon - lonrange[C]) / δlon
    idxs = (
        CartesianIndex(R, C),
        CartesianIndex(R₊₁, C),
        CartesianIndex(R, C₊₁),
        CartesianIndex(R₊₁, C₊₁),
    )
    (; idxs, δr, δc)
end


end