module ItuRP1144

export bilinear_interpolation

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion

const version = ItuRVersion("ITU-R", "P.1144", 12, "(08/2023)")

##### Bi-Linear, Section 1b #####

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

##### Bi-Cubic, Section 2 #####

# This function just provides the interpolation coefficints on a single line based on 4 values at the line points, and on the normalized fractional distance δ′ which is (c - C+1) for columns and (r - R+1) for rows following the notation of Section 2 of Annex 1 of ITU-R P.1144-12
function _Kvals(δ′::Real)
    a = -0.5
    K₁(absδ) = (a + 2) * absδ^3 - (a + 3) * absδ^2 + 1
    K₂(absδ) = a * absδ^3 - 5a * absδ^2 + 8a * absδ - 4a
    (
        K₂(δ′ + 1), # This is associated with j = C
        K₂(δ′),     # This is associated with j = C+1
        K₂(1 - δ′), # This is associated with j = C+2
        K₂(2 - δ′), # This is associated with j = C+3
    )
end

# This function assumes that the input is wrapped around so we don't have to check for reaching the end
function bicubic_itp_inputs(latlon::LatLon, latrange::AbstractRange, lonrange::AbstractRange)
    R₊₁ = searchsortedlast(latrange, latlon.lat)
    C₊₁ = searchsortedlast(lonrange, latlon.lon)
    δlat = step(latrange)
    δlon = step(lonrange)
    δr = (latlon.lat - latrange[R₊₁]) / δlat
    δc = (latlon.lon - lonrange[C₊₁]) / δlon
    idxs = ntuple(4) do Δr
        ntuple(4) do Δc
            CartesianIndex(R₊₁ + Δr - 2, C₊₁ + Δc - 2)
        end
    end
    (; idxs, δr, δc)
end

function bicubic_interpolation(data::Matrix, idxs::NTuple{4, NTuple{4, CartesianIndex}}, δr::Real, δc::Real)
    row_Kvals = _Kvals(δr)
    col_Kvals = _Kvals(δc)
    out = 0.0
    @inbounds for ridx in eachindex(idxs, row_Kvals)
        colidxs = idxs[ridx]
        for cidx in eachindex(colidxs, col_Kvals)
            idx = colidxs[cidx]
            R, C = Tuple(idx)
            val = data[idx]
            Kc = col_Kvals[cidx]
            Kr = row_Kvals[ridx]
            @info "asd" val Kc Kr R C
            out += val * Kc * Kr
        end
    end
    return out
end

function bicubic_interpolation(data::Matrix, latlon::LatLon, latrange::AbstractRange, lonrange::AbstractRange)
    R = searchsortedlast(latrange, latlon.lat) - 1
    C = searchsortedlast(lonrange, latlon.lon) - 1
    δr = (latlon.lat - latrange[R] + 1) / step(latrange)
    δc = (latlon.lon - lonrange[C] + 1) / step(lonrange)
    row_Kvals = _Kvals(δr)
    col_Kvals = _Kvals(δc)
    out = 0.0
    @inbounds for ri in 1:4
        for ci in 1:4
            val = data[R + ri - 1, C + ci - 1]
            Kc = col_Kvals[ci]
            Kr = row_Kvals[ri]
            out += val * Kc * Kr
        end
    end
    return out
end

end