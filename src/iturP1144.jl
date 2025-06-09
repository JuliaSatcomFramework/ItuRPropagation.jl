module ItuRP1144

export bilinear_interpolation

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion, tolatlon

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

@inline function _checkp(p, (pmin, pmax); warn::Bool, kind::String)
    _p = min(max(p, pmin), pmax)
    pmin ≤ p ≤ pmax || !warn || @noinline(@warn "The provided value p = $(p)% is not supported for the interpolation of $kind, results are given for p = $(_p)%")
    return _p
end

### Helper for extracting upper and lower indices for interpolation of a single real variable, mostly used for the exceedance probability
@inline function ccdf_itp_inputs(p::Real, pvec::AbstractVector; warn::Bool, kind::String)
    pmin = first(pvec)
    pmax = last(pvec)
    p = _checkp(p, (pmin, pmax); warn, kind)
    prange = searchsorted(pvec, p)
    pindexbelow = prange.stop
    pindexabove = prange.start

    (; pindexabove, pindexbelow)
end

### Structure for holding data for simple square grid interpolation
abstract type AbstractSquareGridITP end
@kwdef struct SquareGridData{T, R, NT} <: AbstractSquareGridITP
    latrange::R
    lonrange::R
    data::Matrix{T} # Stores the actual data
    id::NT # This stores any extra information
end

function (sd::SquareGridData)(idxs::NTuple{4, CartesianIndex}, δr::Real, δc::Real) 
    vals = ntuple(4) do i
        sd.data[idxs[i]]
    end
    bilinear_interpolation(vals, δr, δc)
end
(sd::SquareGridData)(latlon::LatLon) = sd(bilinear_itp_inputs(latlon, sd.latrange, sd.lonrange)...)

struct SquareGridStatisticalData{SGD}
    ps::Vector{Float64} # Stores the exceedance probabilities
    items::Vector{SGD} # Stores the data
end

# This implements the square interpolation for statistical attenuations, used for P840 and P434
function (sg::SquareGridStatisticalData)(latlon::LatLon, p; warn::Bool, kind::String, kwargs...)
    fd = first(sg.items)
    (; idxs, δr, δc) = bilinear_itp_inputs(latlon, fd.latrange, fd.lonrange)
    (; pindexabove, pindexbelow) = ccdf_itp_inputs(p, sg.ps; warn, kind)
    
    above = sg.items[pindexabove](idxs, δr, δc; kwargs...)
    pindexabove == pindexbelow && return above
    below = sg.items[pindexbelow](idxs, δr, δc; kwargs...)
    psabove = sg.ps[pindexabove]
    psbelow = sg.ps[pindexbelow]
    out = (above - below) / log(psabove/psbelow) * log(p/psbelow) + below
    return out
end


##### Bi-Cubic, Section 2 #####

# This function just provides the interpolation coefficints on a single line based on 4 values at the line points, and on the normalized fractional distance δ′ which is (c - C+1) for columns and (r - R+1) for rows following the notation of Section 2 of Annex 1 of ITU-R P.1144-12
function _Kvals(δ′::Real)
    a = -0.5
    K₁(absδ) = (a + 2) * absδ^3 - (a + 3) * absδ^2 + 1
    K₂(absδ) = a * absδ^3 - 5a * absδ^2 + 8a * absδ - 4a
    (
        K₂(δ′ + 1), # This is associated with j = C
        K₁(δ′),     # This is associated with j = C+1
        K₁(1 - δ′), # This is associated with j = C+2
        K₂(2 - δ′), # This is associated with j = C+3
    )
end

"""
    bicubic_interpolation(data::Matrix, latlon::LatLon, latrange::AbstractRange, lonrange::AbstractRange)

Performs bicubic interpolation at point `latlon` from the gridded values in `data` assuming the grid to be identified by the ranges `latrange` and `lonrange`.

**Note:** This function assumes that the input grid is already including 2 extra rows/columns at the edges to simplify the interpolation.
"""
function bicubic_interpolation(data::Matrix, latlon::LatLon, latrange::AbstractRange, lonrange::AbstractRange)
    R₊₁ = searchsortedlast(latrange, latlon.lat)
    C₊₁ = searchsortedlast(lonrange, latlon.lon)
    δlat = step(latrange)
    δlon = step(lonrange)
    δr = (latlon.lat - latrange[R₊₁]) / δlat
    δc = (latlon.lon - lonrange[C₊₁]) / δlon
    Krows = _Kvals(δr)
    Kcols = _Kvals(δc)
    out = 0.0
    for ridx in eachindex(Krows)
        for cidx in eachindex(Kcols)
            out += data[R₊₁ + ridx - 2, C₊₁ + cidx - 2] * Krows[ridx] * Kcols[cidx]
        end
    end
    out
end

end