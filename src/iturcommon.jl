export LatLon
export IturEnum
export EnumWater, EnumIce
export EnumHeightAndIndex, EnumIndexOnly
export EnumHorizontalPolarization, EnumVerticalPolarization, EnumCircularPolarization
export ItuRVersion

struct LatLon
    lat::Float64
    lon::Float64
    function LatLon(lat, lon)
        if lat < -90.0 || lat > 90.0
            throw(ErrorException("lat=$lat\nlat (latitude) should be between -90 and 90"))
        end
        if lon > 180
            lon -= 360
        end
        if lon < -180 || lon > 180
            throw(ErrorException("lon=$lat\nlon (longitude) should be between -180 and 180"))
        end
        new(lat, lon)
    end
end

Base.show(io::IO, p::LatLon) = print(io, "(", p.lat, ", ", p.lon, ")")

@enum IturEnum begin
    # for ITU-R P.453-14
    EnumWater
    EnumIce

    # for ITU-R P.1511-11 (custom by Kchao)
    EnumHeightAndIndex
    EnumIndexOnly

    # for ITU-R P.838-3
    EnumHorizontalPolarization
    EnumVerticalPolarization
    EnumCircularPolarization
end
function tilt_from_polarization(polarization::IturEnum)
    τ = if polarization == EnumCircularPolarization
        45
    elseif polarization == EnumHorizontalPolarization
        0
    elseif polarization == EnumVerticalPolarization
        90
    else
        throw(ArgumentError("Invalid polarization value in ItuRP838.rainspecificattenuation."))
    end
    return τ
end

struct ItuRVersion
    doctype::String
    recommendation::String
    dashednumber::Int8
    datestring::String
end

"""
    SquareGridInterpolator

Creates an interpolator on a square grid following the bilinear interpolation
method specified in Section 1b of Annex 1 of ITU-R P.1144-12.
"""
struct SquareGridInterpolator{T,R,C}
    rowrange::R
    colrange::C
    data::Matrix{T}
end

SquareGridInterpolator(latrange, lonrange) = SquareGridInterpolator(latrange, lonrange, fill(NaN, length(latrange), length(lonrange)))

@inline function (itp::SquareGridInterpolator)(R::Integer, C::Integer, r, c)
    (; data ) = itp
    δr = r - R
    δc = c - C
    latidxs, lonidxs = axes(data)
    R₊₁ = min(R+1, last(latidxs))
    C₊₁ = min(C+1, last(lonidxs))
    DL, UL, DR, UR = data[R, C], data[R₊₁, C], data[R, C₊₁], data[R₊₁, C₊₁]
    out = (
        DL * (1 - δr) * (1 - δc) +
        UL * δr * (1 - δc) + 
        DR  * (1 - δr) * δc + 
        UR * δr * δc
    )
    return out
end

function (itp::SquareGridInterpolator)(rowval, colval)
    R = searchsortedlast(itp.rowrange, rowval)
    C = searchsortedlast(itp.colrange, colval)
    Δrow = R === length(itp.rowrange) ? 1 : itp.rowrange[R + 1] - itp.rowrange[R]
    Δcol = C === length(itp.colrange) ? 1 : itp.colrange[C + 1] - itp.colrange[C]
    r = R + (rowval - itp.rowrange[R]) / Δrow
    c = C + (colval - itp.colrange[C]) / Δcol
    val = itp(R, C, r, c)
    return val
end

struct AltitudeScaledSquareGridInterpolator{T,R,C,F}
    rowrange::R
    colrange::C
    data::Matrix{T}
    scale::Matrix{T}
    alt::Matrix{T}
    f::F
end

@inline function (itp::AltitudeScaledSquareGridInterpolator)(R::Integer, C::Integer, r, c; alt)
    (; data ) = itp
    δr = r - R
    δc = c - C
    latidxs, lonidxs = axes(data)

    R₊₁ = min(R+1, last(latidxs))
    C₊₁ = min(C+1, last(lonidxs))
    iters = ( # idx, Weight
        (CartesianIndex(R, C), (1 - δr) * (1 - δc)),
        (CartesianIndex(R₊₁, C), δr * (1 - δc)),
        (CartesianIndex(R, C₊₁), (1 - δr) * δc),
        (CartesianIndex(R₊₁, C₊₁), δr * δc)
    )
    # Function to compute Xᵢ
    fun(i) = let itp=itp, alt = alt
        Xᵢ′ = itp.data[i]
        scaleᵢ = itp.scale[i]
        altᵢ = itp.alt[i]
        itp.f(Xᵢ′, alt, scaleᵢ, altᵢ)
    end
    out = 0.0
    @inbounds @simd for el in iters
        idx, weight = el
        out += fun(idx) * weight
    end
    return out
end