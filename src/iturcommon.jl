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

# This functions are used for processing inputs and allow easier extension to support types from other packages
@inline _torad(val::Real) = deg2rad(val)

# We expects length values to be provided in km
@inline _tokm(val::Real) = val

# We expects inputs in GHz
@inline _toghz(val::Real) = val

@inline _tolatlon(x) = convert(LatLon, x)::LatLon # Type assertion on convert is used to help compiler, see https://github.com/JuliaLang/julia/issues/42372 for more details

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