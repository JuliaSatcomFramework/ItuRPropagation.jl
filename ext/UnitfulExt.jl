module UnitfulExt

using ItuRPropagation: ItuRPropagation, _todeg, _toghz, _tokm, LatLon
using Unitful

const Len = Quantity{<:Real,u"𝐋"}
const Freq = Quantity{<:Real,u"𝐓^-1"}
const Deg = Quantity{<:Real,NoDims,typeof(u"°")}
const Rad = Quantity{<:Real,NoDims,typeof(u"rad")}
const Angle = Union{Deg,Rad}

ItuRPropagation.LatLon(lat::Angle, lon::Angle) = LatLon(_todeg(lat), _todeg(lon))
ItuRPropagation.LatLon(lat::Real, lon::Angle) = LatLon(lat, _todeg(lon))
ItuRPropagation.LatLon(lat::Angle, lon::Real) = LatLon(_todeg(lat), lon)

@inline ItuRPropagation._todeg(val::Angle) = uconvert(u"°", val) |> ustrip
@inline ItuRPropagation._toghz(val::Freq) = uconvert(u"GHz", val) |> ustrip
@inline ItuRPropagation._tokm(val::Len) = uconvert(u"km", val) |> ustrip

end