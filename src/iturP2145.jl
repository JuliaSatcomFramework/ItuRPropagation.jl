module ItuRP2145

#=
This Recommendation provides methods to predict the surface total (barometric) pressure, surface
temperature, surface water vapour density and integrated water vapour content required for the calculation of
gaseous attenuation and related effects on terrestrial and Earth-space paths.
=#

using ..ItuRPropagation: ItuRPropagation, LatLon, ItuRVersion, SquareGridInterpolator
using Artifacts: Artifacts, @artifact_str
using Interpolations: Interpolations, linear_interpolation
using DelimitedFiles: DelimitedFiles, readdlm

const version = ItuRVersion("ITU-R", "P.2145", 0, "(08/2022)")

#region initialization

@kwdef mutable struct Initialized
    T::Bool = false
    RHO::Bool = false
    V::Bool = false
    P::Bool = false
    Z_ground::Bool = false
end
const INITIALIZED = Initialized()


const δlat = 0.25
const δlon = 0.25
const latrange = range(-90, 90, step=δlat)
const lonrange = range(-180, 180, step=δlon)
const datasize = (length(latrange), length(lonrange))

const latsize = length(latrange) + 1 # number of latitude points (-90, 90, 0.25) plus one extra row for interpolation
const lonsize = length(lonrange) + 1 # number of longitude points (-180, 180, 0.25) plus one extra column for interpolation

# exceedance probability, section 1 of ITU-R P.836-6
const psannual = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 60, 70, 80, 90, 95, 99]
const npsannual = length(psannual)

# exceedance probability values for reading files
const filespsannual = ["001", "002", "003", "005", "01", "02", "03", "05", "1", "2", "3", "5", "10", "20", "30", "50", "60", "70", "80", "90", "95", "99"]

const AnnualData = let 
    Z_ground = fill(NaN, datasize) # We use one single matrix for all data as altitude is consistent
    f_T(Xᵢ′, scaleᵢ, altᵢ; alt) = Xᵢ′ + scaleᵢ * (alt - altᵢ) # This is the scaling function for temperature
    f_others(Xᵢ′, scaleᵢ, altᵢ; alt) = Xᵢ′ * exp((alt - altᵢ) / scaleᵢ) # This is the scaling function for P, V and RHO
    # This will be a NamedTuple
    (;
        T = (;
            ccdf = map(_ -> fill(NaN, datasize), filespsannual),
            mean = fill(NaN, datasize),
            scale = fill(NaN, datasize), # TSCH.TXT
            Z_ground,
            scale_func = f_T,
        ),
        RHO = (;
            ccdf = map(_ -> fill(NaN, datasize), filespsannual),
            mean = fill(NaN, datasize),
            scale = fill(NaN, datasize), # VSCH.TXT
            Z_ground,
            scale_func = f_others,
        ),
        V = (;
            ccdf = map(_ -> fill(NaN, datasize), filespsannual),
            mean = fill(NaN, datasize),
            scale = fill(NaN, datasize), # VSCH.TXT
            Z_ground,
            scale_func = f_others,
        ),
        P = (;
            ccdf = map(_ -> fill(NaN, datasize), filespsannual),
            mean = fill(NaN, datasize),
            scale = fill(NaN, datasize), # PSCH.TXT
            Z_ground,
            scale_func = f_others,
        ),
    )
end

is_initialized(data::Matrix) = data |> first |> !isnan

function initialize!(nt, kind::Symbol)
    getproperty(INITIALIZED, kind) && return nothing
    # We make sure that Z_ground is also initialized
    @info "P2145: Loading data for $kind variable, this may take some time..."
    if !INITIALIZED.Z_ground
        initialize!(AnnualData.T.Z_ground, :Z_ground, "Z_ground.TXT")
        INITIALIZED.Z_ground = true
    end
    scalename = if kind === :T
        "TSCH.TXT"
    elseif kind === :P
        "PSCH.TXT"
    else 
        "VSCH.TXT"
    end
    initialize!(nt.scale, kind, scalename)
    # We initialize the mean
    initialize!(nt.mean, kind, "$(kind)_mean.TXT")
    # We initialize the ccdf values
    for (i, suffix) in enumerate(filespsannual)
        initialize!(nt.ccdf[i], kind, "$(kind)_$(suffix).TXT")
    end
    setproperty!(INITIALIZED, kind, true)
    return nothing
end

function initialize!(data::Matrix, kind::Symbol, filename::String)
    is_initialized(data) && return nothing
    dir = if kind in (:T, :Z_ground)
        joinpath(artifact"p2145_annual", "T_Annual")
    elseif kind === :RHO
        joinpath(artifact"p2145_annual", "RHO_Annual")
    elseif kind === :V
        joinpath(artifact"p2145_annual", "V_Annual")
    elseif kind === :P
        joinpath(artifact"p2145_annual", "P_Annual")
    end
    data .= readdlm(joinpath(dir, filename), ' ')
    return nothing
end

# This is a helper function to return indices for faster interpolation inside SquareGridInterpolator
@inline function itp_inputs(latlon::LatLon)
    R = searchsortedlast(latrange, latlon.lat)
    C = searchsortedlast(lonrange, latlon.lon)
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
@inline function itp_inputs(p::Real)
    prange = searchsorted(psannual, p)
    pindexbelow = prange.stop
    pindexabove = prange.start

    (; pindexabove, pindexbelow)
end

# This will perform bilinear interpolation of the points stored in `data` assuming the 4 neighboring indices are stored in `idxs` and expecting as input δr = r - R and δc = c - C, where r, R, c, C are the variables used in ITU-R P.1144-12
function bilinear_interpolation(data::Matrix, idxs::NTuple{4, CartesianIndex}, δr::Real, δc::Real)
    vals = ntuple(4) do i
        data[idxs[i]]
    end
    bilinear_interpolation(vals, δr, δc)
end

# This will compute first altitude-based scaling using function `f` (Following ITU-R P.2145 guidelines) to each of the 4 neighboring points and then perform bilinear interpolation as per ITU-R P.1144-12
function bilinear_interpolation(data::Matrix, scale::Matrix, Z::Matrix, f::F, idxs::NTuple{4, CartesianIndex}, δr::Real, δc::Real; alt = 0.0) where F
    vals = ntuple(4) do i
        idx = idxs[i]
        Xᵢ′ = data[idx]
        altᵢ = Z[idx]
        scaleᵢ = scale[idx]
        f(Xᵢ′, scaleᵢ, altᵢ; alt)
    end
    bilinear_interpolation(vals, δr, δc)
end

function bilinear_interpolation(vals::NTuple{4, Real}, δr::Real, δc::Real)
    vals[1] * (1 - δr) * (1 - δc) + # R,C
    vals[2] * δr * (1 - δc) + # R+1,C
    vals[3] * (1 - δr) * δc + # R,C+1
    vals[4] * δr * δc # R+1,C+1
end


const latvalues = [(-90.0 + (i - 1) * 0.25) for i in 1:latsize]
const lonvalues = [(-180.0 + (j - 1) * 0.25) for j in 1:lonsize]

const surfacetemperatureannualdata = zeros(Float64, (npsannual, latsize, lonsize))
const surfacerhodata = zeros(Float64, (npsannual, latsize, lonsize))
const scaleheightrhodata = zeros(Float64, (latsize, lonsize))
const surfaceheightdata = zeros(Float64, (latsize, lonsize))

const initialized = Ref{Bool}(false)

function initialize()
    initialized[] && return nothing
    tempdata = zeros(Float64, (latsize, lonsize))
    for nps in range(1, npsannual)
        read!(
            joinpath(artifact"input-maps", "surfacetemperatureannual_$(string(latsize))_x_$(string(lonsize))_x_$(filespsannual[nps]).bin"),
            tempdata
        )
        @views surfacetemperatureannualdata[nps, :, :] = tempdata
        read!(
            joinpath(artifact"input-maps", "surfacewatervapordensityannual_$(string(latsize))_x_$(string(lonsize))_x_$(filespsannual[nps]).bin"),
            tempdata
        )
        @views surfacerhodata[nps, :, :] = tempdata
    end
    read!(
        joinpath(artifact"input-maps", "scaleheightwatervapordensityannual_$(string(latsize))_x_$(string(lonsize)).bin"),
        scaleheightrhodata
    )

    read!(
        joinpath(artifact"input-maps", "surfaceheightannual_$(string(latsize))_x_$(string(lonsize)).bin"),
        surfaceheightdata
    )
    initialized[] = true
    return nothing
end

#endregion initialization

"""
    surfacetemperatureannual(latlon::LatLon, p::Real, hs::Union{Missing, Real} = missing)

Computes annual surface temperature for a given exceedance probability and altitude based on Section 2.1.

# Arguments
- `latlon::LatLon`: latitude and longitude (degrees)
- `p::Real`: exceedance probability (%)
- `hs::Union{Missing, Real}`: altitude (m), defaults to missing

# Return
- `T::Real`: annual surface temperature (°K)
"""
function surfacetemperatureannual(
    latlon::LatLon,
    p::Real,
)
    initialize()
    prange = searchsorted(psannual, p)
    pindexbelow = prange.stop
    pindexabove = prange.start
    pexact = pindexbelow == pindexabove ? true : false

    latrange = searchsorted(latvalues, latlon.lat)
    lonrange = searchsorted(lonvalues, latlon.lon)
    R = latrange.stop
    C = lonrange.stop

    δg = 0.25
    r = ((90.0 + latlon.lat) / δg) + 1
    c = ((180.0 + latlon.lon) / δg) + 1

    T00a = surfacetemperatureannualdata[pindexabove, R, C]
    T01a = surfacetemperatureannualdata[pindexabove, R, C+1]
    T10a = surfacetemperatureannualdata[pindexabove, R+1, C]
    T11a = surfacetemperatureannualdata[pindexabove, R+1, C+1]

    T00b = surfacetemperatureannualdata[pindexbelow, R, C]
    T01b = surfacetemperatureannualdata[pindexbelow, R, C+1]
    T10b = surfacetemperatureannualdata[pindexbelow, R+1, C]
    T11b = surfacetemperatureannualdata[pindexbelow, R+1, C+1]

    Tabove = (
        T00a * ((R + 1 - r) * (C + 1 - c)) +
        T10a * ((r - R) * (C + 1 - c)) +
        T01a * ((R + 1 - r) * (c - C)) +
        T11a * ((r - R) * (c - C))
    )

    if pexact == true
        return Tabove
    else
        Tbelow = (
            T00b * ((R + 1 - r) * (C + 1 - c)) +
            T10b * ((r - R) * (C + 1 - c)) +
            T01b * ((R + 1 - r) * (c - C)) +
            T11b * ((r - R) * (c - C))
        )
        pslogabove = log(psannual[pindexabove])
        pslogbelow = log(psannual[pindexbelow])
        T = (Tabove - Tbelow) / (pslogabove - pslogbelow) * (log(p) - pslogbelow) + Tbelow
        return T
    end
end

function surfacetemperatureannual2(latlon::LatLon; alt = 0.0)
    nt = AnnualData.T
    INITIALIZED.T || initialize!(nt, :T)
    data = nt.mean
    scale = nt.scale
    Z = nt.Z_ground
    f = nt.scale_func
    (; idxs, δr, δc) = itp_inputs(latlon)
    bilinear_interpolation(data, scale, Z, f, idxs, δr, δc; alt)
end

function surfacetemperatureannual2(latlon::LatLon, p::Real; alt = 0.0)
    (; idxs, δr, δc) = itp_inputs(latlon)
    (; pindexabove, pindexbelow) = itp_inputs(p)

    nt = AnnualData.T
    INITIALIZED.T || initialize!(nt, :T)

    Ts = nt.ccdf

    scale = nt.scale
    Z = nt.Z_ground
    f = nt.scale_func
    Tabove = bilinear_interpolation(Ts[pindexabove], scale, Z, f, idxs, δr, δc; alt)

    pindexabove == pindexbelow && return Tabove
    Tbelow = bilinear_interpolation(Ts[pindexbelow], scale, Z, f, idxs, δr, δc; alt)
    psabove = psannual[pindexabove]
    psbelow = psannual[pindexbelow]
    T = (Tabove - Tbelow) / log(psabove/psbelow) * log(p/psbelow) + Tbelow
    return T
end


"""
    surfacewatervapordensityannual(latlon::LatLon, p::Real, hs::Union{Missing, Real} = missing)

Computes annual surface water vapor density for a given exceedance probability and altitude based on Section 2.1.

# Arguments
- `latlon::LatLon`: latitude and longitude (degrees)
- `p::Real`: exceedance probability (%)
- `hs::Union{Missing, Real}`: altitude (m), defaults to missing

# Return
- `ρ::Real`: annual surface water vapor density (g/m^3)
"""
function surfacewatervapordensityannual(
    latlon::LatLon,
    p::Real,
    hs::Union{Missing,Real}=missing
)
    initialize()
    prange = searchsorted(psannual, p)
    pindexbelow = prange.stop
    pindexabove = prange.start
    pexact = pindexbelow == pindexabove ? true : false

    latrange = searchsorted(latvalues, latlon.lat)
    lonrange = searchsorted(lonvalues, latlon.lon)
    R = latrange.stop
    C = lonrange.stop

    δg = 0.25
    r = ((90.0 + latlon.lat) / δg) + 1
    c = ((180.0 + latlon.lon) / δg) + 1

    if hs === missing
        hs = ItuRP1511.topographicheight(latlon)
    end

    h00 = surfaceheightdata[R, C]
    h01 = surfaceheightdata[R, C + 1]
    h10 = surfaceheightdata[R + 1, C]
    h11 = surfaceheightdata[R + 1, C + 1]

    ρ′00a = surfacerhodata[pindexabove, R, C]
    ρ′01a = surfacerhodata[pindexabove, R, C+1]
    ρ′10a = surfacerhodata[pindexabove, R+1, C]
    ρ′11a = surfacerhodata[pindexabove, R+1, C+1]

    ρ′00b = surfacerhodata[pindexbelow, R, C]
    ρ′01b = surfacerhodata[pindexbelow, R, C+1]
    ρ′10b = surfacerhodata[pindexbelow, R+1, C]
    ρ′11b = surfacerhodata[pindexbelow, R+1, C+1]

    ρ00a = ρ′00a * exp((h00 - hs) / scaleheightrhodata[R, C])
    ρ01a = ρ′01a * exp((h01 - hs) / scaleheightrhodata[R, C+1])
    ρ10a = ρ′10a * exp((h10 - hs) / scaleheightrhodata[R+1, C])
    ρ11a = ρ′11a * exp((h11 - hs) / scaleheightrhodata[R+1, C+1])

    ρ00b = ρ′00b * exp((h00 - hs) / scaleheightrhodata[R, C])
    ρ01b = ρ′01b * exp((h01 - hs) / scaleheightrhodata[R, C+1])
    ρ10b = ρ′10b * exp((h10 - hs) / scaleheightrhodata[R+1, C])
    ρ11b = ρ′11b * exp((h11 - hs) / scaleheightrhodata[R+1, C+1])

    ρabove = (
        ρ00a * ((R + 1 - r) * (C + 1 - c)) +
        ρ10a * ((r - R) * (C + 1 - c)) +
        ρ01a * ((R + 1 - r) * (c - C)) +
        ρ11a * ((r - R) * (c - C))
    )

    if pexact == true
        return ρabove
    else
        ρbelow = (
            ρ00b * ((R + 1 - r) * (C + 1 - c)) +
            ρ10b * ((r - R) * (C + 1 - c)) +
            ρ01b * ((R + 1 - r) * (c - C)) +
            ρ11b * ((r - R) * (c - C))
        )
        pslogabove = log(psannual[pindexabove])
        pslogbelow = log(psannual[pindexbelow])
        ρ = (ρabove - ρbelow) / (pslogabove - pslogbelow) * (log(p) - pslogbelow) + ρbelow
        return ρ
    end
end

end # module ItuRP2145
