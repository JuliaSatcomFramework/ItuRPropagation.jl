@testitem "Errors" begin
    @test_throws ArgumentError ItuRP618.rainattenuation(LatLon(0, 0), 30, 100, 1)
    @test_throws ArgumentError ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 100, 1)
    @test_throws ArgumentError ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 1, 1; efficiency = 101)
    @test_throws "forgot to provide one argument" ItuRP618.attenuations(0, 0, 30, 10; D = 1)

    # LatLon with wrong lat
    @test_throws ArgumentError LatLon(100, 0)
end

@testitem "Warnings" begin
    # Clouds
    @test_logs (:warn, r"between 1 and 200 GHz") match_mode=:any ItuRP840.cloudattenuation(LatLon(0, 0), 1000, 30, 5) 
    @test_logs (:warn, r"between 5 and 90 degrees") match_mode=:any ItuRP840.cloudattenuation(LatLon(0, 0), 100, 1, 1)

    # Interpolations
    @test_logs (:warn, r"liquid water content") match_mode=:any ItuRP840.cloudattenuation(LatLon(0, 0), 100, 30, 99.9) 
    @test_logs (:warn, r"surface temperature") match_mode=:any ItuRP2145.surfacetemperatureannual(LatLon(0, 0), 99.99) 
    @test_logs (:warn, r"surface water vapour density") match_mode=:any ItuRP2145.surfacewatervapourdensityannual(LatLon(0, 0), 99.99) 
    @test_logs (:warn, r"surface total barometric pressure") match_mode=:any ItuRP2145.surfacepressureannual(LatLon(0, 0), 99.99) 
    @test_logs (:warn, r"surface integrated water vapour content") match_mode=:any ItuRP2145.surfacewatervapourcontentannual(LatLon(0, 0), 99.99) 
    @test_logs (:warn, r"wet term of surface refractivity") match_mode=:any ItuRP453.wettermsurfacerefractivityannual(LatLon(0, 0), 99.99) 

    # P618
    @test_logs (:warn, r"between 1 and 55 GHz") match_mode=:any ItuRP618.rainattenuation(LatLon(0, 0), 1000, 30, 3)
    @test_logs (:warn, r"between 0.001% and 5%") match_mode=:any ItuRP618.rainattenuation(LatLon(0, 0), 30, 1, 6)

    @test_logs (:warn, r"between 4 and 55 GHz") match_mode=:any ItuRP618.scintillationattenuation(LatLon(0, 0), 1000, 30, 3)
    @test_logs (:warn, r"between 0.01% and 50%") match_mode=:any ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 1, 0.001)
    @test_logs (:warn, r"between 5 and 90 degrees") match_mode=:any ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 1, 1)
    @test_logs (:warn, r"between 0 and 100") match_mode=:any ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 20, 1; efficiency = .5)

    # Total Attenuations, we just test that scintillation warn is not triggered when p < 0.01 but function is called through attenuations
    @test_logs ItuRP618.attenuations(LatLon(0, 0), 30, 10, 0.001; D = 1)

    # Test suppression of warnings
    ItuRPropagation.SUPPRESS_WARNINGS[] = true
    try
        @test_logs ItuRP840.cloudattenuation(LatLon(0, 0), 1000, 30, 5)
        @test_logs ItuRP840.cloudattenuation(LatLon(0, 0), 100, 1, 1)

        # Interpolations
        @test_logs ItuRP840.cloudattenuation(LatLon(0, 0), 100, 30, 99.9)
        @test_logs ItuRP2145.surfacetemperatureannual(LatLon(0, 0), 99.99)
        @test_logs ItuRP2145.surfacewatervapourdensityannual(LatLon(0, 0), 99.99)
        @test_logs ItuRP2145.surfacepressureannual(LatLon(0, 0), 99.99)
        @test_logs ItuRP2145.surfacewatervapourcontentannual(LatLon(0, 0), 99.99)
        @test_logs ItuRP453.wettermsurfacerefractivityannual(LatLon(0, 0), 99.99)

        # P618
        @test_logs ItuRP618.rainattenuation(LatLon(0, 0), 1000, 30, 3)
        @test_logs ItuRP618.rainattenuation(LatLon(0, 0), 30, 1, 6)
        @test_logs ItuRP618.scintillationattenuation(LatLon(0, 0), 1000, 30, 3)
        @test_logs ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 1, 0.001)
        @test_logs ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 1, 1)
        @test_logs ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 20, 1; efficiency = .5)
    finally
        ItuRPropagation.SUPPRESS_WARNINGS[] = false
    end
end

@testitem "Separate Lat, Lon arguments" begin
    function randll()
        lat, lon = rand() * 180 - 90, rand() * 360 - 180
        lls = (lat, lon)
        ll = LatLon(lat, lon)
        return lls, ll
    end
    lls, ll = randll()
    p = rand() * 3 + 1
    f = 28
    el = 20

    ###### P453 ######
    @test ItuRP453.wettermsurfacerefractivityannual(lls..., p) == ItuRP453.wettermsurfacerefractivityannual(ll, p)
    @test ItuRP453.wettermsurfacerefractivityannual_50(lls...) == ItuRP453.wettermsurfacerefractivityannual_50(ll)

    ###### P618 ######
    @test ItuRP618.rainattenuation(lls..., f, el, p) == ItuRP618.rainattenuation(ll, f, el, p)
    @test ItuRP618.scintillationattenuation(lls..., f, el, p) == ItuRP618.scintillationattenuation(ll, f, el, p)
    @test ItuRP618.attenuations(lls..., f, el, p; D = 1) == ItuRP618.attenuations(ll, f, el, p; D = 1)

    ###### P676 ######
    @test ItuRP676.gaseousattenuation(lls..., f, el, p) == ItuRP676.gaseousattenuation(ll, f, el, p)

    ###### P837 ######
    @test ItuRP837.rainfallrate001(lls...) == ItuRP837.rainfallrate001(ll)

    ###### P839 ######
    @test ItuRP839.isothermheight(lls...) == ItuRP839.isothermheight(ll)
    @test ItuRP839.rainheightannual(lls...) == ItuRP839.rainheightannual(ll)

    ###### P840 ######
    @test ItuRP840.liquidwatercontent(lls..., p) == ItuRP840.liquidwatercontent(ll, p)
    @test ItuRP840.cloudattenuation(lls..., f, el, p) == ItuRP840.cloudattenuation(ll, f, el, p)

    ###### P1511 ######
    @test ItuRP1511.topographicheight(lls...) == ItuRP1511.topographicheight(ll)

    ###### P2145 ######
    ## Mean values
    @test ItuRP2145.surfacetemperatureannual(lls...) == ItuRP2145.surfacetemperatureannual(ll)
    @test ItuRP2145.surfacewatervapourdensityannual(lls...) == ItuRP2145.surfacewatervapourdensityannual(ll)
    @test ItuRP2145.surfacepressureannual(lls...) == ItuRP2145.surfacepressureannual(ll)
    @test ItuRP2145.surfacewatervapourcontentannual(lls...) == ItuRP2145.surfacewatervapourcontentannual(ll)

    ## CCDF Values
    @test ItuRP2145.surfacetemperatureannual(lls..., p) == ItuRP2145.surfacetemperatureannual(ll, p)
    @test ItuRP2145.surfacewatervapourdensityannual(lls..., p) == ItuRP2145.surfacewatervapourdensityannual(ll, p)
    @test ItuRP2145.surfacepressureannual(lls..., p) == ItuRP2145.surfacepressureannual(ll, p)
    @test ItuRP2145.surfacewatervapourcontentannual(lls..., p) == ItuRP2145.surfacewatervapourcontentannual(ll, p)
end

@testitem "Implement Interface" begin
    # We create a custom struct that also stores location (in m) and lat/lon in radians
    struct LLA
        lat::Float64
        lon::Float64
        alt::Float64
    end
    lla = LLA(π/6, π/4, 1200)
    Base.convert(::Type{LatLon}, lla::LLA) = LatLon(lla.lat |> rad2deg, lla.lon |> rad2deg)
    ItuRPropagation.altitude_from_location(lla::LLA) = lla.alt / 1e3

    alt_1511 = ItuRP1511.topographicheight(lla)
    @test ItuRPropagation.altitude_from_location(lla) != alt_1511

    atts = ItuRP618.attenuations(lla, 30, 20, 1; D = 1)
    bench_wrongalt = ItuRP618.attenuations(LatLon(30, 45), 30, 20, 1; D = 1)
    bench_correctalt = ItuRP618.attenuations(lla, 30, 20, 1; D = 1, alt = 1.2)
    @test atts.At ≉ bench_wrongalt.At
    @test atts.At ≈ bench_correctalt.At
end

@testitem "Unitful Extension" begin
    using Unitful

    alt_m = 1200u"m"
    f_hz = 20e9u"Hz"
    hr_m = 100u"m"
    el_d = 20u"°"
    el_rad = uconvert(u"rad", el_d)

    alt = uconvert(u"km", alt_m)
    f = uconvert(u"GHz", f_hz)
    hr = uconvert(u"km", hr_m) 
    el = uconvert(u"°", el_rad)
    p = 1
    D = 1

    att_nounit = ItuRP618.attenuations(LatLon(0, 0), f, el, p; D = D, hᵣ = hr)
    att_unit_deg = ItuRP618.attenuations(LatLon(0, 0), f_hz, el_d, p; D = D, hᵣ = hr_m)
    att_unit_rad = ItuRP618.attenuations(LatLon(0, 0), f_hz, el_rad, p; D = D, hᵣ = hr_m)
    # This is providing lat and lon as seprate numbers, not a LatLon object
    att_unit_split = ItuRP618.attenuations(0, 0u"°", f, el, p; D = D)

    @test all(propertynames(att_nounit)) do fld
        getproperty(att_nounit, fld) ≈ getproperty(att_unit_deg, fld) ≈ getproperty(att_unit_rad, fld)
    end

    # Test LatLon constructor with Unitful numbers as inputs
    @test LatLon(0, 0) == LatLon(0u"°", 0u"°") == LatLon(0u"°", 0u"rad") == LatLon(0u"rad", 0u"°") == LatLon(0u"rad", 0u"rad") == LatLon(0, 0u"°") == LatLon(0u"°", 0)

    # We also test implicit conversion works
end

@testitem "CoordRefSystems Extension" begin
    using CoordRefSystems: CoordRefSystems, LatLonAlt, LatLon as CLL
    using Unitful

    ll = CLL(30, 45)
    lla = LatLonAlt(ll.lat, ll.lon, 1200u"m")

    alt_1511 = ItuRP1511.topographicheight(ll)
    @test ItuRPropagation.altitude_from_location(lla) != alt_1511

    atts_lla = ItuRP618.attenuations(lla, 30, 20, 1; D = 1)
    atts_ll = ItuRP618.attenuations(ll, 30, 20, 1; D = 1)
    bench_1511alt = ItuRP618.attenuations(LatLon(30, 45), 30, 20, 1; D = 1)
    bench_customalt = ItuRP618.attenuations(lla, 30, 20, 1; D = 1, alt = 1.2)
    @test atts_lla.At ≉ bench_1511alt.At
    @test atts_ll.At ≈ bench_1511alt.At
    @test atts_lla.At ≈ bench_customalt.At
end