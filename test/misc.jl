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
    @test_throws ArgumentError ItuRP618.rainattenuation(LatLon(0, 0), 30, 100, 1)

    @test_logs (:warn, r"between 4 and 55 GHz") match_mode=:any ItuRP618.scintillationattenuation(LatLon(0, 0), 1000, 30, 3)
    @test_logs (:warn, r"between 0.01% and 50%") match_mode=:any ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 1, 0.001)
    @test_logs (:warn, r"between 5 and 90 degrees") match_mode=:any ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 1, 1)
    @test_throws ArgumentError ItuRP618.scintillationattenuation(LatLon(0, 0), 30, 100, 1)



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