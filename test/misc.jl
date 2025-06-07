@testitem "Warnings" begin
    # Clouds
    @test_logs (:warn, r"between 1 and 200 GHz") ItuRP840.cloudattenuation(LatLon(0, 0), 1000, 30, 5) match_mode=:any
    @test_logs (:warn, r"between 5 and 90 degrees") ItuRP840.cloudattenuation(LatLon(0, 0), 100, 1, 1) match_mode=:any

    # Interpolations
    @test_logs (:warn, r"liquid water content") ItuRP840.cloudattenuation(LatLon(0, 0), 100, 30, 99.9) match_mode=:any
    @test_logs (:warn, r"surface temperature") ItuRP2145.surfacetemperatureannual(LatLon(0, 0), 99.99) match_mode=:any
    @test_logs (:warn, r"surface water vapour density") ItuRP2145.surfacewatervapourdensityannual(LatLon(0, 0), 99.99) match_mode=:any
    @test_logs (:warn, r"surface total barometric pressure") ItuRP2145.surfacepressureannual(LatLon(0, 0), 99.99) match_mode=:any
    @test_logs (:warn, r"surface integrated water vapour content") ItuRP2145.surfacewatervapourcontentannual(LatLon(0, 0), 99.99) match_mode=:any
    @test_logs (:warn, r"wet term of surface refractivity") ItuRP453.wettermsurfacerefractivityannual(LatLon(0, 0), 99.99) match_mode=:any

    # Test suppression of warnings
    ItuRPropagation.SUPPRESS_WARNINGS[] = true
    try
        @test_nowarn ItuRP840.cloudattenuation(LatLon(0, 0), 1000, 30, 5)
        @test_nowarn ItuRP840.cloudattenuation(LatLon(0, 0), 100, 1, 1)

        # Interpolations
        @test_nowarn ItuRP840.cloudattenuation(LatLon(0, 0), 100, 30, 99.9)
        @test_nowarn ItuRP2145.surfacetemperatureannual(LatLon(0, 0), 99.99)
        @test_nowarn ItuRP2145.surfacewatervapourdensityannual(LatLon(0, 0), 99.99)
        @test_nowarn ItuRP2145.surfacepressureannual(LatLon(0, 0), 99.99)
        @test_nowarn ItuRP2145.surfacewatervapourcontentannual(LatLon(0, 0), 99.99)
        @test_nowarn ItuRP453.wettermsurfacerefractivityannual(LatLon(0, 0), 99.99)
    finally
        ItuRPropagation.SUPPRESS_WARNINGS[] = false
    end
end