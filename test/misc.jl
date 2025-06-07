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
    finally
        ItuRPropagation.SUPPRESS_WARNINGS[] = false
    end
end