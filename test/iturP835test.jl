@testitem "P.835-7 - Standard atmospheric models" begin
    # Here we only do coverage testing as there are no official validation tests for 835
    ρ = map(range(0, 100, 100)) do Z
        ItuRP835.standardwatervapordensity(Z)
    end
    @test ρ[1] == 7.5
end





