@testsnippet setup_2145 begin
    using ItuRPropagation
    using Test
    using XLSX
    validation_file = joinpath(@__DIR__, "CG-3M3J-13-ValEx-Rev8.3.0.xlsx")
end

@testitem "P.2145-0 - Surface temperature model" setup=[setup_2145] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.2145 CLIMATIC_MAPS_GAS_ATT")
        map(eachrow(sheet["C23:N101"])) do row
            (;
                ll = LatLon(row[1], row[2]),
                alt = row[3],
                p = row[4],
                RHO_mean = row[5],
                P_mean = row[7],
                T_mean = row[7],
                V_mean = row[8],
                RHO_ccdf = row[9],
                P_ccdf = row[10],
                T_ccdf = row[11],
                V_ccdf = row[12],
            )
        end
    end

    error_tolerance = 1e-7

    @testset "Average Temperature T̄ₛ" begin
        for entry in entries
            (; ll, alt, T_mean) = entry
            T̄ₛ = ItuRP2145.surfacetemperatureannual2(ll; alt)
            @test T̄ₛ ≈ T_mean rtol=error_tolerance
        end
    end

    @testset "Average Surface Water Vapour Density ρ̄ₛ" begin
        for entry in entries
            (; ll, alt, RHO_mean) = entry
            ρ̄ₛ = ItuRP2145.surfacewatervapordensityannual2(ll; alt)
            @test ρ̄ₛ ≈ RHO_mean rtol=error_tolerance
        end
    end

    @testset "CCDF Temperature Tₛ(p)" begin
        for entry in entries
            (; ll, alt, p, T_ccdf) = entry
            Tₛ = ItuRP2145.surfacetemperatureannual2(ll, p; alt)
            @test Tₛ ≈ T_ccdf rtol=error_tolerance
        end
    end
end
