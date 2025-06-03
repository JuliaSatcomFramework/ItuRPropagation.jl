@testitem "P.1511-3 - Topographic height model" setup=[setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.1511 TOPO")
        map(eachrow(sheet["C23:E37"])) do row
            (;
                ll = LatLon(row[1], row[2]),
                alt = row[3] / 1e3, # input is in m but we output km
            )
        end
    end

    error_tolerance = 1e-7

    @testset "Topographic height" begin
        for entry in entries
            (; ll, alt) = entry
            @test ItuRP1511.topographicheight(ll) ≈ alt rtol=error_tolerance
        end
    end
end