@testitem "P.837-7 - Characteristics of precipitation for propagation modelling" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.837-7 Rp")
        map(eachrow(sheet["C69:F76"])) do row
            (;
                ll=LatLon(row[1], row[2]),
                Rp=row[4],
            )
        end
    end
    for entry in entries
        (; ll, Rp) = entry
        @test ItuRP837.rainfallrate001(ll) ≈ Rp rtol = error_tolerance
    end
end
