@testitem "P.453-14 - Wet term of surface refractivity Interpolation" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.453-14 Nwet")
        map(eachrow(sheet["C21:F28"])) do row
            (;
                ll=LatLon(row[1], row[2]),
                Nwet = row[4],
            )
        end
    end
    for entry in entries
        (; ll, Nwet) = entry
        direct = ItuRP453.wettermsurfacerefractivityannual_50(ll)
        itp = ItuRP453.wettermsurfacerefractivityannual(ll, 50)
        @test itp ≈ direct
        @test Nwet ≈ direct rtol = error_tolerance 
    end
end

