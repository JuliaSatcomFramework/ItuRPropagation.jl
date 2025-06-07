@testitem "P.839-4 - Rain Height Model" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.839-4 Rain_Height")
        map(eachrow(sheet["C22:G29"])) do row
            (;
                ll=LatLon(row[1], row[2]),
                h0=row[4],
                hR=row[5],
            )
        end
    end
    for entry in entries
        (; ll, h0, hR) = entry
        @test ItuRP839.isothermheight(ll) ≈ h0 rtol = error_tolerance
        @test ItuRP839.rainheightannual(ll) ≈ hR rtol = error_tolerance
    end
end