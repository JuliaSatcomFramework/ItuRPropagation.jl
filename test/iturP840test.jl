@testitem "P.840-9 - Attenuation due to clouds and fog" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.840-9 A_Clouds")
        map(eachrow(sheet["C21:M52"])) do row
            ll = LatLon(row[1], row[2])
            p = row[3]
            f = row[4]
            el = row[5]
            K = (;
                ϵ′=row[6],
                ϵ′′=row[7],
                η=row[8],
                K_L=row[9],
            )
            L = row[10]
            Ac = row[11]
            (; ll, p, f, el, K, L, Ac)
        end
    end
    for entry in entries
        (; ll, p, f, el, L, Ac) = entry
        out = ItuRP840._K_L(f)
        @test all(propertynames(entry.K)) do fld
            validation = getproperty(entry.K, fld)
            computed = getproperty(out, fld)
            valid = isapprox(computed, validation; rtol=error_tolerance)
            if !valid
                @warn "Kₗ entry $n, $fld: $computed != $validation"
            end
            valid
        end
        L = ItuRP840.liquidwatercontent(ll, p)
        Ac = ItuRP840.cloudattenuation(ll, f, el, p)
        @test L ≈ L rtol = error_tolerance
        @test Ac ≈ Ac rtol = error_tolerance
    end
end

@testitem "P.840-9 - Liquid Water Content Interpolation" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P840_9_L")
        map(eachrow(sheet["C24:G59"])) do row
            (;
                ll=LatLon(row[1], row[2]),
                p=row[3],
                L=row[5],
            )
        end
    end
    for entry in entries
        (; ll, p, L) = entry
        @test ItuRP840.liquidwatercontent(ll, p) ≈ L rtol = error_tolerance
    end
end