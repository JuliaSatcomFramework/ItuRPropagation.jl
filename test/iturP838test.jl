@testitem "P.838-3 - Specific attenuation model for rain for use in prediction methods" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.838-3 Sp.Att")
        map(eachrow(sheet["C20:P83"])) do row
            (;
                inputs=(;
                    ll=LatLon(row[1], row[2]),
                    el = row[4],
                    f=row[5],
                    R=row[6],
                    polarization_angle=row[7],
                ),
                out=(;
                    kₕ=row[8],
                    kᵥ=row[9],
                    αₕ=row[10],
                    αᵥ=row[11],
                    k=row[12],
                    α=row[13],
                    γᵣ=row[14],
                ),
            )
        end
    end
    for n in eachindex(entries)
        entry = entries[n]
        (; ll, f, el, R, polarization_angle) = entry.inputs
        out = ItuRP838._rainspecificattenuation(f, el; R, polarization_angle)
        @test all(propertynames(entry.out)) do fld
            validation = getproperty(entry.out, fld)
            computed = getproperty(out, fld)
            valid = isapprox(computed, validation; rtol=error_tolerance)
            if !valid
                @warn "Specific rain attenuation entry $n, $fld: $computed != $validation"
            end
            valid
        end
    end
end