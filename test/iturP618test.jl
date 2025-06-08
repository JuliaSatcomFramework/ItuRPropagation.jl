@testitem "P.618-14 - Scintillation attenuation" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.618-14 A_Scint")
        map(eachrow(sheet["C21:R68"])) do row
            (;
                inputs=(;
                    ll=LatLon(row[1], row[2]),
                    alt=row[3],
                    f=row[5],
                    el=row[6],
                    D=row[7],
                    η=row[8] * 100,
                    p=row[9],
                ),
                out=(;
                    Nwet=row[10],
                    σref=row[11],
                    L=row[12],
                    x=row[13],
                    g=row[14],
                    σ=row[15],
                    As=row[16],
                ),
            )
        end
    end
    for n in eachindex(entries)
        entry = entries[n]
        (; ll, alt, f, el, D, η, p) = entry.inputs
        out = ItuRP618._scintillationattenuation(ll, f, el, p; D, η)
        @test all(propertynames(entry.out)) do fld
            validation = getproperty(entry.out, fld)
            computed = getproperty(out, fld)
            valid = isapprox(computed, validation; rtol=error_tolerance)
            if !valid
                @warn "Scintillation attenuation entry $n, $fld: $computed != $validation"
            end
            valid
        end
    end
end

@testitem "P.618-14 - Rain attenuation" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.618-14 A_Rain")
        map(eachrow(sheet["C23:U86"])) do row
            (;
                inputs=(;
                    ll=LatLon(row[1], row[2]),
                    alt=row[4],
                    f=row[5],
                    el=row[6],
                    polarization_angle=row[7],
                    p=row[8],
                    R001=row[10],
                ),
                out=(;
                    hᵣ = row[9],
                    Lₛ = row[11],
                    Lg = row[12],
                    γᵣ = row[13],
                    r001 = row[14],
                    v001 = row[15],
                    Lₑ = row[16],
                    A001 = row[17],
                    β = row[18],
                ),
                Ap = row[19],
            )
        end
    end
    for n in eachindex(entries)
        entry = entries[n]
        (; ll, alt, f, el, polarization_angle, p, R001) = entry.inputs
        # We first test that the input R001 is within 1e-3 of the one computed based on the specific r001 map
        @test R001 ≈ ItuRP837.rainfallrate001(ll) rtol=1e-3
        # Then we test all the other intermediate outputs
        out = ItuRP618._rainattenuation(ll, f, el, p; polarization_angle, alt, R001)
        @test all(propertynames(entry.out)) do fld
            validation = getproperty(entry.out, fld)
            computed = getproperty(out, fld)
            valid = isapprox(computed, validation; rtol=error_tolerance)
            if !valid
                @warn "Rain attenuation entry $n, $fld: $computed != $validation"
            end
            valid
        end
        # Finally we test the call with the public function
        @test entry.Ap ≈ ItuRP618.rainattenuation(ll, f, el, p; polarization_angle, alt, R001) rtol=error_tolerance
    end
end