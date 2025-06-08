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

@testitem "P.618-14 - Total Attenuations" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.618-14 Att_Tot")
        map(eachrow(sheet["C26:S73"])) do row
            (;
                inputs=(;
                    ll=LatLon(row[1], row[2]),
                    alt=row[3],
                    f=row[5],
                    el=row[6],
                    D = row[7],
                    η = row[8] * 100,
                    polarization_angle=row[9],
                    p=row[10],
                ),
                gas_intermediate=(;
                    Aox = row[11],
                    Awv = row[12],
                ),
                out=(;
                    Ag = row[13],
                    Ac = row[14],
                    Ar = row[15],
                    As = row[16],
                    At = row[17],
                ),
            )
        end
    end
    for n in eachindex(entries)
        entry = entries[n]
        (; ll, alt, f, el, polarization_angle, D, η, p) = entry.inputs
        out = ItuRP618.attenuations(ll, f, el, p; polarization_angle, alt, D, η)
        @test all(propertynames(entry.out)) do fld
            validation = getproperty(entry.out, fld)
            computed = getproperty(out, fld)
            rtol = fld in (:At, :Ar) ? 5e-4 : error_tolerance # We have to use higher tolerance for rain (and total) as that uses the more precise method for R001
            valid = isapprox(computed, validation; rtol)
            if !valid
                @warn "Total attenuation entry $n, $fld: $computed != $validation"
            end
            valid
        end
    end
end

@testitem "P.618-14 - Total Attenuation CCDF Samples" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.618-14 Att_Tot")
        inputs = let 
            row = sheet["C81:K81"]
            (;
                ll=LatLon(row[1], row[2]),
                alt=row[3],
                f=row[5],
                el=row[6],
                D = row[7],
                η = row[8] * 100,
                polarization_angle=row[9],
            )
        end
        map(eachrow(sheet["C88:K107"])) do row
            (;
                inputs = (;
                    inputs...,
                    p = row[1]
                ),
                out=(;
                    Ag = row[7],
                    Ac = row[6],
                    Ar = row[5],
                    As = row[4],
                    At = row[9],
                ),
            )
        end
    end
    for n in eachindex(entries)
        entry = entries[n]
        (; ll, alt, f, el, polarization_angle, D, η, p) = entry.inputs
        out = ItuRP618.attenuations(ll, f, el, p; polarization_angle, alt, D, η)
        @test all(propertynames(entry.out)) do fld
            validation = getproperty(entry.out, fld)
            computed = getproperty(out, fld)
            rtol = fld in (:At, :Ar) ? 5e-4 : error_tolerance # We have to use higher tolerance for rain (and total) as that uses the more precise method for R001
            valid = isapprox(computed, validation; rtol)
            if !valid
                @warn "Total attenuation entry $n, $fld: $computed != $validation"
            end
            valid
        end
    end
end

@testitem "P.618-14 - Total Attenuations fast computation" setup = [setup_common] begin
    D = 1
    for n in 1:20
        f = rand(18:50)
        p = rand() * 49.999 + .001
        latlon = LatLon(rand() * 180 - 90, rand() * 360 - 180)
        # Compute the stable kwargs for faster computation. Elevation is irrelevant here as we only need the zenith values which are contained in the `kwargs` field
        (; kwargs) = ItuRP618._attenuations(latlon, f, rand() * 50 + 10, p; D)
        for _ in 1:20
            el = rand() * 85 + 5 # Elevation is irrelevant here as we only need the zenith values
            out_standard = ItuRP618.attenuations(latlon, f, el, p; D)
            out_fast = ItuRP618._attenuations(latlon, f, el, p; D, kwargs...).attenuations
            @test all(propertynames(out_standard)) do fld
                isapprox(getproperty(out_standard, fld), getproperty(out_fast, fld))
            end
        end
    end

end
