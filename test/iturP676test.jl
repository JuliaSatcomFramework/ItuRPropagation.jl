@testitem "P.676-13 - Specific Attenuation Coefficients" setup = [setup_common] begin

    @testset "Oxygen coefficients" begin
        entries = XLSX.openxlsx(validation_file) do wb
            sheet = XLSX.getsheet(wb, "P.676-13 SpAtt")
            fs = map(sheet["E19:I19"] |> vec) do fstr
                parse(Float64, fstr[1:end-4])
            end
            map(zip(fs, eachcol(sheet["E20:I28"]))) do (f, col)
                (;
                    P=col[1],
                    e=col[2],
                    T=col[3],
                    ρ=col[4],
                    S₁=col[5],
                    F₁=col[6],
                    d=col[7],
                    Nppₒ=col[8],
                    γₒ=col[9],
                    f,
                )
            end
        end
        tbl1_row = ItuRP676.TABLE1_DATA[1]
        for entry in entries
            (; P, e, T, ρ, f) = entry
            θ = 300 / T
            @test entry.S₁ ≈ ItuRP676._Sₒ(tbl1_row, θ, P) rtol = error_tolerance
            @test entry.F₁ ≈ ItuRP676._Fₒ(tbl1_row, f, θ, P, e) rtol = error_tolerance
            calc = ItuRP676._gammaoxygen(f, T, P, ρ)
            for fld in (:γₒ, :Nppₒ, :d)
                validation = getproperty(entry, fld)
                computed = getproperty(calc, fld)
                @test validation ≈ computed rtol = error_tolerance
            end
        end
    end

    @testset "Oxygen coefficients" begin
        entries = XLSX.openxlsx(validation_file) do wb
            sheet = XLSX.getsheet(wb, "P.676-13 SpAtt")
            fs = map(sheet["N19:R19"] |> vec) do fstr
                parse(Float64, fstr[1:end-4])
            end
            map(zip(fs, eachcol(sheet["N20:R28"]))) do (f, col)
                (;
                    P=col[1],
                    e=col[2],
                    T=col[3],
                    ρ=col[4],
                    S₁=col[5],
                    F₁=col[6],
                    # d = col[7],
                    Nppᵥ=col[8],
                    γᵥ=col[9],
                    f,
                )
            end
        end
        tbl2_row = ItuRP676.TABLE2_DATA[1]
        for entry in entries
            (; P, e, T, ρ, f) = entry
            θ = 300 / T
            @test entry.S₁ ≈ ItuRP676._Sᵥ(tbl2_row, θ, e) rtol = error_tolerance
            @test entry.F₁ ≈ ItuRP676._Fᵥ(tbl2_row, f, θ, P, e) rtol = error_tolerance
            calc = ItuRP676._gammawater(f, T, P, ρ)
            for fld in (:γᵥ, :Nppᵥ)
                validation = getproperty(entry, fld)
                computed = getproperty(calc, fld)
                @test validation ≈ computed rtol = error_tolerance
            end
        end
    end


    @testset "Oxygen + Water Vapour coefficients" begin
        entries = XLSX.openxlsx(validation_file) do wb
            sheet = XLSX.getsheet(wb, "P.676-13 SpAtt")
            map(eachrow(sheet["C34:F383"])) do row
                (;
                    f=row[1],
                    γₒ=row[2],
                    γᵥ=row[3],
                    γ=row[4]
                )
            end
        end
        # These 3 values of pressure, temperature and water vapour density are not specifically given in the subtable, so taking them from the previous two tables
        P = 1013.25
        T = 288.15
        ρ = 7.5
        for entry in entries
            calc = ItuRP676._gamma(entry.f, T, P, ρ)
            for fld in (:γ, :γₒ, :γᵥ)
                validation = getproperty(entry, fld)
                computed = getproperty(calc, fld)
                @test validation ≈ computed rtol = error_tolerance
            end
        end
    end
end

@testitem "P.676-13 - Slant Path Layered Attenuation" setup = [setup_common] begin
    expected_attenuation = 0.0
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.676-13 A_Gas_A1_2.2.1a")
        global expected_attenuation = sheet["AB23"]
        map(eachrow(sheet["H23:Y944"])) do row
            layer = (;
                δ = row[1],
                r = row[2],
                h = row[4],
                P = row[6],
                T = row[7],
                ρ = row[8],
                Pd = row[9],
                e = row[10],
                n = row[11],
            )
            outs = (;
                β = row[12],
                a = row[14],
                γₒ = row[16],
                γᵥ = row[17],
                γ = row[18],
            )
            (; layer, outs)
        end
    end
    layers = ItuRP676.STANDARD_LAYERS
    f = 28
    el = 30
    @testset "Per Layer results" begin
        for (n, entry) in enumerate(entries)
            # Test the basic layer-specific data
            @test all(propertynames(entry.layer)) do fld
                implemented = getproperty(layers[n], fld)
                validation = getproperty(entry.layer, fld)
                isapprox(implemented, validation; rtol = error_tolerance)
            end
            # Test the actual values of attenuation computation per layer
            outs = ItuRP676.layerattenuation(layers[n], f, el)
            @test all(propertynames(entry.outs)) do fld
                validation = getproperty(entry.outs, fld)
                implemented = getproperty(outs, fld)
                isapprox(implemented, validation; rtol = error_tolerance)
            end
        end
    end
    @testset "Total Attenuation" begin
        @test expected_attenuation ≈ ItuRP676._gasattenuation_layers(layers, f, el) rtol = error_tolerance
    end
end

@testitem "P.676-13 - Slant Path Statistical (Annual) Attenuation" setup = [setup_common] begin
    entries = XLSX.openxlsx(validation_file) do wb
        sheet = XLSX.getsheet(wb, "P.676-13 A_Gas_A2_STAT")
        map(eachrow(sheet["C24:AP100"])) do row
            aux = (;
                latlon = LatLon(row[1], row[2]),
                alt = row[3],
                f = row[5], # GHz
                P̄ = row[6], # hPa
                el = row[7], # °
                ρ̄ = row[8], # g/m^3
                T̄ = row[9], # K
                ē = row[10], # hPa
                P̄d = row[11], # hPa
                p = row[12], # %
                P = row[13], # hPa
                T = row[14], # K
                ρ = row[15], # g/m^3
            )
            oxygen = (;
                γₒ = row[16], # dB/km
                aₒ = row[17],
                bₒ = row[18],
                cₒ = row[19],
                dₒ = row[20],
                hₒ = row[21],
                Aₒ_zenith = row[22], # dB
                Aₒ = row[23], # dB
            )
            (; aux, oxygen)
        end
    end

    @testset "Auxiliary Values" begin
        for entry in entries
            (; latlon, f, alt, P̄, T̄, ρ̄, ē, P̄d, p, P, T, ρ) = entry.aux
            @test alt ≈ ItuRP1511.topographicheight(latlon) rtol = error_tolerance
            # We check the 2145 values
            @test P̄ ≈ ItuRP2145.surfacepressureannual(latlon; alt) rtol = error_tolerance
            @test T̄ ≈ ItuRP2145.surfacetemperatureannual(latlon; alt) rtol = error_tolerance
            @test ρ̄ ≈ ItuRP2145.surfacewatervapourdensityannual(latlon; alt) rtol = error_tolerance
            @test P ≈ ItuRP2145.surfacepressureannual(latlon, p; alt) rtol = error_tolerance
            @test T ≈ ItuRP2145.surfacetemperatureannual(latlon, p; alt) rtol = error_tolerance
            @test ρ ≈ ItuRP2145.surfacewatervapourdensityannual(latlon, p; alt) rtol = error_tolerance
        end
    end

    @testset "Oxygen" begin
        for (n, entry) in enumerate(entries)
            (; latlon, f, alt, el, p) = entry.aux
            out = ItuRP676._slantoxygenattenuation(latlon, f, el, p; alt)
            @test all(propertynames(entry.oxygen)) do fld
                validation = getproperty(entry.oxygen, fld)
                computed = getproperty(out, fld)
                valid = isapprox(computed, validation; rtol = error_tolerance)
                if !valid
                    @warn "Oxygen entry $n, $fld: $computed != $validation"
                end
                valid
            end
        end
    end

end
