@testitem "P.676-13 - Specific Attenuation Coefficients" setup=[setup_common] begin

    @testset "Oxygen coefficients" begin
        entries = XLSX.openxlsx(validation_file) do wb
            sheet = XLSX.getsheet(wb, "P.676-13 SpAtt")
            fs = map(sheet["E19:I19"] |> vec) do fstr
                parse(Float64, fstr[1:end-4])
            end
            map(zip(fs, eachcol(sheet["E20:I28"]))) do (f, col)
                (;
                    P = col[1],
                    e = col[2],
                    T = col[3],
                    ρ = col[4],
                    S₁ = col[5],
                    F₁ = col[6],
                    d = col[7],
                    Nppₒ = col[8],
                    γₒ = col[9],
                    f,
                )
            end
        end
        tbl1_row = ItuRP676.TABLE1_DATA[1]
        for entry in entries
            (; P, e, T, ρ, f) = entry
            θ = 300 / T
            @test entry.S₁ ≈ ItuRP676._Sₒ(tbl1_row, θ, P) rtol=error_tolerance
            @test entry.F₁ ≈ ItuRP676._Fₒ(tbl1_row, f, θ, P, e) rtol=error_tolerance
            calc = ItuRP676._gammaoxygen(f, T, P, ρ)
            for fld in (:γₒ, :Nppₒ, :d)
                validation = getproperty(entry, fld)
                computed = getproperty(calc, fld)
                @test validation ≈ computed rtol=error_tolerance
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
                    P = col[1],
                    e = col[2],
                    T = col[3],
                    ρ = col[4],
                    S₁ = col[5],
                    F₁ = col[6],
                    # d = col[7],
                    Nppᵥ = col[8],
                    γᵥ = col[9],
                    f,
                )
            end
        end
        tbl2_row = ItuRP676.TABLE2_DATA[1]
        for entry in entries
            (; P, e, T, ρ, f) = entry
            θ = 300 / T
            @test entry.S₁ ≈ ItuRP676._Sᵥ(tbl2_row, θ, e) rtol=error_tolerance
            @test entry.F₁ ≈ ItuRP676._Fᵥ(tbl2_row, f, θ, P, e) rtol=error_tolerance
            calc = ItuRP676._gammawater(f, T, P, ρ)
            for fld in (:γᵥ, :Nppᵥ)
                validation = getproperty(entry, fld)
                computed = getproperty(calc, fld)
                @test validation ≈ computed rtol=error_tolerance
            end
        end
    end


    @testset "Oxygen + Water Vapour coefficients" begin
        entries = XLSX.openxlsx(validation_file) do wb
            sheet = XLSX.getsheet(wb, "P.676-13 SpAtt")
            map(eachrow(sheet["C34:F383"])) do row
                (;
                    f = row[1],
                    γₒ = row[2],
                    γᵥ = row[3],
                    γ = row[4]
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
                @test validation ≈ computed rtol=error_tolerance
            end
        end
    end
end
