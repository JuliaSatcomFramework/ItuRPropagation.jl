using TestItemRunner

@testsnippet setup_common begin
    using ItuRPropagation
    using Test
    using XLSX
    validation_file = joinpath(@__DIR__, "CG-3M3J-13-ValEx-Rev8.3.0.xlsx")

    error_tolerance = 1e-7
end
@testitem "Aqua" begin
    using Aqua
    using ItuRPropagation
    Aqua.test_all(ItuRPropagation)
end

@testitem "P.618-14 - Propagation data and prediction methods required for the design of Earth-space telecommunication systems" begin
    include("iturP618test.jl")
end

@testitem "P.835-7 - Standard atmospheric models" begin
    include("iturP835test.jl")
end

@testitem "P.837-7 - Characteristics of precipitation for propagation modelling" begin
    include("iturP837test.jl")
end

@testitem "P.838-3 - Specific attenuation model for rain for use in prediction methods" begin
    include("iturP838test.jl")
end

@run_package_tests verbose=true