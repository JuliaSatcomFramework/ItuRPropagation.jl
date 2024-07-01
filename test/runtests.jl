using TestItemRunner

@testitem "P.618-14 - Propagation data and prediction methods required for the design of Earth-space telecommunication systems" begin
    include("iturP618test.jl")
end

@testitem "P.676-13 - Attenuation by atmospheric gases and related effects" begin
    include("iturP676test.jl")
end

@testitem "P.837-7 - Characteristics of precipitation for propagation modelling" begin
    include("iturP837test.jl")
end

@testitem "P.838-3 - Specific attenuation model for rain for use in prediction methods" begin
    include("iturP838test.jl")
end

@testitem "P.839-4 - Rain Height Model" begin
    include("iturP839test.jl")
end

@testitem "P.840-9 - Attenuation due to clouds and fog" begin
    include("iturP840test.jl")
end

@testitem "P.1511-2 - Topography for Earth-space propagation modelling
" begin
    include("iturP1511test.jl")
end

@run_package_tests verbose=true