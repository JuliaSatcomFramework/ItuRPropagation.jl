using Pkg: Artifacts
using Downloads
using SHA
using ZipArchives
using DelimitedFiles

function sha256sum(tarball_path)
    return open(tarball_path, "r") do io
        return bytes2hex(sha256(io))
    end
end
downloads_dir = joinpath(@__DIR__, "..", "downloads")
isdir(downloads_dir) || mkdir(downloads_dir)

assets_dir = joinpath(@__DIR__, "..", "assets")
isdir(assets_dir) || mkdir(assets_dir)

artifact_toml = joinpath(@__DIR__, "..", "Artifacts.toml")

release_root_url = "https://github.com/JuliaSatcomFramework/ItuRPropagation.jl/releases/download/artifact_releases/"

function parseline(str::AbstractString, ::Type{T} = Float64) where T
    cleaned = replace(strip(str), r" +" => ' ')
    return map(x -> parse(T, x), split(cleaned, ' '))
end

function parsematrix(file::AbstractString, ::Type{T} = Float64) where T
    return stack(s -> parseline(s, T), eachline(file)) |> permutedims
end

