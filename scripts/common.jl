using Pkg.Artifacts
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

artifact_toml = joinpath(@__DIR__, "..", "Artifacts.toml")
