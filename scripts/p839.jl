include(joinpath(@__DIR__, "common.jl"))

function create_p839_artifact()
    artifact_name = "p839"
    _artifact_hash = Artifacts.artifact_hash(artifact_name, artifact_toml)
    
    if isnothing(_artifact_hash) || !Artifacts.artifact_exists(_artifact_hash)
        _artifact_hash = Artifacts.create_artifact() do artifact_folder
            # This is the URL for direct download of the zip file containing the components of the ITU-R P.839-4 recommendation
            url = "https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.839-4-201309-I!!ZIP-E.zip"

            zip_path = joinpath(downloads_dir, "p839.zip")
            if !isfile(zip_path)
                @info "Downloading raw zip file for ITU-R P.839-4 from ITU Database"
                Downloads.download(url, zip_path)
                @info "Download completed"
            end

            latres = 1.5
            lonres = 1.5

            latrange = range(-90, 90, step=latres)
            lonrange = range(-180, 180, step=lonres)

            matsize = (length(latrange), length(lonrange))

            # Add a README
            open(joinpath(artifact_folder, "README"), "w") do io
                println(io, "This folder contains Maps for the mean annual isotherm height above mean sea level, for use in the ITU-R P.839-4 recommendation, directly in binary format.")
                println(io)
                println(io, "It was automatically generated from the TXT files of the ITU-R P.839-4 recommendation available at the following URL:")
                println(io)
                println(io, url)
                println(io, "The matrices stored into each of the binary files corresponds to a lat/lon grid with $(latres)°/$(lonres)° corresponding resolution and have a size of $(matsize) elements.")
                println(io, "The top-left corner has negative latitude and longitude.")
                println(io)
                println(io, "This artifact was automatically generated using the script at the following URL:")
                println(io, "https://github.com/JuliaSatcomFramework/ItuRPropagation.jl/blob/ITU_artifacts/scripts/p837_r001.jl")
            end

            archive = ZipReader(read(zip_path))

            name = "h0.txt"
            @info "Converting file $name"
            filecontent = zip_readentry(archive, name)
            binfile = joinpath(artifact_folder, replace(name, ".txt" => ".bin"))
            data = readdlm(filecontent, '\t')
            # We have to permute over lat as we want negative latitude at the top
            data = reverse(data; dims=1)
            # Now we want to circshift longitude as we want it from -180 to 180
            data = circshift(data, (0, -(last(matsize) ÷ 2)))
            size(data) == matsize || error("Unexpected size: $(size(data)) instead of $matsize for file $name")
            open(binfile, "w") do io
                write(io, data)
            end
        end
    end


    asset_name = "p839.tar.gz"
    tarball_path = joinpath(assets_dir, asset_name)
    tarball_sha = if !isfile(tarball_path)
        @info "Creating the artifact tarball"
        Artifacts.archive_artifact(_artifact_hash, tarball_path)
    else
        sha256sum(tarball_path)
    end
    release_url = release_root_url * asset_name
    @info "Updating the Artifacts.toml file"
    Artifacts.bind_artifact!(artifact_toml, artifact_name, _artifact_hash; force=true, download_info=[(release_url, tarball_sha)])
end

create_p839_artifact()