include(joinpath(@__DIR__, "common.jl"))

function create_p453_nwet_annual_artifact()
    artifact_name = "p453_nwet_annual"
    _artifact_hash = Artifacts.artifact_hash(artifact_name, artifact_toml)

    if isnothing(_artifact_hash) || !Artifacts.artifact_exists(_artifact_hash)
        _artifact_hash = Artifacts.create_artifact() do artifact_folder
            # This is the URL for direct download of the zip file containing the components of the ITU-R P.453-14 recommendation
            url = "https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.453-14-201908-I!!ZIP-E.zip"

            zip_path = joinpath(downloads_dir, "p453_components.zip")
            if !isfile(zip_path)
                @info "Downloading raw zip file for ITU-R P.453-14 components from ITU Database"
                Downloads.download(url, zip_path)
                @info "Download completed"
            end

            latres = 0.75
            lonres = 0.75

            latrange = range(-90, 90, step=latres)
            lonrange = range(0, 360, step=lonres)

            matsize = (length(latrange), length(lonrange))

            # Add a README
            open(joinpath(artifact_folder, "README"), "w") do io
                println(io, "This folder contains annual Maps for the wet term of surface refractivity to be used for the ITU-R P.453-14 recommendation, directly in binary format.")
                println(io)
                println(io, "It was automatically generated from the TXT files of components file of the ITU-R P.453-14 recommendation available at the following URL:")
                println(io)
                println(io, url)
                println(io, "The matrices stored into each of the binary files corresponds to a lat/lon grid with $(latres)°/$(lonres)° corresponding resolution and have a size of $(matsize) elements.")
                println(io)
                println(io, "This artifact was automatically generated using the script at the following URL:")
                println(io, "https://github.com/JuliaSatcomFramework/ItuRPropagation.jl/blob/ITU_artifacts/scripts/p453_nwet_annual.jl")
            end

            top = ZipReader(read(zip_path))
            archive = let 
                mid = ZipReader(zip_readentry(top, "P.453_NWET_Maps.zip"))
                bottom = ZipReader(zip_readentry(mid, "P.453_NWET_Maps_Annual.zip"))
            end
            # We copy the original docx readme for reference
            itu_readme = "Readme_P.453.docx"
            open(joinpath(artifact_folder, "ITU_" * itu_readme), "w") do io
                write(io, zip_readentry(top, itu_readme))
            end


            for name in zip_names(archive)
                startswith(name, "NWET") || continue
                @info "Converting file $name"
                endswith(name, "TXT") || error("Unexpected file type: $name")
                filecontent = zip_readentry(archive, name)
                binfile = joinpath(artifact_folder, replace(name, ".TXT" => ".bin"))
                data = readdlm(filecontent, ' ')
                size(data) == matsize || error("Unexpected size: $(size(data)) instead of $matsize for file $name")
                open(binfile, "w") do io
                    write(io, data)
                end
            end
        end
    end


    asset_name = "p453_nwet_annual.tar.gz"
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

create_p453_nwet_annual_artifact()
