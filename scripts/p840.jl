include(joinpath(@__DIR__, "common.jl"))

function create_p840_annual_artifact()
    artifact_name = "p840_annual"
    _artifact_hash = Artifacts.artifact_hash(artifact_name, artifact_toml)

    if isnothing(_artifact_hash) || !Artifacts.artifact_exists(_artifact_hash)
        _artifact_hash = Artifacts.create_artifact() do artifact_folder
            # This is the URL for direct download of the zip file containing the Part 1 of the ITU-R P.840-8 recommendation (Annual Data)
            url = "https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.840Part01-0-202308-I!!ZIP-E.zip"

            zip_path = joinpath(downloads_dir, "p840_2023_Part1_annual.zip")
            if !isfile(zip_path)
                @info "Downloading raw zip file for ITU-R P.840-8 Part 1 (Annual Data) from ITU Database"
                Downloads.download(url, zip_path)
                @info "Download completed"
            end


            latres = 0.25
            lonres = 0.25

            latrange = range(-90, 90, step=latres)
            lonrange = range(0, 360, step=lonres)

            matsize = (length(latrange), length(lonrange))

            # Add a README
            open(joinpath(artifact_folder, "README"), "w") do io
                println(io, "This folder contains data for the ITU-R P.840-8 recommendation (Annual Data of Integrated Liquid Water Content) directly in binary format.")
                println(io)
                println(io, "It was automatically generated from the TXT files of the Part 1 annex of the ITU-R P.840-8 recommendation available at the following URL:")
                println(io)
                println(io, url)
                println(io, "The matrices stored into each of the binary files corresponds to a lat/lon grid with $(latres)° x $(lonres)° resolution and have a size of $(matsize) elements.")
                println(io)
                println(io, "This artifact was automatically generated using the script at the following URL:")
                println(io, "https://github.com/JuliaSatcomFramework/ItuRPropagation.jl/blob/ITU_artifacts/scripts/p840.jl")
            end

            archive = ZipReader(read(zip_path))
            for name in zip_names(archive)
                endswith(lowercase(name), "txt") || throw(ArgumentError("Unexpected file type: $name"))
                @info "converting file $name"
                filecontent = zip_readentry(archive, name)
                binfile = joinpath(artifact_folder, replace(name, ".TXT" => ".bin"))
                data = readdlm(filecontent, ' ')
                size(data) == matsize || error("Unexpected size: $(size(data)) instead of $matsize for file $name")
                open(binfile, "w") do io
                    write(io, data)
                end
            end

            # We also add the original ITU Readme which is stored in a separate url (Part 15)
            readme_url = "https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.840Part15-0-202308-I!!MSW-E.docx"
            itu_readme_path = joinpath(artifact_folder, "ITU_README.docx")
            Downloads.download(readme_url, itu_readme_path)
        end
    end


    asset_name = "p840_2023_Part1_annual.tar.gz"
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

create_p840_annual_artifact()
