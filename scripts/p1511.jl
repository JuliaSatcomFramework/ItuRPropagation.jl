include(joinpath(@__DIR__, "common.jl"))

function create_p1511_artifact()
    artifact_name = "p1511"
    _artifact_hash = Artifacts.artifact_hash(artifact_name, artifact_toml)

    if isnothing(_artifact_hash) || !Artifacts.artifact_exists(_artifact_hash)
        _artifact_hash = Artifacts.create_artifact() do artifact_folder
            # This is the URL for direct download of the zip file containing the files applicable for ITU-R P.1511-3 recommendation
            url = "https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.1511-3-202408-I!!ZIP-E.zip"

            zip_path = joinpath(downloads_dir, "p1511_data.zip")
            if !isfile(zip_path)
                @info "Downloading raw zip file for ITU-R P.1511-3 recommendation from ITU Database"
                Downloads.download(url, zip_path)
                @info "Download completed"
            end


            # Add a README
            open(joinpath(artifact_folder, "README"), "w") do io
                println(io, "This folder contains data for the ITU-R P.1511-3 recommendation.")
                println(io)
                println(io, "It was automatically generated from the TXT/dat files of the ITU-R P.1511-3 recommendation available at the following URL:")
                println(io)
                println(io, url)
                println(io)
                println(io, "It also still contains the readmes originally present in the ITU data")
                println(io)
                println(io, "The matrices stored in this artifact contains values in meters and have slightly different grids but with a fixed resolution of 1/12 degress in both lat and lon.")
                println(io, "For both files/variables, the grid _wraps_ around the edges of the grid to simplify the bi-cubic interpolation.")
                println(io, "- For the TOPO.dat file, the grid goes from -90.125 to 90.125 degrees in latitude and from -180.125 to 180.125 degrees in longitude. This results in a matrix of 2164 x 4324 elements.")
                println(io, "- For the EGM2008.txt file, the grid goes from -(90 + 2/12) to (90 + 2/12) degrees in latitude and from -(180 + 2/12) to (180 + 2/12)  degrees in longitude. This results in a matrix of 2165 x 4325 elements.")
                println(io)
                println(io, "The data is stored with negative latitude/longitude in the upper left corner and positive latitude/longitude in the lower right corner. This is the opposite of how the raw data is stored in the text files downloaded directly from the ITU database.")
                println(io)
                println(io, "This artifact was automatically generated using the script at the following URL:")
                println(io, "https://github.com/JuliaSatcomFramework/ItuRPropagation.jl/blob/ITU_artifacts/scripts/p1511.jl")
            end

            archive = ZipReader(read(zip_path))
            for name in zip_names(archive)
                if endswith(name, "zip")
                    @info "Processing zip file $name"
                    # We have to recursively process this zip and extract it in a subfolder
                    subarchive = ZipReader(zip_readentry(archive, name))
                    for subname in zip_names(subarchive)
                        if subname in ("TOPO.dat", "EGM2008.txt")
                            @info "Converting file $subname"
                            matsize = subname === "TOPO.dat" ? (2164, 4324) : (2165, 4325)
                            filecontent = zip_readentry(subarchive, subname)
                            binfile = joinpath(artifact_folder, replace(subname, r"txt|dat" => "bin"))
                            data = stack(parseline, eachline(IOBuffer(filecontent))) |> permutedims |> reverse # The data is permuted because stack creates a transposed matrix in this case. It is also reversed (U/D and L/R) because the data has positive lat and positive lon in the upper left corner, while for indexing we want the opposite.
                            size(data) == matsize || error("Unexpected size: $(size(data)) instead of $matsize for file $subname")
                            open(binfile, "w") do io
                                write(io, data)
                            end
                        else
                            @info "Copying file $subname"
                            open(joinpath(artifact_folder, subname), "w") do io
                                write(io, zip_readentry(subarchive, subname))
                            end
                        end
                    end
                else
                    # We don't know what to do with this file
                    error("Unknown file type: $name")
                end
            end
        end

    end


    asset_name = "p1511_data.tar.gz"
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

create_p1511_artifact()