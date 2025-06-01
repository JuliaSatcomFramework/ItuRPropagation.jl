include(joinpath(@__DIR__, "common.jl"))

function create_p2145_annual_artifact()
    p2145_annual_hash = artifact_hash("p2145_annual", artifact_toml)

    if isnothing(p2145_annual_hash) || !artifact_exists(p2145_annual_hash)
        p2145_annual_hash = create_artifact() do artifact_folder
            # This is the URL for direct download of the zip file containing the Part 1 of the ITU-R P.2145-0 recommendation (Annual Data)
            url = "https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.2145Part01-0-202208-I!!ZIP-E.zip"

            zip_path = joinpath(downloads_dir, "p2145_2022_Part1_annual.zip")
            if !isfile(zip_path)
                @info "Downloading raw zip file for ITU-R P.2145-0 Part 1 (Annual Data) from ITU Database"
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
                println(io, "This folder contains data for the ITU-R P.2145-0 recommendation (Annual Data) directly in binary format.")
                println(io)
                println(io, "It was automatically generated from the TXT files of the Part 1 annex of the ITU-R P.2145-0 recommendation available at the following URL:")
                println(io)
                println(io, url)
                println(io, "The matrices stored into each of the binary files corresponds to a lat/lon grid with $(latres)° x $(lonres)° resolution and have a size of $(matsize) elements.")
                println(io)
                println(io, "This artifact was automatically generated using the script at the following URL:")
                println(io, "https://github.com/JuliaSatcomFramework/ItuRPropagation.jl/blob/ITU_artifacts/scripts/p2145.jl")
            end

            archive = ZipReader(read(zip_path))
            for name in zip_names(archive)
                if endswith(lowercase(name), "txt")
                    # We copy this to the artifact folder
                    open(joinpath(artifact_folder, name), "w") do io
                        write(io, zip_readentry(archive, name))
                    end
                elseif endswith(name, "zip")
                    @info "Processing zip file $name"
                    # We have to recursively process this zip and extract it in a subfolder
                    foldername = joinpath(artifact_folder, replace(name, ".zip" => ""))
                    isdir(foldername) || mkdir(foldername)
                    subarchive = ZipReader(zip_readentry(archive, name))
                    for subname in zip_names(subarchive)
                        @info "Converting file $subname"
                        endswith(subname, "TXT") || error("Unexpected file type: $subname")
                        filecontent = zip_readentry(subarchive, subname)
                        binfile = joinpath(foldername, replace(subname, ".TXT" => ".bin"))
                        data = readdlm(filecontent, ' ')
                        size(data) == matsize || error("Unexpected size: $(size(data)) instead of $matsize for file $subname")
                        open(binfile, "w") do io
                            write(io, data)
                        end
                    end
                else
                    # We don't know what to do with this file
                    error("Unknown file type: $name")
                end
            end
        end
        bind_artifact!(artifact_toml, "p2145_annual", p2145_annual_hash)
    end

    tarball_path = joinpath(assets_dir, "p2145_2022_Part1_annual.tar.gz")
    tarball_sha = Pkg.Artifacts.archive_artifact(p2145_annual_hash, tarball_path)
    return p2145_annual_hash

end

create_p2145_annual_artifact()
