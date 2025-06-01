include(joinpath(@__DIR__, "common.jl"))

# This is the URL for direct download of the zip file containing the Part 1 of the ITU-R P.2145-0 recommendation (Annual Data)
url = "https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.2145Part01-0-202208-I!!ZIP-E.zip"

zip_path = joinpath(downloads_dir, "p2145_2022_Part1_annual.zip")
if !isfile(zip_path)
    Downloads.download(url, zip_path)
end

p1245_annual_hash = artifact_hash(artifact_toml, "p2145_annual")

if isnothing(p1245_annual_hash) || !artifact_exists(p1245_annual_hash)
    p2145_annual_hash = create_artifact() do artifact_folder
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
            println(io, "This artifact was automatically generated using the script")
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
