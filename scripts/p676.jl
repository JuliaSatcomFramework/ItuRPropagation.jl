include(joinpath(@__DIR__, "common.jl"))

function create_p676_artifact()
    artifact_name = "p676"
    _artifact_hash = Artifacts.artifact_hash(artifact_name, artifact_toml)

    if isnothing(_artifact_hash) || !Artifacts.artifact_exists(_artifact_hash)
        _artifact_hash = Artifacts.create_artifact() do artifact_folder
            part1_url = "https://www.itu.int/dms_pub/itu-r/oth/11/01/R11010000020001TXTM.txt"
            part2_url = "https://www.itu.int/dms_pub/itu-r/oth/11/01/R11010000020002TXTM.txt"

            part1_txt = joinpath(artifact_folder, "R11010000020001TXTM.txt")
            part2_txt = joinpath(artifact_folder, "R11010000020002TXTM.txt")

            Downloads.download(part1_url, part1_txt)
            Downloads.download(part2_url, part2_txt)

            fvrange = range(1, 350; step = 0.5) |> collect
            forange = let 
                v = copy(fvrange)
                idx = searchsortedfirst(v, 118.75)
                insert!(v, idx, 118.75)
                v
            end

            parts = [
                (; file = part1_txt, name = "Part1", frange = forange),
                (; file = part2_txt, name = "Part2", frange = fvrange)
            ]

            # Add a README
            open(joinpath(artifact_folder, "README"), "w") do io
                println(io, "This folder contains data for the ITU-R P.676-13 recommendation.")
                println(io)
                println(io, "It was automatically generated from the TXT files of the ITU-R P.676-13 recommendation available at the following URL:")
                println(io, "https://www.itu.int/oth/R1101000002/en")
                println(io)
                println(io, "This artifact contains the raw txt data as well as the processed binary data obtained by parsing the txt files")
                println(io)
                println(io, "For both Part1 and Part2, the data is stored as a n x 5 column matrix where the columns are:")
                println(io, "1. Frequency in MHz")
                println(io, "2. aₒ (or aᵥ for Part 2)")
                println(io, "3. bₒ (or bᵥ for Part 2)")
                println(io, "4. cₒ (or cᵥ for Part 2)")
                println(io, "5. dₒ (or dᵥ for Part 2)")
                println(io)
                println(io, "For both parts, the input frequency goes from 1 MHz to 350 MHz in steps of 0.5 MHz. For the Part 1, an additional line at 118.75 MHz is added, making the size different between Part1 and Part2 by 1")
                println(io)
                println(io, "This artifact was automatically generated using the script at the following URL:")
                println(io, "https://github.com/JuliaSatcomFramework/ItuRPropagation.jl/blob/ITU_artifacts/scripts/p676.jl")
            end

            for part in parts
                @info "Processing $(part.name) raw file ($(basename(part.file)))"
                data = parsematrix(part.file)
                binfile = joinpath(artifact_folder, part.name * ".bin")
                open(binfile, "w") do io
                    write(io, data)
                end
            end
        end
    end


    asset_name = "p676_data.tar.gz"
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

create_p676_artifact()