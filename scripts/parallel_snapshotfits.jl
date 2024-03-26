using Distributed
@everywhere begin
    const filedir = @__DIR__
end

using Pkg; Pkg.activate(filedir)
using Distributed
using DelimitedFiles
using Serialization

@everywhere begin
    using Pkg;Pkg.activate(filedir)
    @info "Done Activating"
end

using EHTAIS
@everywhere using EHTAIS


@info "Done Loading Modules"

using Comrade
using Comonicon
using DataFrames
using CSV

function loaddir(file)
    open(file) do f
        return readlines(f)
    end
end

function load_obs(uvname)
    return deserialize(uvname)
end

function load_data(uvname, ::Type{<:Closures})
    obs = load_obs(uvname)
    return Closures(obs...)
end


"""
Runs snapshot fitting on the list of files passed

# Arguments

- `imfile`: The file containing the paths to all the GRMHD we will to analyze
- `readme`: The file with the GRMHD readme describing the dimensions
- `uvfile`: The path to the data we are going to fit, this must be serialized data
            Please run `convert2com.jl` to convert the uvfits data to the correct format.
- `outfile`: The name of the output file where the results are saved.

# Options

- `-f, --fevals=<int>`: The number of evaluations of the loglikelihood.
- `-y, --year=<string>`: The year of the data. Options are 2017 and 2018

"""
@main function main(imfile::String, readme::String, uvfile::String,
                    outfile::String="snapshot_fitresults.csv";
                    fevals::Int=10_000,
                    )

    @info "Image files path: $(imfile)"
    @info "Readme path: $(readme)"
    @info "Fitting data: $(uvfile)"
    @info "outputting results to $(outfile)"

    data = load_data(uvfile, Closures)
    imfiles = loaddir(imfile)
    res = pmap(imfiles) do f
        @info "Fitting $f"
        return snapshot_fit(f, readme, data; fevals)
    end
    CSV.write(outfile, DataFrame(res))
end
