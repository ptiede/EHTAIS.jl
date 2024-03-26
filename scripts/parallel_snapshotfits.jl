using Distributed
@everywhere begin
    const filedir = @__DIR__
end


using Pkg; Pkg.activate(filedir)
using Distributed
using DelimitedFiles

@everywhere begin
    using Pkg;Pkg.activate(filedir)
end
# using EHTAIS
# @everywhere 1 using Pyehtim
# @everywhere using EHTAIS

using Comrade
using Comonicon
using DataFrames
using CSV

function loaddir(file)
    open(file) do f
        return readlines(f)
    end
end

function load_obs(uvname; cutzbl=true, fracnoise=0.01)
    obs = ehtim.obsdata.load_uvfits(uvname)
    obsavg = scan_average(obs)
    if cutzbl
        obsavg = obsavg.flag_uvdist(uv_min=0.1e9)
    end

    if fracnoise > 0.0
        obsavg = obsavg.add_fractional_noise(fracnoise)
    end
    return obsavg
end

function load_data(uvname, ::Type{<:Closures}; snrcut=3.0, cutzbl=true, fracnoise=00.01)
    obs = load_obs(uvname; cutzbl, fracnoise)
    return Closures(extract_table(obs, LogClosureAmplitudes(;snrcut), ClosurePhases(;snrcut))...)
end


"""
Runs snapshot fitting on the list of files passed

# Arguments

- `imfile`: The file containing the paths to all the GRMHD we will to analyze
- `readme`: The file with the GRMHD readme describing the dimensions
- `uvfile`: The path to the data we are going to fit
- `outfile`: The name of the output file where the results are saved.

# Options

- `-f, --fevals=<int>`: The number of evaluations of the loglikelihood.
- `-y, --year=<string>`: The year of the data. Options are 2017 and 2018

"""
@main function main(imfile::String, readme::String, uvfile::String,
                    outfile::String="snapshot_fitresults.csv";
                    fevals::Int=250_000,
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
