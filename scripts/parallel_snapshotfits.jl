
@everywhere begin
    const filedir = @__DIR__
end


using Pkg; Pkg.activate(filedir)
using Distributed
using DelimitedFiles

Pkg.add(path="../")
# Pkg.add(["Comonicon", "DataFrames", "CSV"])
Pkg.instantiate()
Pkg.precompile()


@everywhere begin
    using Pkg;Pkg.activate(filedir)
end
using CondaPkg
ENV["JULIA_CONDAPKG_OFFLINE"] = "yes"
using EHTAIS
# Now turn off CondaPkg
@everywhere begin
    using CondaPkg
    ENV["JULIA_CONDAPKG_OFFLINE"] = "yes"
end
@everywhere using EHTAIS

using Comonicon
using DataFrames
using CSV

function loaddir(file)
    string.(reshape(readdlm(file), :))
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

# Flags

- `-a, --amp`: A flag that we should fit amp+cp and not closures. Warning amps are 10x slower.
"""
@main function main(imfile::String, readme::String, uvfile::String,
                    outfile::String="snapshot_fitresults.csv";
                    fevals::Int=250_000, amp::Bool=false)

    @info "Image files path: $(imfile)"
    @info "Readme path: $(readme)"
    @info "Fitting data: $(uvfile)"
    @info "Fit amps: $amp"
    @info "outputting results to $(outfile)"

    if amp
        data = load_data(uvfile, AmpCP)
        distamp = station_tuple(data.amp, Normal(0.0, 0.1); LM=Normal(0.0, 1.0))
    else
        data = load_data(uvfile, Closures)
        distamp = nothing
    end
    imfiles = loaddir(imfile)
    res = pmap(imfiles) do f
        amp && return snapshot_fit(f, readme, data, distamp; lbfgs=true, fevals)
        return snapshot_fit(f, readme, data, distamp; lbfgs=false, fevals)
    end
    CSV.write(outfile, DataFrame(res))
end
