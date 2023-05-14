@info "Initiating Setup"

using Distributed

if nworkers() > 1
    @error "Please run setup with only a single process"
end

using Pkg; Pkg.activate(@__DIR__)
using Distributed
using DelimitedFiles


Pkg.add(path="../")
Pkg.instantiate()
Pkg.precompile()

using CondaPkg
ENV["JULIA_CONDAPKG_OFFLINE"] = "yes"
using EHTAIS

ENV["JULIA_CONDAPKG_OFFLINE"] = "yes"

@info "Finished setup have ready to run snapshot fits"
