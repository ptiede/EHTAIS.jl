@info "Initiating Setup"

using Distributed

if nworkers() > 1
    @error "Please run setup with only a single process"
end

using Pkg; Pkg.activate(@__DIR__)


Pkg.add(path="../")
Pkg.instantiate()
Pkg.precompile()

using EHTAIS

using PreferenceTools
pkg"preference add CondaPkg offline=true"

@info "Finished setup have ready to run snapshot fits"
