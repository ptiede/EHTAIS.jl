module EHTAIS

using Comrade
using ComradeOptimization
using DataFrames
using CSV
using DelimitedFiles
using Distributions
using DistributionsAD
using Zygote
using OptimizationMetaheuristics
using OptimizationOptimJL
using ForwardDiff

export load_image, load_readme, load_data,
       Closures, AmpCP, Vis, station_tuple,
       snapshot_fit

abstract type DataProds end
keys(c::DataProds) = fieldnames(typeof(c))
values(c::DataProds) = getfield.(Ref(c), fieldnames(typeof(c)))

struct Closures{A,C} <: DataProds
    lcamp::A
    cphase::C
end

struct AmpCP{A, C} <: DataProds
    amp::A
    cphase::C
end

struct Vis{V} <: DataProds
    vis::V
end

function load_obs(uvname; cutzbl=true, fracnoise=0.01)
    obs = load_ehtim_uvfits(uvname)
    obsavg = scan_average(obs)
    if cutzbl
        obsavg = obsavg.flag_uvdist(uv_min=0.1e9)
    end

    if fracnoise > 0.0
        obsavg = obsavg.add_fractional_noise(fracnoise)
    end
    return obsavg
end

"""
    load_data(uvname, dataprod::Type{<:Closures}; snrcut=3.0, cutzbl=true, fracnoise=0.01)
    load_data(uvname, dataprod::Type{<:AmpCP};    snrcut=3.0, cutzbl=true, fracnoise=0.01)
    load_data(uvname, dataprod::Type{<:Vis}; s    nrcut=0.0, cutzbl=true, fracnoise=0.01)

Loads the data from `uvname` and extracts the data products `dataprod`.cpus

# Arguments
  - `uvname`: The file path for the uvfits data you would like to fit
  - `dataprod`: The data products you would like to extract, see `subtypes(EHTAIS.DataProds)` for a list
  - `snrcut`: The SNR cut use for the data products
  - `cutzbl`: If true we flag any baselines whose uv distance < 0.1e9
  - `fracnoise`: How much multiplicative fractional error to add to the data
"""
function load_data end

function load_data(uvname, ::Type{<:Closures}; snrcut=3.0, cutzbl=true, fracnoise=00.01)
    obs = load_obs(uvname; cutzbl, fracnoise)
    return Closures(extract_lcamp(obs; snrcut), extract_cphase(obs; snrcut))
end

function load_data(uvname, ::Type{<:AmpCP}; snrcut=3.0, cutzbl=true, fracnoise=00.01)
    obs = load_obs(uvname; cutzbl, fracnoise)
    return AmpCP(extract_amp(obs; debias=true), extract_cphase(obs; snrcut))
end

function load_data(uvname, ::Type{<:Vis}; snrcut=0.0, cutzbl=true, fracnoise=00.01)
    obs = load_obs(uvname; cutzbl, fracnoise)
    return Vis(extract_vis(obs; snrcut, debias=true))
end






"""
    load_readme(readme::String)

Loads the README used for traditional AIS. The output is a `NamedTuple` with
field of view in the `x/y` direction (`fovx/fovy`), the number of pixels (`nx/ny`) in the
`x/y` directions, and the simulations M/D in radians.
"""
function load_readme(readme::String)
    open(readme, "r") do io
        readline(io)
        readline(io)
        readline(io)
        lx = readline(io)
        fovx = split(last(split(lx, "="))) |> first |> x->parse(Float64, x) |> μas2rad
        ly = readline(io)
        fovy = split(last(split(ly, "="))) |> first |> x->parse(Float64, x) |> μas2rad
        lnx = readline(io)
        nx = parse(Int, last(split(lnx)))
        lny = readline(io)
        ny = parse(Int, last(split(lny)))

        # spin
        readline(io)

        # mass
        lma = readline(io)
        m = parse(Float64, split(lma)[3])

        # distance
        ld = readline(io)
        d = parse(Float64, split(ld)[3])

        mod = 6.67430e-11*m*1.9891e30/(3.086e16*d*(2.99792458e8)^2)

        return (;fovx, fovy, nx, ny, mod)
    end
end


"""
    load_image(file, readme)

Loads the simulations
"""
function load_image(file, readme)
    r = load_readme(readme)
    img = @view readdlm(file, dims=(r.nx*r.ny, 7))[:,4]

    grid = imagepixels(r.fovx, r.fovy, r.nx, r.ny)
    return IntensityMap(reshape(img, r.nx, r.ny)[end:-1:begin, :]./sum(img), grid), r.mod
end

function score_closures(θ, metadata)
    (;mod, pa,) = θ
    (;mimg, mod0) = metadata
    m = modify(mimg, Stretch(mod/mod0, mod/mod0), Rotate(pa))
    #jT = jonesStokes(exp.(lgamp), gcache)
    return m#JonesModel(jT, m)
end

function create_post(mimg::Comrade.ModelImage, mod0::Float64, data::Closures, distamp=nothing)
    metadata = (;mimg, mod0)
    lklhd = RadioLikelihood(score_closures, metadata, values(data)...)
    prior = (
                mod   = Uniform(μas2rad(0.1), μas2rad(6.0)),
                pa    = Uniform(-π/2, 3π/2),
            )
    post = Posterior(lklhd, prior)
    return post
end


function score_ampcp(θ, metadata)
    (;mod, pa, flux, lgamp) = θ
    (;mimg, mod0, gcache) = metadata
    m = modify(mimg, Stretch(mod/mod0, mod/mod0), Rotate(pa), Renormalize(flux))
    jT = jonesStokes(exp.(lgamp), gcache)
    return JonesModel(jT, m)
end

function create_post(mimg::Comrade.ModelImage, mod0::Float64, data::AmpCP, distamp)
    gcache = jonescache(data.amp, ScanSeg())
    metadata = (;mimg, mod0, gcache)
    lklhd = RadioLikelihood(score_ampcp, metadata, values(data)...)
    gcache = jonescache(data.amp, ScanSeg())
    prior = (
                mod   = Uniform(μas2rad(0.1), μas2rad(6.0)),
                pa    = Uniform(-π/2, 3π/2),
                flux  = Uniform(0.1, 3.0),
                lgamp = CalPrior(distamp, gcache)
            )
    post = Posterior(lklhd, prior)
    return post
end




"""
    snapshot_fit(img::IntensityMap, mod0::Float64, data::DataProds, distamp=nothing; fevals=250_000, lbfgs=true)
    snapshot_fit(imgname::String, readme::String, data::DataProds, distamp=nothing; fevals=250_000, lbfgs=true)

Computes a snapshot fit required for `ais`.

# Arguments
  - `img/imgname`: If `img::IntensityMap` this should be the GRMHD snapshot you want to fit.
                   If `imgname::String` this is the path to the file containing the GRMHD snapshot.
  - `mod0/readme`: If `mod0:Float64` this is the intrinsic M/D of the simulation
  - `data`: The data products you want to fit. See [`load_data`](@ref) to see more information about
            how to load data
  - `distamp`: The station priors for the log-gain amplitudes. This is only required when fitting `AmpCP` or `Vis`
  - `fevals`: The maximum number of posterior evaluations for the evoluationary optimizer
  - `lbfgs`: After use the evolutionary optimizer, use LBFGS to zoom to the peak.

# Example
```julia
# For closures
data = load_data(uvname, Closures)
res = snapshot_fit(imgname, readme, data; lbfgs=false)

# For amplitudes + closure phases
data = load_data(uvname, AmpCP)
distamp = station_tuple(data.amp, Normal(0.0, 0.1); LM = Normal(0.0, 1.0))
res = snapshot_fit(imgname, readme, data, distamp; lbfgs=true)
```
"""
function snapshot_fit(img::IntensityMap, mod0::Float64, data::DataProds, distamp=nothing; fevals=250_000, lbfgs=true)
    mimg = modelimage(ContinuousImage(img, DeltaPulse()), FFTAlg())
    post = create_post(mimg, mod0, data, distamp)

    tpost = asflat(post)
    ndim = dimension(tpost)

    fpost = OptimizationFunction(tpost, Optimization.AutoForwardDiff())
    prob0 = OptimizationProblem(fpost, rand(ndim) .- 0.5, nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))
    sol0 = solve(prob0, ECA(N=500, options=Metaheuristics.Options(f_calls_limit=fevals, f_tol=1e-3)))
    if lbfgs
        prob = OptimizationProblem(fpost, sol0.u, nothing)
        sol = solve(prob, LBFGS(); maxiters=2_000, g_tol=1e-1, f_tol=1e-7)
    else
        sol = sol0
    end

    score = post.lklhd.model
    xopt = Comrade.transform(tpost, sol.u)
    mopt = score(xopt)
    r2 = map(x->chi2(mopt, x)/(length(x) - ndim), values(data))
    r2data = NamedTuple{map(x->Symbol(:chi2_, x), keys(data))}(r2)
    rchi2 = chi2(mopt, values(data)...)/(sum(length, values(data)) - ndim)
    chi2data = merge(r2data, (chi2 = rchi2,))
    @info "chi2 $(rchi2)"
    return merge(xopt, chi2data)
end



function snapshot_fit(fname::String, readme::String, data::DataProds, distamp=nothing; fevals=250_000, lbfgs=true)
    img, mod = load_image(fname, readme)
    res = snapshot_fit(img, mod, data, distamp; fevals, lbfgs)
    return merge(res, (image_filename = fname,))
end


end # module EHTAIS
