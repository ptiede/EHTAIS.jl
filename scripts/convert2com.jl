using Pkg;Pkg.activate(@__DIR__)
using Pyehtim
using Comrade
using Comonicon
using Serialization

"""
Converts from uvfits to the serialized data format used by Comrade

# Arguments

- `uvfile`: The UVFITS file to convert
- `out`: The output file to save the data to

# Options
- `-s, --snrcut=<float>`: The signal to noise ratio cut to apply to the data
- `-u, --uvmin=<float>`:  Flag any data below this limit
- `-f, --fracnoise=<float>`: Add this fractional noise to the data
"""
@main function main(uvfile, out; snrcut=3.0, uvmin=0.1e9, fracnoise=0.01)
    obs = ehtim.obsdata.load_uvfits(uvfile)
    obsavg = scan_average(obs)
    obsavg = obsavg.flag_uvdist(uvmin)
    obsavg = obsavg.add_fractional_noise(fracnoise)
    data = extract_table(obsavg, LogClosureAmplitudes(;snrcut), ClosurePhases(;snrcut))
    serialize(out, data)
end
