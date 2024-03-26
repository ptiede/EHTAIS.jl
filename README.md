# EHT Average Image Scoring

This is a set of scripts that will perform the EHT AIS procedure. Currently we only have the snapshot
fitting implemented natively in Julia.

To run the snapshot scoring we recommend that most users use the `scripts/parallel_snapshotfits.jl` 
script. To run this you will first need to install julia on your machine. I recommend doing this with [`juliaup`](https://github.com/JuliaLang/juliaup).

To get setup your environment then move to the script directory and then run

- Instantiate the environment `julia setup.jl`
- Convert your data to Comrade serialized format `julia convert_data.jl /path/to/eht/uvfits /path/to/eht/comrade/data`

Finally to run the scoring do
```
> readlink -f /path/to/grmhd/images > imfiles
> julia -p NCORES parallel_snapshotfits.jl "imfiles" "/path/to/grmhd/README.txt" "path/to/eht/comrade/data" "path/to/output/file"
```

Note that julia is a JIT language so there is some warmup time as the compiler runs. If this is really long please let me know and we can try to compile a special Julia image.