# EHT Average Image Scoring

This is a set of scripts that will perform the EHT AIS procedure. Currently we only have the snapshot
fitting implemented natively in Julia.

To run the snapshot scoring we recommend that most users use the `scripts/parallel_snapshotfits.jl` 
script. To run this you just need to call

```
> cd scripts 
> readlink -f /path/to/grmhd/images > imfiles
> julia -p NCORES parallel_snapshotfits.jl "imfiles" "/path/to/grmhd/README.txt" "path/to/eht/uvfits/data" "path/to/output/file"
```

The first time calling this will be a little slow because Julia will be installing all the dependencies. This should speed up considerably after that. On my machine with 16 cores I am able to analyze 3000 GRMHD snapshot with the EHT 2017 data in ~ 5min. If you experience is a lot slower please let me know.
