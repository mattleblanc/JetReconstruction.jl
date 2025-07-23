# JetReconstruction.jl Examples - GhostedArea

## `area_calculation.jl`

This example can be run as:

```sh
julia --project area_calculation.jl  --maxevents=100 --resolution=100  --algorithm=Kt --data-file=../../test/data/ghosted_area_example.hepmc.zst
```

This is an example of GhostedArea that reads from a HepMC3 file
and prints the areas of the inputted/reconstructed jets.