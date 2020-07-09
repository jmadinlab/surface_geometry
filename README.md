# A geometric basis for surface habitat complexity and biodiversity

The `analysis.R` code will reproduce the analyses and figures for the research posted currently at [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.02.03.929521v1). The paper has been provisionally accepted at *Nature Ecology & Evolution*.

The `example.R` code contains a generic example of how the code might be used to calculate the surface descriptors for from any DEM. There are essentially three steps:

1. Select the size of the patch for your calculations (`L`)
2. Select the scales across which fractal dimension will be calculated (`scl`), the smallest scale automatically is the resolution (`L0`)
3. Load a DEM file as `data` (geotif)
4. Pick the bottom left point of the patch (`x0` and `y0`)
5. Run the `height_variation()` function, which is essentially a wrapper function (see `R/functions.R`). This function only requires the variables mentioned so far (`L`, `scl`, `L0`, `data`, `x0` and `y0`), so make sure there're assigned.  The output of the `height_variation()` function is the DEM height range the prescribed scales. For the largest scale `H` (or 2x2m) there is only one value for the entire patch, and the next scale down (1x1m) there are four values, and so on. Best to assign the output to a variable and save it somewhere (see the Example), because the smallest scales are very time-consuming.
6. Calculate the surface descriptor metrics using the `rdh()` function. This function only requires the output from `height_variation()`. `rdh()` sends back a list with several metrics:

Variable | Description
--- | ---
D | Fractal dimension from model fit
D_ends | Fractal dimension only considering the largest (L) and smallest (L0) scale
D_theory | Fractal dimension calculated from theory (i.e., from R and H)
R | Surface rugosity calculated using `surfaceArea` function in R
R_theory | Surface rugosity calculated from theory ()
H | The height range (or height range at L)

### ToDo

- Apply to uses with 3D mesh
- Look at issues with D when large drop-offs in DEM (can cause D to go below 2)
