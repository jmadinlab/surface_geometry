# Surface geometry and biodiversity
source("R/functions.R")

# Scope (extent), scales of variation, and resolution (grain)
L <- 0.5 # Scope, 2 by 2 m reef patches
scl <- L / c(1, 2, 4, 8) # , 16)#, 32, 64, 128) # Scales, aim for 2 orders of magnitude
L0 <- min(scl) # Grain, resolution of processing ~ 6 cm

# Example surface (an 8x8m section of Horseshoe from Lizard Island)
output <- "noaa" # For housekeeping

# Load example geotif
# data <- raster("data/example/horseshoe.tif")
data <- raster("data/noaa/HAW_4263.tif")
plot(data)

rep <- 1
# Choose patch in which to calculate RDH (rugosity, fractal D and height range).
xb <- data@extent[1]
yb <- data@extent[3]

for (i in seq(0, 2, 0.5)) {
  for (j in seq(0, 0.5, 0.5)) {
    
    x0 <- xb + i
    y0 <- yb + j
    
    rect(x0, y0, x0+L, y0+L, border="white", lty=2)
    
    # Calulate height variation at different scales (scl) within patch, and save output (because a time-consuming step)
    example <- height_variation(write=TRUE, return=TRUE)
    
    rep <- rep + 1
  }
}


# Load the file if starting here:
example <- read.csv(paste0("output/", output, "/var_", names(data), "_0010.csv"), as.is=TRUE)

# Calculate rugosit, fractal dimension and height range (rdh function)
rdh(example)

