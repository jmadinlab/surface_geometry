# Process the megaplot (Trimodal large plot) data
source("R/functions.R")
output <- "megaplot"

# Load geotif
data <- raster("data/megaplot/trimodal3set_photoscan_2chunk_DEM_7mm_proj_clip.tif")
# These are midpoints for 2x2m squares, which are our "patch"-level samples
mids <- readOGR("data/megaplot/trimodal_patch_grid_midpnts.shp")
# These are the annotated coral colonies, ID'ed to species level
anno <- readOGR("data/megaplot/trimodal_ann_aligned_cleaned.shp")

# Scope (extent), scales of variation, and resolution (grain)
L <- 2 # Scope
scl <- L / c(1, 2, 4, 8, 16, 32) # Scales, aim for 2 orders magnitude
L0 <- min(scl) # Grain

# plot(data)
# rect(mids$x-1, mids$y-1, mids$x+1, mids$y+1, border="white")
# points(anno, col=rgb(1, 0, 0, 0.2), cex=0.1)

store <- data.frame()

for (mid in 1:length(mids)) {
  rep <- mid
  # Get lower corner of 2x2m bounding box; alternatively, prescribe these values some other way.
  x0 <- mids@data$x[mid] - L/2
  y0 <- mids@data$y[mid] - L/2

  fname <- paste0("output/megaplot/var_", names(data), "_", sprintf("%04d", rep), ".csv")
  if (file.exists(fname)) {
    temp <- read.csv(fname, as.is=TRUE)
  } else {
    temp <- height_variation(write=TRUE, return=TRUE)
  }

  tax <- crop(anno, extent(x0, x0 + 2, y0, y0 + 2))
  # points(tax, pch=20, cex=0.3)
  spp <- length(unique(tax$Species))
  abd <- length(tax$Species)
  pie <- 0
  if (abd > 0) {
    pie <- 1 - sum((table(tax$Species) / abd)^2)
  }
  
  # Calculate rugosit, fractal dimension and height range (rdh function)
  store <- rbind(store, data.frame(rec=names(data), rep=mid, rdh(temp), spp=spp, abd=abd, pie=pie, site="megaplot"))
  
}

write.csv(store, "output/megaplot.csv", row.names=FALSE)
