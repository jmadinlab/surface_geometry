# Data preparation

# Load geotifs
data <- raster("data/megaplot/trimodal3set_photoscan_2chunk_DEM_7mm_proj_clip.tif")
mids <- readOGR("data/megaplot/trimodal_patch_grid_midpnts.shp")
anno <- readOGR("data/megaplot/trimodal_ann_aligned_cleaned.shp", stringsAsFactors=FALSE)

# Load reef records and megaplot (Trimodal large) dataset
records <- read.csv("output/records.csv", as.is=TRUE)
megaplot <- read.csv("output/megaplot.csv", as.is=TRUE)

# Add empty columns for records
records$spp <- NA
records$abd <- NA
records$pie <- NA

# Merge datasets
dat <- rbind(records, megaplot)

# Surface descriptors and transformations
dat$R2_log10 <- log10(dat$R_theory^2 - 1)
dat$R2_log10_sq <- dat$R2_log10^2

dat$HL0_log10 <- log10(dat$H / (L0 * sqrt(2)))
dat$HL0_log10_sq <- dat$HL0_log10^2

dat$D_theory_sq <- dat$D_theory^2

write.csv(dat, "output/master_20200709.csv", row.names = FALSE)
