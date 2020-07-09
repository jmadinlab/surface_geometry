# Reef record test
source("R/functions.R")
output <- "record"

# Reef records
files <- dir("data/records")
files <- unique(gsub("\\.tfw|\\.tif", "", files))

# Scope (extent), scales of variation, and resolution (grain)
L <- 2 # Scope
scl <- L / c(1, 2, 4, 8, 16, 32) # Scales, aim for 2 orders magnitude
L0 <- min(scl) # Grain

store <- data.frame()

for (rec in files) {
  # Load geotif for reef record
  data <- raster(paste0("data/records/", rec, ".tif"))
  # Get lower corner of 8x8m bounding box; alternatively, prescribe these values.
  xb <- mean(data@extent[1:2]) - 4
  yb <- mean(data@extent[3:4]) - 4

  # Iterate through  2x2m quadrats (reps = 16) in reef record
  rep <- 1
  for (i in c(0, 2, 4, 6)) {
    for (j in c(0, 2, 4, 6)) {
      x0 <- xb + i
      y0 <- yb + j
      fname <- paste0("output/records/var_", names(data), "_", sprintf("%04d", rep), ".csv")
      if (file.exists(fname)) {
        temp <- read.csv(fname, as.is=TRUE)
      } else {
        temp <- height_variation(write=TRUE, return=TRUE)
      }
      store <- rbind(store, data.frame(rec=rec, rep=rep, rdh(temp)))
      rep <- rep + 1
    }
  }
}

store$site[store$rec=="rr201611_004_Osprey_dem_low"] <- "Osprey"
store$site[store$rec=="rr201611_007_CooksPath_dem_low"] <- "Cooks Path"
store$site[store$rec=="rr201611_018_Resort_dem_low"] <- "Resort"
store$site[store$rec=="rr201611_023_NorthReef03_dem_low"] <- "Mermaid Cove"
store$site[store$rec=="rr201611_033_ConerBeach_dem_low"] <- "Corner Beach"
store$site[store$rec=="rr201611_037_southeast_dem_low"] <- "Southeast"
store$site[store$rec=="rr201611_039_EasterPoint_dem_low"] <- "Easter Point"
store$site[store$rec=="rr201611_040_NoMansLand_dem_low"] <- "No Mans Land"
store$site[store$rec=="rr201611_041_NorthOfParadise_dem_low"] <- "North of Paradise"
store$site[store$rec=="rr201611_042_GnarlyTree_dem_low"] <- "Gnarly Tree"
store$site[store$rec=="rr201611_046_NorthReef02_dem_low"] <- "North Reef 2"
store$site[store$rec=="rr201611_050_horsehoe_DEM_low"] <- "Horseshoe"
store$site[store$rec=="rr201611_051_Vickis_dem_low"] <- "Vickis"
store$site[store$rec=="rr201611_052_SouthIsland_dem_low"] <- "South Island"
store$site[store$rec=="rr201611_053_Trimodal_dem_low"] <- "Trimodal"
store$site[store$rec=="rr201611_054_lagoon02_dem_low"] <- "Lagoon 2"
store$site[store$rec=="rr201611_055_Lagoon01_dem_low"] <- "Lagoon 1"
store$site[store$rec=="rr201611_055_LizardHead_dem_low"] <- "Lizard Head"
store$site[store$rec=="rr201611_056_TurtleBeach_dem_low"] <- "Turtle Beach"
store$site[store$rec=="rr201611_045_WashingMachine_dem_low"] <- "Washing Machine"
store$site[store$rec=="rr201611_049_NorthReef01_dem_low"] <- "North Reef 1"

write.csv(store, paste0("output/records.csv"), row.names=FALSE)
