# Functions and packages

library(raster)
library(parallel) # For faster processing
library(rgdal)
#library(rgl)
library(plot3D)
library(mgcv)
library(sgeostat)
library(jpeg)
library(fields)

# Function for calculating variation in windows
fd_func <- function(x, y, s) {
  bx <- extent(cbind(c(x0 + x, y0 + y), c(x0 + x + s, y0 + y + s)))
  return(diff(range(getValues(crop(data, bx)), na.rm=TRUE)))
}

rescale <- function(x, x0, xm, n) {
  (x - x0)/(xm - x0)*n
}

count_filled_cells <- function(dat, xmin, xmax, ymin, ymax, zmin, zmax, ngrid) {
  ret <- table(ceiling(rescale(dat[,1], xmin, xmax, ngrid)), ceiling(rescale(dat[,2], ymin, ymax, ngrid)), ceiling(rescale(dat[,3], zmin, zmax, ngrid)))
  sum(ret > 0)
}

R_func <- function(H0, L0) {
  sqrt((H0^2) / (2 * L0^2) + 1)
}

D_func <- function(H, R, L, L0) {
  3 - log10(H / (sqrt(2) * L0 * sqrt(R^2 - 1))) / log10(L / L0)
}

HL0_func <- function(D, R, L, L0) {
  (3 - D) * log10(L/L0) + 0.5 * log10(R^2 - 1)
}

# Main functions

height_variation <- function(write=TRUE, return=FALSE) {
  # Fractal dimension, D
  temp <- data.frame()
  for (s in scl) {
    inc <- seq(0, L-s, s)
    x <- rep(inc, L/s)
    y <- rep(inc, each=L/s)
    # If you can't get multicore (parallel) package working, change "mcmapply" to "mapply" below:
    temp <- rbind(temp, data.frame(L0=s, x=x, y=y, H0=mcmapply(fd_func, x, y, s)))
  }
  # This variation method is time-consuming, and so save the result to avoid reprocessing if recalculting RDH
  if (write) {
    write.csv(temp, paste0("output/", output, "/var_", names(data), "_", sprintf("%04d", rep), ".csv"), row.names=FALSE)
  }
  print(paste0("Complete: ", names(data), "_", sprintf("%04d", rep)))
  # You can return the data and assign to variable if wish
  if (return) {
    return(temp)
  }
}

rdh <- function(hvar) {
  # log10 transform
  hvar$H0 <- log10(hvar$H0)
  hvar$L0 <- log10(hvar$L0)
  # plot(H0 ~ L0, hvar)
  # Mean of scales to avoid biased sampling at smaller scales
  hvar_m <- aggregate(H0 ~ L0, hvar, mean)
  # points(H0 ~ L0, hvar_m, col="red")
  
  # Find the height ranges at both ends of the scale
  H <- 10^hvar_m$H0[hvar_m$L0==log10(L)]
  H0 <- 10^hvar_m$H0[hvar_m$L0==log10(L0)]

  # Re-centering, probably unnessary
  hvar_m$H0 <- hvar_m$H0 - hvar_m$H0[hvar_m$L0==log10(L)]
  hvar_m$L0 <- hvar_m$L0 - log10(L)

  # Calculate slopes and minus from 3 (i.e., to get D from S)
  mod <- lm(H0 ~ L0, hvar_m)
  D <- 3 - coef(mod)[2]
  mod_ends <- lm(H0 ~ L0, hvar_m[c(1, nrow(hvar_m)),])
  D_ends <- 3 - coef(mod_ends)[2]

  # Calculate rugosity from theory
  R_theory <- R_func(H0, L0)		
  
  # Calculate rugosity from theory (integral)
  HL0 <- 10^hvar$H0[hvar$L0==min(hvar$L0)]
  R_theory2 <-  sum(R_func(HL0, L0)) / (2/L0)^2
  
  R <- NA
  # Optional test: calculating R using another method
  temp3 <- crop(data, extent(x0, x0 + L, y0, y0 + L))
  temp3 <- as(temp3, 'SpatialGridDataFrame')
  R <- surfaceArea(temp3) / L^2
  
  # Calcualte D from theory
  D_theory <- D_func(H, R_theory, L, L0)

  return(list(D=D, D_ends=D_ends, D_theory=D_theory, R=R, R_theory=R_theory, R_theory2=R_theory2, H=H))
}

