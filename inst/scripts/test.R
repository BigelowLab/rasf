library(raster)
library(terra)
library(rasf)

R <- volcano_raster(klass = "BasicRaster")
T <- volcano_raster(klass = "SpatRaster")


example_filepairs() %>%
 paired_quality_scores() %>%
  paired_ee_threshold(sample_names  = c("foo", "bar"), filename = "~/my_ee_stuff.csv")