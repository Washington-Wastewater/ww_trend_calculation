# This script contains simple metric functions for calculating 
# geometric mean, and others in the core dataset.
# WVAL functions are in a separate R file due to their length and complexity.

# load libraries
library(tidyverse)

## Function: geo_mean()
# - this is for generating 7-day rolling geometric means of normalized concentrations
geo_mean <- function(x){
  exp(mean(log(x), na.rm=TRUE))
}


