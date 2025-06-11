# needed libraries
library(magrittr)
library(dplyr)
library(tidyselect)

# load data
setwd("G:/Hangkai/Anttarctic Vegetation Dynamic/RF_trainning_data/cleared_data/")
df <- read.csv("merged_all_data_VPD_updated.csv") %>%
  # make an indicator variable for vegetation area
  dplyr::mutate(veg_binary = ifelse(vegetation_area_ratio == 0, 0, 1)) %>%
  # define pixel to be a unique site indicator
  dplyr::mutate(pixel = as.numeric(as.factor(paste0(latitude, ", ", longitude))))

# helper_functions
# percent zero vs positive pixels
zero_perc <- function(var) {
  group_counts <- table(var)
  total_count <- sum(group_counts)
  group_percentages <- (group_counts / total_count) * 100
  print(group_percentages)
}
# correlation coefficient check function
cor_check <- function(threshold = 0.7, df, var_list) {
  correlation_matrix <- df %>%
    dplyr::select(all_of(var_list)) %>%
    cor(use = "pairwise.complete.obs")
  correlation_tidy <- as.data.frame(as.table(correlation_matrix)) %>%
    dplyr::rename(Var1 = Var1, Var2 = Var2, Correlation = Freq) %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::filter(abs(Correlation) > threshold) %>%
    dplyr:: mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pair = paste0(pmin(Var1, Var2), "-", pmax(Var1, Var2))) %>% 
    dplyr::distinct(pair, .keep_all = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-pair) 
  # display the significant correlations
  print(correlation_tidy, n = Inf)
}

#----------------------------------------------------------------------------------------------------------------------
# stratify by subregion to get split of zeroes down to 50/50 or some other target
# ok this dataset appears to only be for one subregion right now, Hangkai is checking

veg_by_pixel <- df %>%
  dplyr::select(pixel, year, veg_binary) %>%
  # compute the number of years out of 22 where there is vegetation in each pixel
  dplyr::group_by(pixel) %>%
  dplyr::summarize(total_veg = sum(veg_binary)) %>%
  dplyr::ungroup()
hist(veg_by_pixel$total_veg, xlab = "years with vegetation", main = "Histogram of years with vegetation across pixels")

zero_perc(df$veg_binary)
df <- df %>%
  # merge this back to df by pixel and filter out pixels that only contain 1 year of vegetation
  dplyr::left_join(veg_by_pixel, by = "pixel") %>%
  dplyr::filter(total_veg > 1)
zero_perc(df$veg_binary)

#----------------------------------------------------------------------------------------------------------------------
# spatial and temporal decomposition

# some of these covariates are constant across time and we should not decompose them
no_decomp_vars <- c("Elevation",
                    "Slope",
                    "Aspect")      
# these covariates vary across space and time and we should decompose them
decomp_vars <- c("temperature_2m",
                 "icefree_area_ratio",
                 "uv_radiation",
                 "volumetric_soil_water",
                 "snowmelt",
                 "skin_temperature",
                 "runoff",
                 "vapour_pressure",
                 "solar_radiation",
                 "X10m_wind_speed",
                 "precipitation",
                 "vapor_pressure_deficit")
df <- df %>% 
  # center without scaling
  dplyr::mutate(across(all_of(decomp_vars), ~ scale(., scale = FALSE)[,1])) %>% 
  # compute the spatial component (mean by pixel across all years of centered variables)
  dplyr::group_by(pixel) %>%
  dplyr::mutate(dplyr::across(tidyselect::all_of(decomp_vars), ~ mean(., na.rm = TRUE), .names = "{.col}_spatial")) %>% 
  dplyr::ungroup() %>%
  # compute the temporal component (mean by year across all pixels of centered variables)
  dplyr::group_by(year) %>%
  dplyr::mutate(dplyr::across(tidyselect::all_of(decomp_vars), ~ mean(., na.rm = TRUE), .names = "{.col}_temporal")) %>% 
  dplyr::ungroup() %>%
  # compute residual for each site i and year j as centered variable value i,j - spatial mean i - temporal mean j
  dplyr::mutate(dplyr::across(tidyselect::all_of(decomp_vars), 
    ~ . - get(paste0(cur_column(), "_spatial")) - get(paste0(cur_column(), "_temporal")), .names = "{.col}_residual")) %>%
  dplyr::arrange(pixel, year)

#----------------------------------------------------------------------------------------------------------------------
# check correlations
# set the threshold for highlighting high correlations
cor_check(threshold = .7, df = df, var_list = c(decomp_vars, no_decomp_vars))
# now do this for the decomposed variables
decomp_spatial <- paste0(decomp_vars, "_spatial")
decomp_temporal <- paste0(decomp_vars, "_temporal")
decomp_residual <- paste0(decomp_vars, "_residual")
# combine all decomposed variable names
all_decomp_vars <- c(decomp_spatial, decomp_temporal, decomp_residual)
cor_check(threshold = .7, df = df, var_list = c(all_decomp_vars, no_decomp_vars))

# remove highly correlated variables and recheck
rf_vars <- c("Aspect", 
             "icefree_area_ratio",
             "temperature_2m", 
             "uv_radiation", 
             "volumetric_soil_water",
             "runoff",
             "X10m_wind_speed",
             "precipitation",
             "vapor_pressure_deficit")
cor_check(threshold = .7, df = df, var_list = rf_vars)
rf_vars <- c("temperature_2m", 
             "uv_radiation", 
             "icefree_area_ratio",
             "volumetric_soil_water",
             "runoff",
             "X10m_wind_speed",
             "precipitation",
             "vapor_pressure_deficit")
decomp_spatial <- paste0(rf_vars, "_spatial")
decomp_temporal <- paste0(rf_vars, "_temporal")
decomp_residual <- paste0(rf_vars, "_residual")
rf_vars_decomp <- c("Aspect", decomp_spatial, decomp_temporal, decomp_residual)
cor_check(threshold = .7, df = df, var_list = rf_vars)
rf_vars_decomp

#----------------------------------------------------------------------------------------------------------------------
# do a simple training test split
all_df <- df %>%
  dplyr::group_by(Subregions) %>%
  dplyr::sample_frac(size = 1.0, replace = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(c(Regions, Subregions, pixel, year, latitude, longitude, vegetation_area_ratio, all_of(rf_vars_decomp))) %>%
  dplyr::arrange(Regions, Subregions, pixel, year)
write.csv(all_df, "decomp_df.csv")