#This is the collaborative NWSS CoE translation of CDC's python code to calculate wastewater viral activity level (WVAL) for SARS-CoV-2, Influenza, and RSV with the NEW updated methodology published in August 2025. 

#WVAL categorizes viral concentration into minimal, low, moderate, high, and very high to help indicate the risk of infection. 

#From CDC's methodology: 

#Calculating the Wastewater Viral Activity Level

#Major Changes: 
#- No longer using normalized data
#- Baseline lookback period is 24 months for all respiratory targets (SARS-CoV-2, Flu A, RSV)
#- Requires 8 weeks of data 
#- Biannual recalculation dates for SARS-CoV-2 changed to April 1 and October 1
#- Updated WVAL level cut off points

#Baseline Calculation:
#-For each combination of site, data submitter, PCR target, lab methods, and normalization method, a baseline is established. The “baseline” is the 10th percentile of the log-transformed concentration data within a specific time frame. Details on the baseline calculation by pathogen are below:

#SARS-CoV-2
#-For site and method combinations (as listed above) with over 6 months of data, baselines are re-calculated every six calendar months (April 1st and October 1st) using the past 24 months of data.
#For sites and method combinations with less than 6 months of data, baselines are computed every time there is a new sample until reaching six months, after which they remain unchanged until the next April 1st or October 1st, at which time baselines are re-calculated.

#Influenza A and RSV
#-For site and method combinations (as listed above) with over 12 months of data, baselines are re-calculated every August 1st using all available data in the previous 24 months.
#For sites and method combinations with less than 12 months of data, baselines are computed weekly until reaching twelve months, after which they remain unchanged until the next August 1st, at which time baselines are re-calculated.

#-The standard deviation for each site and method combination is calculated using the same time frame as the baseline.

#Wastewater Viral Activity Level Calculation:
#-The number of standard deviations that each log-transformed concentration value deviates from the baseline (positive if above, negative if below) is calculated.
#This value (x) is then converted back to a linear scale (by calculating e^x) to form the Wastewater Viral Activity Level for the site and method combination.
#The Wastewater Viral Activity Levels from a site are averaged by week for all figures.


#This reflects the updated cut points: 
#Wastewater Viral Activity Level Categories for SARS-CoV-2 (https://www.cdc.gov/nwss/about-data.html#data-method, Feb 2025): 
#The current Wastewater Viral Activity Level for each state and territory is categorized into very low, low, moderate, high, or very high as follows:
#SARS-CoV-2:
#Very Low: Up to 2
#Low: Greater than 2 and up to 3.4	
#Moderate: Greater than 3.4 and up to 5.3	
#High: Greater than 5.3 and up to 7.8	
#Very High: Greater than 7.8

#Influenza A:
#Very Low: Up to 2.7
#Low: Greater than 2.7 and up to 6.2	
#Moderate: Greater than 6.2 and up to 11.2	
#High: Greater than 11.2 and up to 17.6	
#Very High: Greater than 17.6


#RSV:
#Very Low: Up to 2.5	
#Low: Greater than 2.5 and up to 5.2	
#Moderate: Greater than 5.2 and up to 8	
#High: Greater than 8 and up to 11	
#Very High: Greater than 11


#Aggregation for National, Regional, and State Levels:
#-We calculate the median Wastewater Viral Activity Levels among sites at national, regional, and state levels, excluding data from site/method combinations with less than 8 weeks of data for SARS-CoV-2, Influenza A, and RSV.

#Data Inclusion Criteria – SARS-COV-2, Influenza A, RSV Wastewater Viral Activity Level
#SARS-CoV-2: New wastewater sampling sites, or sites with a substantial change in laboratory methods are included in national, regional, state, or territorial median values once there are at least 8 weeks of samples reported for that location.

#Influenza A and RSV: New wastewater sampling sites, or sites with a substantial change in laboratory methods, are included in national and state or territorial median values beginning on August 1st of each year once there are at least 8 weeks of samples reported for that pathogen. Data must be reported by October 1st of that year to be included in national and state or territorial median values for that respiratory virus season. If data are reported after October 1st of that year, they will not be displayed until August 1st of the following year.


## Load dependencies
# load libraries
library(pacman)
p_load(tidyverse, lubridate, here, REDCapR, sf, zoo, ggformula, slider, stringr)


## Helper functions
#################################################
# Quantile approximation to mimic PySpark
#################################################
quantile_approx <- function(x, p) {
  x <- sort(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  if (p >= 1) return(x[length(x)])
  i_x <- p * (length(x) - 1) + 1
  x[floor(i_x)]
}


##########################################################
# Function for handling n1/n2 geometric mean to mimic CDC
##########################################################
# CDC calculates a geometric mean of samples for the same pcr_target from the same site/day/MLM
# This function applies for sars-cov-2 samples from the same sample ID but different gene targets
# Results below the LOD are replaced with the geometric mean of LOD/2 for the samples
# LOD sewage for both n1/n2 should be the same, so geometric mean of LOD/2 is just LOD/2
calc_n1n2_avg <- function(data) {
  samples_with_both <- data %>%
    # identify sample IDs that have BOTH n1 and n2
    # filter to sars n1/n2 gene targets
    filter(pcr_target == 'sars-cov-2', 
           pcr_gene_target %in% c('n1', 'n2')) %>%
    # group by sample info/date/lab/method
    group_by(sample_id) %>%
    summarise(n_gene_targets = n_distinct(pcr_gene_target), .groups = 'drop') %>%
    filter(n_gene_targets == 2) %>%
    pull(sample_id)
  
  # if no paired samples, just add the column and return
  if (length(samples_with_both) == 0) {
    return(data %>% mutate(pcr_gene_target_agg = as.character(pcr_gene_target)))
  }
  
  # process samples with BOTH n1 and n2 and calculate geometric mean of them
  sars_averaged <- data %>%
    # filter to samples with both n1 and n2 with same sample ID 
    filter(sample_id %in% samples_with_both) %>%
    group_by(sample_id) %>%
    # calculate intermediate values
    reframe(
      # average LOD across n1 and n2
      avg_lod = mean(lod_sewage, na.rm = TRUE),
      
      # count how many values are above LOD
      n_above_lod = sum(pcr_target_avg_conc >= lod_sewage, na.rm = TRUE),
      
      # calculate the geometric mean of the n1/n2 concentrations, else set to LOD/2
      pcr_target_avg_conc = if_else(
        n_above_lod > 0,
        10^(mean(log10(pcr_target_avg_conc), na.rm = TRUE)),
        mean(lod_sewage, na.rm = TRUE) / 2
      ),
      
      # update lod_sewage to averaged value
      lod_sewage = avg_lod,
      
      # create pcr_gene_target_agg variable to match CDC nwss
      pcr_gene_target_agg = 'n2 and n1'
    ) %>%
    
    ungroup() %>%
    
    # keep only one row per sample (remove the duplicate N1/N2 rows)
    distinct(sample_id, .keep_all = TRUE) %>%
    
    # remove intermediate calculation columns
    select(-avg_lod, -n_above_lod)
  
  # get samples with only individual results (keep all other samples as is)
  sars_individual <- data %>%
    filter(!sample_id %in% samples_with_both) %>%
    mutate(pcr_gene_target_agg = pcr_gene_target)
  
  # combine averaged and individual results
  final_result <- bind_rows(sars_averaged, sars_individual)
  
  return(final_result)
}

################################################################
# Functions for calculating days to nearest recalculation target
################################################################
# Baseline recalculates on April 1st and Oct 1st of every year for sars-cov-2
# Baseline recalculates on August 1st of every year for FLUAV/RSV
# These functions create a variable that is a count of the number of days to nearest recalculation date
# Use before baseline calculations

# Days to nearest Apr 1 / Oct 1 (SARS-CoV-2)
calc_days_to_target_sars <- function(date) {
  date    <- as.Date(date)
  y       <- year(date)
  april_1 <- as.Date(sprintf("%d-04-01", y))
  oct_1   <- as.Date(sprintf("%d-10-01", y))
  m <- month(date)
  if (m >= 4 && m <= 9) abs(as.integer(difftime(date, april_1, units = "days")))
  else                  abs(as.integer(difftime(date, oct_1,   units = "days")))
}

# Days to nearest Aug 1 (FLU/RSV)
calc_days_to_target_flu_rsv <- function(date) {
  date <- as.Date(date)
  y    <- year(date)
  aug1 <- as.Date(sprintf("%d-08-01", y))
  augN <- as.Date(sprintf("%d-08-01", y + 1))
  pmin(abs(as.integer(difftime(date, aug1, units = "days"))),
       abs(as.integer(difftime(date, augN, units = "days"))))
}



##################################################
## Function for handling outliers
##################################################
# Function to remove outliers with z-score > 4
# data - the pre-processed ww data containing the site_id_with_pcr_source_mlm
# and pcr_target_avg_conc_ln
# The output is a df with outliers removed

remove_outliers <- function(data) {
  ## calculate mean and standard deviation for each site_lab_mlm combo
  mean_std_ln <- data %>%
    group_by(site_id_with_pcr_source_mlm) %>%
    summarise(
      mean_avg_conc_ln = mean(pcr_target_avg_conc_ln, na.rm = TRUE),
      stddev_avg_conc_ln = sd(pcr_target_avg_conc_ln, na.rm = TRUE)) 
  
  clean_df <- left_join(data, mean_std_ln, by = "site_id_with_pcr_source_mlm")
  
  # set z-score threshold
  z_score_threshold <- 4
  
  ## calculate and filter the Z-score, handling NA stddev and division by zero
  clean_df <- clean_df %>%
    mutate(z_score = ifelse((is.na(stddev_avg_conc_ln) | stddev_avg_conc_ln == 0), 0, ### assign a z-score to ensure rows won't be removed
                            abs(pcr_target_avg_conc_ln - mean_avg_conc_ln) / stddev_avg_conc_ln))
  
  ## store outliers
  outliers <- clean_df %>%
    filter(z_score > z_score_threshold)
  
  ## filter out outliers based on the z-score threshold
  clean_df <- clean_df %>% 
    filter((z_score <= z_score_threshold) | is.na(z_score)) %>%
    select(-mean_avg_conc_ln, -stddev_avg_conc_ln, -z_score)
  return(clean_df)
}

##################################################
## Function for categorizing wval 
##################################################
# This function categorizes wastewater index values into wastewater viral activity levels by pathogen
# The current Wastewater Viral Activity Level for each state and territory is categorized into very low, low, moderate, high, or very high as follows:
# SARS-CoV-2:
# Very Low: Up to 2
# Low: Greater than 2 and up to 3.4	
# Moderate: Greater than 3.4 and up to 5.3	
# High: Greater than 5.3 and up to 7.8	
# Very High: Greater than 7.8

# Influenza A:
# Very Low: Up to 2.7
# Low: Greater than 2.7 and up to 6.2	
# Moderate: Greater than 6.2 and up to 11.2	
# High: Greater than 11.2 and up to 17.6	
# Very High: Greater than 17.6

# RSV:
# Very Low: Up to 2.5	
# Low: Greater than 2.5 and up to 5.2	
# Moderate: Greater than 5.2 and up to 8	
# High: Greater than 8 and up to 11	
# Very High: Greater than 11 

categorize_wval <- function(ww_index, target) {
  t <- tolower(target)
  out <- rep(NA_character_, length(ww_index))
  
  is_sars <- stringr::str_detect(t, "sars")
  is_flu  <- stringr::str_detect(t, "fluav")
  is_rsv  <- stringr::str_detect(t, "rsv")
  
  out[is_sars] <- as.character(cut(
    ww_index[is_sars],
    breaks = c(-Inf, 2.0, 3.4, 5.3, 7.8, Inf),
    labels = c("Very Low", "Low", "Moderate", "High", "Very High"),
    right  = TRUE
  ))
  
  out[is_flu] <- as.character(cut(
    ww_index[is_flu],
    breaks = c(-Inf, 2.7, 6.2, 11.2, 17.6, Inf),
    labels = c("Very Low", "Low", "Moderate", "High", "Very High"),
    right  = TRUE
  ))
  
  out[is_rsv] <- as.character(cut(
    ww_index[is_rsv],
    breaks = c(-Inf, 2.5, 5.2, 8, 11, Inf),
    labels = c("Very Low", "Low", "Moderate", "High", "Very High"),
    right  = TRUE
  ))
  
  out
}


## Function for generating WVAL for sars-cov-2, FLU A, RSV
wval_v2 <- function(data, target) { 
  
  #########################
  # Data Prep
  #########################
  # Function Inputs: 
  # - data: raw wastewater dataset
  # - pcr_target: pathogen (SARS-CoV-2, FLUAV, RSV)
  
  # Dataset Variables Needed 
  # - totalpop: sewershed population
  # - lod_sewage: limit of detection
  # - pcr_target_avg_conc: SARS-CoV-2 raw concentration in gene copies/L
  # - sample_site_id: site ID 
  # - funding: the source of the data submitted to 1CDP (NWSS, WWS, CDC)
  # - micro_lab_name: ID assigned to the testing lab
  # - sample_site_name: name of wastewater treatment facility sample collection location
  # - county: county containing wastewater treatment facility sample collection location
  # - pcr_gene_target: specific pathogen gene target being tested (i.e., n for SARS-CoV-2)
  # - concentration_method: method used to concentrate the sample prior to analysis of the concentrate
  # - extraction_method: method used for nucleic acid extraction from the sample
  # - major_lab_method: A number used to distinguish major lab methods at the reporting jurisdiction level
  
  # process n1/n2 sars-cov-2 data
  data <- calc_n1n2_avg(data)
  
  # prepare dataset and apply filters
  all_final <- data %>% 
    # filter to pcr target
    filter(pcr_target == target) %>% 
    
    
    # create dummy variable for below LOD yes/no
    mutate(
      pcr_target_below_lod = if_else(pcr_target_avg_conc < lod_sewage, 'yes', 'no'),
      # pcr target detection yes/no
      pcr_target_detect = if_else(pcr_target_avg_conc < lod_sewage, 'no', 'yes'),
      # if target avg concentration is below LOD, set to 1/2*LOD
      # otherwise leave as is
      pcr_target_avg_conc_lin = if_else(pcr_target_below_lod == 'yes',
                                        lod_sewage/2, 
                                        pcr_target_avg_conc),
      # take the log10 of the target avg concentration
      pcr_target_avg_conc_log10 = if_else(pcr_target_below_lod == 'yes',
                                          log10(lod_sewage/2),
                                          log10(pcr_target_avg_conc)),
      # calculate natural log of linear avg conc
      pcr_target_avg_conc_ln = log(pcr_target_avg_conc_lin),
      # calculate year
      year = year(sample_collect_date)
    )
  
  all_final1 <- all_final %>%
    # rename variables to match dcipher/1CDP code
    rename(
      site_id = sample_site_id, 
      source = funding,
      lab_id = micro_lab_name,
      wwtp_name = sample_site_name,
      population_served = totalpop,
      # county dashboard names
      county_names = county_main
    ) %>%
    
    # drop NAs 
    drop_na(site_id) %>%
    drop_na(sample_collect_date) %>%
    drop_na(population_served) %>%
    
    # select variables of interest
    select(
      sample_collect_date, 
      site_id, 
      source, 
      lab_id, 
      wwtp_name,
      population_served,
      pcr_target,
      pcr_gene_target,
      concentration_method, 
      extraction_method, 
      pcr_target_detect, 
      pcr_target_avg_conc, 
      pcr_target_avg_conc_log10,
      pcr_target_avg_conc_lin,
      pcr_target_avg_conc_ln,
      lod_sewage,
      pcr_target_below_lod,
      county_names,
      major_lab_method
    ) %>%
    # add new column that combines site_id and pcr_target
    mutate(
      site_id_with_pcr = paste(site_id, pcr_target, sep = "_"),
      # add new column that combines site_id, pcr_target, source
      site_id_with_pcr_source = paste(site_id, pcr_target, source, sep = "_"),
      # add new column that combines site_id, pcr_target, source, and major lab method
      site_id_with_pcr_source_mlm = paste(site_id_with_pcr, source, major_lab_method, sep = "_")) %>%
    # group by site, pcr target, and lab method to calculate days since first sample
    group_by(site_id_with_pcr_source_mlm) %>%
    # calculate days since first sample and first sample collection date
    mutate(
      days_since_first_sample = as.integer(difftime(sample_collect_date, min(sample_collect_date), units = "days")),
      first_collection_date = min(sample_collect_date)) %>%
    ungroup()
  
  # remove outliers function call
  all_final2 <- remove_outliers(data=all_final1)
  
  
  #####################################################
  # Calculating days to nearest recalculation target
  #####################################################
  
  # create dataframe with days to nearest target variables based on pathogen type
  if(target == 'sars-cov-2'){
    all_final3 <- all_final2 %>% 
      # fill in missing sample dates
      arrange(sample_collect_date, .by_group = TRUE) %>%
      complete(sample_collect_date = seq.Date(min(sample_collect_date, na.rm = TRUE),
                                              max(sample_collect_date, na.rm = TRUE),
                                              by = 'day')) %>%
      # days to nearest target for sars-cov-2 based on april and october recalculation dates
      mutate(half_year = if_else(month(sample_collect_date) >= 4 & month(sample_collect_date) <= 9, 1, 2),
             year = year(sample_collect_date),
             days_to_nearest_target = sapply(sample_collect_date, calc_days_to_target_sars)) %>%
      group_by(site_id_with_pcr_source_mlm, half_year, year) %>%
      # calculate minimum days to target for sars-cov-2
      mutate(min_days_to_target = min(days_to_nearest_target, na.rm = TRUE)) %>% 
      ungroup() %>%
      # drop helper variable
      select(-half_year)
    
    # influenza A and RSV recalculate on august 1
  } else {
    # calculate days to nearest target for flu/rsv
    all_final3 <- all_final2 %>%
      # fill in missing sample dates
      arrange(sample_collect_date, .by_group = TRUE) %>%
      complete(sample_collect_date = seq.Date(min(sample_collect_date, na.rm = TRUE),
                                              max(sample_collect_date, na.rm = TRUE),
                                              by = 'day')) %>%
      # add flu season year: august 1st recalculation
      mutate(flu_season_year = if_else(month(sample_collect_date) >= 8,
                                       year(sample_collect_date),
                                       year(sample_collect_date) - 1L),
             # days to nearest august 1st recalculation target date
             days_to_nearest_target = sapply(sample_collect_date, calc_days_to_target_flu_rsv)) %>%
      group_by(site_id_with_pcr_source_mlm, flu_season_year) %>%                       
      mutate(min_days_to_target = min(days_to_nearest_target, na.rm = TRUE)) %>%
      ungroup() %>%
      # drop helper variable
      select(-flu_season_year)
  }
  
  ##################################################
  ## Baseline Calculation
  ##################################################
  # Calculate baseline (10th percentile) and standard deviation for last 24 months 
  # First 6 months or since the nearest to Apr 1/Oct 1 for sars-cov-2
  # First 12 months or since nearest Aug 1st for FLU A/RSV
  # Input dataset with outliers removed
  # create half year - specific to SARS
  # Output dataset with baseline, sd, and half year variables
  
  if(target == "sars-cov-2") {
    all_final4 <- all_final3 %>%
      group_by(site_id_with_pcr_source_mlm) %>%
      mutate(baseline_flag = (days_since_first_sample <= 182) | (days_to_nearest_target == min_days_to_target), 
             temp_baseline = if_else(
               baseline_flag,
               slide_index_dbl(pcr_target_avg_conc_ln, sample_collect_date,
                               .before = 730, .after = 0, .complete = FALSE,
                               .f = ~ quantile_approx(.x, 0.10)),
               NA_real_
             ),
             baseline_avg_conc_ln = na.locf(temp_baseline, na.rm= FALSE),
             temp_std = if_else(
               baseline_flag,
               slide_index_dbl(pcr_target_avg_conc_ln, sample_collect_date,
                               .before = 730, .after = 0, .complete = FALSE,
                               .f = ~ sd(.x, na.rm = TRUE)),
               NA_real_
             ),
             std_avg_conc_ln = na.locf(temp_std, na.rm = FALSE)
      ) %>%
      select(-c(temp_std, temp_baseline, baseline_flag)) %>%
      drop_na(pcr_target) %>%
      ungroup()
  } 
  else {
    all_final4 <- all_final3 %>%
      group_by(site_id_with_pcr_source_mlm) %>%
      mutate(baseline_flag = (days_since_first_sample <= 365) | (days_to_nearest_target == min_days_to_target),
             temp_baseline = if_else(
               baseline_flag,
               slide_index_dbl(pcr_target_avg_conc_ln, sample_collect_date,
                               .before = 730, .after = 0, .complete = FALSE,
                               .f = ~ quantile_approx(.x, 0.10)),
               NA_real_
             ),
             baseline_avg_conc_ln = na.locf(temp_baseline, na.rm = FALSE),
             temp_std = if_else(
               baseline_flag,
               slide_index_dbl(pcr_target_avg_conc_ln, sample_collect_date,
                               .before = 730, .after = 0, .complete = FALSE,
                               .f = ~ sd(.x, na.rm = TRUE)),
               NA_real_
             ),
             std_avg_conc_ln = na.locf(temp_std, na.rm = FALSE)
      ) %>%
      select(-temp_baseline, -temp_std, -baseline_flag) %>%
      drop_na(pcr_target) %>%
      ungroup()
  }
  
  # saving last calculated baseline and standard deviation
  last_vals <- all_final4 %>%
    filter(!is.na(std_avg_conc_ln) & !is.na(baseline_avg_conc_ln)) %>%
    group_by(site_id_with_pcr_source_mlm) %>%
    slice_tail(n = 1) %>%
    summarise(
      last_baseline_ln = baseline_avg_conc_ln,
      last_std_ln = std_avg_conc_ln
    )
  
  #################################################
  # Calculate wastewater index values
  #################################################
  
  # create result dataframe and join last calculated baseline/standard deviation
  result <- all_final4 %>%
    left_join(last_vals, by = "site_id_with_pcr_source_mlm") %>%
    group_by(site_id_with_pcr_source_mlm) %>%
    # create variable for the number of days that site has been active
    mutate(
      days_site_online = max(days_since_first_sample, na.rm = TRUE)) %>% 
    ungroup() %>% 
    # create wastewater index variables, excluding sites that have been online for less than 56 days for sars-cov-2   
    mutate(
      ww_index_ln = case_when(
        days_site_online < 56 ~ NA_real_, TRUE ~ (pcr_target_avg_conc_ln - last_baseline_ln) / last_std_ln),
      ww_index_ln_lin = exp(ww_index_ln),
      # create day of week and week end variables
      day_of_week = wday(sample_collect_date),
      week_end = sample_collect_date + days(7-day_of_week)
    ) %>%
    group_by(site_id_with_pcr_source_mlm, week_end) %>%
    arrange(sample_collect_date, .by_group = TRUE) %>%
    ungroup() %>%
    distinct()
  
  # create intermediate dataset
  intermediate_sars <<- result
  
  
  ##################################################
  ## MLM Overlap Handling
  ##################################################
  # if a site_id_with_pcr has more than one MLM reporting in a week,
  # keep data only if ww_index_ln_lin is non-null and it is from a more recent MLM
  result <- result %>%
    group_by(site_id_with_pcr_source_mlm) %>%
    mutate(min_sample_date_for_mlm = min(sample_collect_date, na.rm = TRUE)) %>%
    ungroup() %>%
    # determines latest min_sample_date_for_mlm across all site_id_with_pcr_source_mlm
    # that share the same site_id_with_pcr_source
    # using rank descending in case of a situation where these min dates happen to be the same (in which case they will both be 1)
    group_by(site_id_with_pcr, week_end) %>%
    arrange(desc(min_sample_date_for_mlm), desc(major_lab_method)) %>%
    mutate(latest_min_date_rank = rank(-as.numeric(min_sample_date_for_mlm), ties.method = "min")) %>%
    ungroup() %>%
    
    # filter to retain only rows from the latest site_id_with_pcr_source_mlm with non-null ww_index_ln_lin
    filter(latest_min_date_rank == 1, !is.na(ww_index_ln_lin))
  
  
  #################################################
  # Categorize WVAL
  #################################################
  # at the result level
  unaggregated_result <- result %>%
    mutate(wval_result = categorize_wval(ww_index_ln_lin, target)) %>% distinct()
  
  unaggregated_result <<- unaggregated_result
  
  
  
  #################################################
  # Calculate the site weekly average for WVAL 
  #################################################
  # function for aggregating by week, site, and target
  # creates categories for categorizing pathogen ww index values
  site_weekly_average <- function(data) {
    site_ag <- data %>%
      # remove NA ww_index_ln_lin, arrange by sample_collect_date and group by site and week
      filter(!is.na(ww_index_ln_lin)) %>%
      arrange(desc(sample_collect_date)) %>%
      group_by(site_id_with_pcr, week_end) %>%
      
      reframe(
        site_id = first(site_id),
        mean_pcr_target_avg_conc_lin = mean(pcr_target_avg_conc_lin, na.rm = TRUE),
        all_pcr_target_avg_conc_lin = paste(pcr_target_avg_conc_lin, collapse = ", "),
        wval_site_wk = mean(ww_index_ln_lin, na.rm = TRUE),
        all_sample_collect_dates = paste(sample_collect_date, collapse = ", "),
        days_since_first_sample = last(days_since_first_sample),
        baseline_avg_conc_ln  = last(baseline_avg_conc_ln),
        last_baseline_ln = last(last_baseline_ln),
        last_std_ln = last(last_std_ln),
        std_avg_conc_ln = last(std_avg_conc_ln),
        pcr_target_below_lod = first(pcr_target_below_lod),
        pcr_target = first(pcr_target),
        all_ww_index_ln_lin = paste(ww_index_ln_lin, collapse = ", "),
        source = source,
        county_names = county_names,
        .groups = "drop"
      ) %>% 
      # categorize wval
      mutate(wval_site_wk_level = categorize_wval(wval_site_wk, target)
      )%>% distinct()
  }
  
  
  ############################
  ## Return aggregated results 
  ############################            
  
  site_ag <- site_weekly_average(unaggregated_result)
  
  return(site_ag)
  
}




### NEW FUNCTION FOR COUNTY LEVEL
# Calculates the WVAL for the county level aggregation
# takes the median of the site-week aggregated WVALs in a county
# data - data must contain site-week aggregated WVAL values
# Outputs a df containing county-week aggregated WVAL values  
wval_county_wk_agg <- function(data) {
  # aggregate by county 
  county_ag <- data %>% 
    # group by county and week and target
    group_by(county_names, week_end, pcr_target) %>%
    # aggregate
    reframe(week_end = first(week_end),
            wval_county_wk = median(wval_site_wk)
    ) %>%
    
    unique() %>%
    arrange(county_names, week_end) %>% ungroup()
  
  
  # calculate wval categories
  county_ag <- county_ag %>%
    mutate(wval_county_wk_level = categorize_wval(wval_county_wk, pcr_target))
  
  return(county_ag)
}

### NEW FUNCTION FOR STATE LEVEL
# Calculates the WVAL for the state level aggregation
# takes the median of the county-week aggregated WVALs in the state
# data - data must contain county-week aggregated WVAL values
# Outputs a df containing county-week aggregated WVAL values  
wval_state_wk_agg <- function(data) {
  # aggregate by state
  state_ag <- data %>% 
    # group by county and week and target
    group_by(week_end, pcr_target) %>%
    # aggregate
    reframe(week_end = first(week_end),
            wval_state_wk = median(wval_site_wk)
    ) %>%
    
    unique() %>%
    arrange(week_end) %>% ungroup()
  
  
  # calculate wval categories
  state_ag <- state_ag %>%
    mutate(wval_state_wk_level = categorize_wval(wval_state_wk, pcr_target))
  
  return(state_ag)
}
