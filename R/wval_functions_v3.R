## Overview: 
#- This is the collaborative NWSS CoE R translation of CDC's 2026 python code to calculate wastewater viral activity level (WVAL) for SARS-CoV-2, Influenza, and RSV with the NEW updated methodology published in August 2025. 

#- WVAL categorizes viral concentration into minimal, low, moderate, high, and very high to help indicate the risk of infection. 


## From CDC's methodology: 
  
### Calculating the Wastewater Viral Activity Level
  
### Major Changes (August 2026): 
#- Outlier handling switched from a Z-score drop to a 99th percentile Winsorization (capping extreme values instead of dropping values with Z-scores > 4).
#- Dual lab-based and site-specific baselines: 
# - Uses lab-based for first 6-12 months until site has established enough history, then switches to site-specific.
# - Lab-based baselines are grouped by pathogen and lab instead of site.
#- Corrected calendar target windows:
#  - For Flu/RSV: Updates are only triggered if closest sample to August 1st falls between August and September
#- For SARS-CoV-2: Splits tracking into half-year 1: April window, and half-year 2: October window
#- Handles sampling gaps of over 182 days by shifting major lab method ID
#- Discrete quantile calculation used by Polars instead of PySpark
#- Updated WVAL level cut off points

### Baseline Calculation:
#- For each combination of site, data submitter, PCR target, lab methods, and normalization method, a baseline is established. The “baseline” is the 10th percentile of the log-transformed concentration data within a specific time frame. Details on the baseline calculation by pathogen are below:
  
#- There are now 2 baselines that will be created: Lab based and site (group) based. 
#** This only impacts new sites. The lab based method will be used for new sites until they reach a certain threshold 
#(365 days for flua and rsv and 182 days for sars) after which they will use the site/ group based method. 

#- SARS-CoV-2
#- For site and method combinations (as listed above) with over 6 months of data, baselines are re-calculated every six calendar months (April 1st and October 1st) using the past 24 months of data.
#- For sites and method combinations with less than 6 months of data, baselines are computed every time there is a new sample until reaching six months, after which they remain unchanged until the next April 1st or October 1st, at which time baselines are re-calculated.

#- Influenza A and RSV
#- For site and method combinations (as listed above) with over 12 months of data, baselines are re-calculated every August 1st using all available data in the previous 24 months.
#- For sites and method combinations with less than 12 months of data, baselines are computed weekly until reaching twelve months, after which they remain unchanged until the next August 1st, at which time baselines are re-calculated.

#- The standard deviation for each site and method combination is calculated using the same time frame as the baseline.

### Wastewater Viral Activity Level Calculation:
#- The number of standard deviations that each log-transformed concentration value deviates from the baseline (positive if above, negative if below) is calculated.
#- This value (x) is then converted back to a linear scale (by calculating e^x) to form the Wastewater Viral Activity Level for the site and method combination.
#- The Wastewater Viral Activity Levels from a site are averaged by week for all figures.

#- This reflects the updated cut points: 
# - Wastewater Viral Activity Level Categories for SARS-CoV-2 (https://www.cdc.gov/nwss/about-data.html#data-method, Feb 2025): 
# - The current Wastewater Viral Activity Level for each state and territory is categorized into very low, low, moderate, high, or very high as follows:
#                                                                 
#SARS-CoV-2:
#- Very Low: Up to 2.6
#- Low: Greater than 2.6 and up to 4.9	
#- Moderate: Greater than 4.9 and up to 7.9	
#- High: Greater than 7.9 and up to 11.6	
#- Very High: Greater than 11.6
                                                               
#Influenza A:
#- Very Low: Up to 2.4
#- Low: Greater than 2.4 and up to 5.5	
#- Moderate: Greater than 5.5 and up to 10.2	
#- High: Greater than 10.2 and up to 15.6	
#- Very High: Greater than 15.6
                                                               
#RSV:
#- Very Low: Up to 1.7	
#- Low: Greater than 1.7 and up to 3.4	
#- Moderate: Greater than 3.4 and up to 5.4	
#- High: Greater than 5.4 and up to 8.1	
#- Very High: Greater than 8.1 
                                                               
                                                               
### Aggregation for National, Regional, and State Levels:
#- We calculate the median Wastewater Viral Activity Levels among sites at national, regional, and state levels, excluding data from site/method combinations with less than 8 weeks of data for SARS-CoV-2, Influenza A, and RSV.
                                                               
### Data Inclusion Criteria – SARS-COV-2, Influenza A, RSV Wastewater Viral Activity Level
#- SARS-CoV-2: New wastewater sampling sites, or sites with a substantial change in laboratory methods are included in national, regional, state, or territorial median values once there are at least 8 weeks of samples reported for that location.
#- Influenza A and RSV: New wastewater sampling sites, or sites with a substantial change in laboratory methods, are included in national and state or territorial median values beginning on August 1st of each year once there are at least 8 weeks of samples reported for that pathogen. Data must be reported by October 1st of that year to be included in national and state or territorial median values for that respiratory virus season. If data are reported after October 1st of that year, they will not be displayed until August 1st of the following year.
                                                               
## Data Needs
                                                               
### Input Data
#- 1CDP downloaded MLM analytics file OR 
#- Raw data from your jurisdiction 
                                                               
## Dataset Variables Needed (using 1CDP-ready column names)
#- wwtp_jurisdiction
#- reporting_jurisdiction
#- pcr_type
#- source: the source of the data submitted to 1CDP (NWSS, WWS, CDC)
#- lab_id: ID assigned to the testing lab
#- site_id: site ID 
#- wwtp_name: name of wastewater treatment facility sample collection location
#- population_served: sewershed population
#- pcr_target: pathogen (sars-cov-2, fluav, rsv)
#- pcr_gene_target: specific pathogen gene target being tested (i.e., n for SARS-CoV-2)
#- pcr_target_avg_conc: SARS-CoV-2 raw concentration in gene copies/L
#- lod_sewage: limit of detection
#- county_names: county containing wastewater treatment facility sample collection location
#- zipcode
#- concentration_method: method used to concentrate the sample prior to analysis of the concentrate
#- extraction_method: method used for nucleic acid extraction from the sample
#- major_lab_method: A number used to distinguish major lab methods at the reporting jurisdiction level
#- sample_location: Sample collection location in the wastewater system, whether at a wastewater treatment plant or upstream in the wastewater system.
#- sample_location_specify: If 'sample_location' is "upstream", specify the collection location in the wastewater system; an arbitrary name may be used if you do not wish to disclose the real name. 
#- sample_matrix: Wastewater matrix from which the sample was collected.
                                                               
#- NOTE: make sure the pcr_target for all sars-cov-2 data is 'sars-cov-2'
                                                               
## Dependencies
#- Libraries: pacman, tidyverse, lubridate, here, REDCapR, sf, zoo, ggformula, slider, stringr, ggplot2
                                                               
  
# load libraries
library(pacman)
p_load(tidyverse, lubridate, here, REDCapR, sf, zoo, ggformula, slider, stringr, ggplot2, data.table)

## Sars-cov-2 Adjustments

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
         avg_lod / 2
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


## Helper functions
  ########################################################
  # Function to handle sampling timeline gaps of >182 days
  # by shifting the major lab method ID 
  ########################################################
  adjust_mlm_for_gaps <- function(data) {
  keys <- c("site_id", "source", "major_lab_method")
  data <- data %>%
    arrange(across(all_of(c(keys, "sample_collect_date")))) %>%
    group_by(across(all_of(keys))) %>%
    mutate(
      days_diff = as.integer(sample_collect_date - lag(sample_collect_date)),
      is_gap = if_else(is.na(days_diff), FALSE, days_diff > 182),
      gap_id = cumsum(as.integer(is_gap))
    ) %>%
    ungroup() %>%
    mutate(
      major_lab_method = if_else(gap_id > 0 & !is.na(as.numeric(major_lab_method)),
                                 formatC(round(as.numeric(major_lab_method) + (gap_id * 0.1), 1), format = "f", digits = 1),
                                 as.character(major_lab_method))

    ) %>%
    select(-gap_id, -days_diff, -is_gap)
  return(data)
  }

  #################################################
  # Quantile approximation to mimic Polars
  # Type 6 in base R quantile uses linear interpolation
  # method similar to Polars
  #################################################
  polars_quantile <- function(x, probs, na.rm = TRUE) {
    x <- x[!is.na(x) & !is.nan(x)]
    n <- length(x)
    if (n == 0) return(NA_real_)
    x <- sort(x)
    pos <- probs * (n - 1)
    idx <- floor(pos) + 1L   # 0-indexed floor, converted to R's 1-indexing
    x[idx]
  }

  #################################################
  # Quantile approximation to mimic Polars for winsorsation
  # method similar to Polars
  #################################################
  quantile_nearest <- function(x, q) {
    x <- x[!is.na(x) & !is.nan(x)]
    n <- length(x)
    if (n == 0) return(NA_real_)
    x <- sort(x)
    pos <- q * (n - 1)
    idx <- round(pos) + 1L        # round-half-to-even, matches R's round() and Polars' tie behavior
    idx <- max(1L, min(n, idx))
    x[idx]
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
    date <- as.Date(date)
    y <- year(date)
    m <- month(date)
    target_date <- if (m >= 4 && m <= 9) {
    as.Date(sprintf("%d-04-01", y))
    } else {
    as.Date(sprintf("%d-10-01", y))
  }
  as.integer(difftime(date, target_date, units = "days"))
  }

  # Days to nearest Aug 1 (FLU/RSV)
  calc_days_to_target_flu_rsv <- function(date) {
    date <- as.Date(date)
    y <- year(date)
    aug1 <- as.Date(sprintf("%d-08-01", y))
    as.integer(difftime(date, aug1, units = "days"))
  }

  # min() over only the "candidate" rows; NA if there are no candidates in the group
  safe_min_candidate <- function(days, is_candidate) {
    vals <- days[is_candidate]
    if (length(vals) == 0 || all(is.na(vals))) return(NA_real_)
    min(vals, na.rm = TRUE)
  }


  ##################################################
  ## Function for handling outliers
  ##################################################
  # Uses a 99th percentile Winsorization method
  # data - the pre-processed ww data containing the GROUP_KEY
  # and pcr_target_avg_conc_ln
  # The output is a df with outliers above 99th percentile clipped and a df of outliers

  get_outliers <- function(data) {
  # calculate 99th percentile cut off 
    percentile_cutoffs <- data %>%
      group_by(GROUP_KEY) %>%
      summarise(
      winsor_perc_99_cutoff = quantile_nearest(pcr_target_avg_conc_ln, 0.99),
      .groups = "drop") 
  
  clean_df <- left_join(data, percentile_cutoffs, by = "GROUP_KEY")
  
  # winsorize only if a site has been sampling over 56 days and is over the 99th percentile cutoff 
  should_winsorize <- (clean_df$days_since_first_sample > 56) & (clean_df$pcr_target_avg_conc_ln >   clean_df$winsor_perc_99_cutoff)
  
  # NA-guard: any NA in the comparison inputs makes should_winsorize NA in R, 
  # which if_else() would otherwise propagate as NA instead of treating as FALSE. 
  # Matches Polars' `.fill_null(False)` on the equivalent boolean. 
  should_winsorize[is.na(should_winsorize)] <- FALSE 
  
  # calculate winsorized values and identify outliers
  clean_df <- clean_df %>%
    mutate(
      pcr_target_avg_conc_ln_winsor = if_else(should_winsorize, winsor_perc_99_cutoff, pcr_target_avg_conc_ln),
      winsor_outlier = if_else(should_winsorize, "yes", "no"))
  
  # store outliers but DO NOT DROP 
  outliers_identified <- clean_df %>% filter(winsor_outlier == "yes")
  winsorized <- clean_df
  
  # return both dfs
  return(list(outliers_identified = outliers_identified, winsorized = winsorized))
  }
  ##################################################
  ## Function for categorizing wval Updated July 2026
  ##################################################
  # This function categorizes wastewater index values into wastewater viral activity levels by pathogen
  # The current Wastewater Viral Activity Level for each state and territory is categorized into very low, low, moderate, high, or very high as follows:
  # SARS-CoV-2:
    # Very Low: Up to 2.6
    # Low: Greater than 2.6 and up to 4.9	
    # Moderate: Greater than 4.9 and up to 7.9	
    # High: Greater than 7.9 and up to 11.6	
    # Very High: Greater than 11.6

  # Influenza A:
    # Very Low: Up to 2.4
    # Low: Greater than 2.4 and up to 5.5	
    # Moderate: Greater than 5.5 and up to 10.2	
    # High: Greater than 10.2 and up to 15.6	
    # Very High: Greater than 15.6

  # RSV:
    # Very Low: Up to 1.7	
    # Low: Greater than 1.7 and up to 3.4	
    # Moderate: Greater than 3.4 and up to 5.4	
    # High: Greater than 5.4 and up to 8.1	
    # Very High: Greater than 8.1 

  categorize_wval <- function(ww_index, target) {
    out <- rep(NA_character_, length(ww_index))
  
  is_sars <- stringr::str_detect(target, "sars")
  is_flu  <- stringr::str_detect(target, "FLUAV")
  is_rsv  <- stringr::str_detect(target, "RSV")
  
  out[is_sars] <- as.character(cut(
    ww_index[is_sars],
    breaks = c(-Inf, 2.6, 4.9, 7.9, 11.6, Inf),
    labels = c("Very Low", "Low", "Moderate", "High", "Very High"),
    right  = TRUE
  ))
  
  out[is_flu] <- as.character(cut(
    ww_index[is_flu],
    breaks = c(-Inf, 2.4, 5.5, 10.2, 15.6, Inf),
    labels = c("Very Low", "Low", "Moderate", "High", "Very High"),
    right  = TRUE
  ))
  
  out[is_rsv] <- as.character(cut(
    ww_index[is_rsv],
    breaks = c(-Inf, 1.7, 3.4, 5.4, 8.1, Inf),
    labels = c("Very Low", "Low", "Moderate", "High", "Very High"),
    right  = TRUE
  ))
  
  out
  }

  ## Function for generating WVAL for sars-cov-2, FLU A, RSV
  wval_v3 <- function(data, target) { 
  
  # Prep the data
  data <- data %>%
    # ***Optional code (not necessary if your columns already match 1CDP data)
    # rename variables to match dcipher/1CDP code
    rename(
      site_id = sample_site_id,
      source = funding,
      lab_id = micro_lab_name,
      wwtp_name = sample_site_name,
      population_served = totalpop,
      county_names = county_main
    ) %>%
    # Convert sample_collect_date to date  
    mutate(sample_collect_date = as.Date(lubridate::parse_date_time(sample_collect_date, orders = c("mdY", "Ymd")))) %>%
    
    # Designate reporting jurisdiction for lab based key
    mutate(
          reporting_jurisdiction = "WA") %>%
    
    # Filter for the respiratory virus targets that are validated for WVAL
    filter(pcr_target %in% c("FLUAV", "RSV", "sars-cov-2")) %>%
    filter(!is.na(pcr_target_avg_conc))
  
  
  KEEP_COLS <- c(
    "sample_collect_date", "site_id", "sample_id", "sample_location", "source", "lab_id", "pcr_target_avg_conc", "pcr_target_avg_conc_lin", "reporting_jurisdiction", "wwtp_jurisdiction", "wwtp_name", "population_served", "pcr_target", "pcr_type","concentration_method", "extraction_method", "sample_matrix", "major_lab_method", "pcr_target_detect", "pcr_gene_target_agg", "lod_sewage", "county_names", "zipcode","pcr_target_below_lod")
  
  
  # Call the functions
  data <- data %>% 
    # process n1/n2 sars-cov-2 data
    # *** Only run this if it applies to your jurisdiction
    calc_n1n2_avg()
  
  # prepare dataset and apply filters
  all_final <- data %>% 
    # filter to pcr target and populations over 3000
    filter(pcr_target == target) %>% 
    select(any_of(KEEP_COLS)) %>% 
    
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
                                          log10(pcr_target_avg_conc_lin)),
      # calculate natural log of linear avg conc
      pcr_target_avg_conc_ln = log(pcr_target_avg_conc_lin),
      
      # calculate date variables
      year = year(sample_collect_date), 
      day_of_week = wday(sample_collect_date), 
      week_end = sample_collect_date + days(7 - day_of_week), 
      
      # create PCR Type2 variable
      pcr_type2 = case_when(pcr_type == "qpcr" ~ "qpcr", pcr_type %in% c("ddpcr", "dPCR", "dpcr", "qiagen dpcr", "life technologies dpcr", "dpcr,ddpcr") ~ "dpcr/ddpcr", TRUE ~ NA_character_)) %>% 
    
    # filter missing site IDs 
    filter(!is.na(site_id))
  
  
  # Adjust major lab method when there is a gap of more than 182 days within the same site/source/lab method series
  all_final1 <- adjust_mlm_for_gaps(all_final) %>%
    
    # ***Uncomment these next lines if they are applicable to your jurisdiction
    # filter(institution_type == "not institution specific",
    # dashboard_ignore == "no",
    # mlm_adjusted == "No",
    # pcr_target_detect !="") %>%
    
    # add new column that combines site_id and pcr_target
    mutate(
      site_id_with_pcr = paste(site_id, pcr_target, sep = "_"),
      # add new column that combines site_id, pcr_target, source
      site_id_with_pcr_source = paste(site_id, pcr_target, source, sep = "_"),
      
      # dynamic site group key assignment to mimic NWSS code, use group_key instead of site_id_with_pcr_source_mlm from here on out
      # add new column that combines site_id, pcr_target, source, and major lab method = Group Key
      GROUP_KEY = paste(site_id_with_pcr, source, major_lab_method, sep = "_"), 
      # add Lab key
      LAB_BASED_KEY = paste(reporting_jurisdiction, lab_id, pcr_target, source, pcr_type2, major_lab_method, sep = "_")) %>%
    
    #Group by GROUP_KEY 
    # *** Some jurisdictions may need this group_by line, if you are seeing a lot of discordant results try to uncomment this line
    group_by(GROUP_KEY) %>% 
    
    # calculate days since first sample and first sample collection date
    mutate(
      days_since_first_sample = as.integer(difftime(sample_collect_date, min(sample_collect_date), units = "days")),
      first_collection_date = min(sample_collect_date)) %>% 
    ungroup() 
  
  
  # get outliers function call
  outlier_result <- get_outliers(all_final1)
  
  all_final2 <- outlier_result$winsorized
  outliers_identified <- outlier_result$outliers_identified
  
  
  #####################################################
  # Calculating days to nearest recalculation target
  #####################################################
  
  # # create dataframe with days to nearest target variables based on pathogen type
  if (target == 'sars-cov-2') {
    all_final3 <- all_final2 %>%
      arrange(sample_collect_date) %>%
      mutate(
        half_year = if_else(month(sample_collect_date) >= 4 & month(sample_collect_date) <= 9, 1, 2),
        year = year(sample_collect_date),
        days_to_nearest_target = sapply(sample_collect_date, calc_days_to_target_sars),
        apr_may_window = (half_year == 1) & (month(sample_collect_date) %in% 4:5),
        oct_nov_window = (half_year == 2) & (month(sample_collect_date) %in% 10:11),
        target_candidate = (days_to_nearest_target >= 0) & (apr_may_window | oct_nov_window)
      ) %>%
      group_by(GROUP_KEY, half_year, year) %>%
      mutate(site_min_days_to_target = safe_min_candidate(days_to_nearest_target, target_candidate)) %>%
      ungroup() %>%
      group_by(LAB_BASED_KEY, half_year, year) %>%
      mutate(lab_min_days_to_target = safe_min_candidate(days_to_nearest_target, target_candidate)) %>%
      ungroup() %>%
      select(-half_year, -apr_may_window, -oct_nov_window)  
    
  } else {
    all_final3 <- all_final2 %>%
      arrange(sample_collect_date) %>%
      mutate(
        year_val = year(sample_collect_date),
        days_to_nearest_target = sapply(sample_collect_date, calc_days_to_target_flu_rsv),
        target_candidate = (month(sample_collect_date) %in% 8:9) & (days_to_nearest_target >= 0)
      ) %>%
      group_by(GROUP_KEY, year_val) %>%
      mutate(site_min_days_to_target = safe_min_candidate(days_to_nearest_target, target_candidate)) %>%
      ungroup() %>%
      group_by(LAB_BASED_KEY, year_val) %>%
      mutate(lab_min_days_to_target = safe_min_candidate(days_to_nearest_target, target_candidate)) %>%
      ungroup()
  }
  ##################################################
  ## Baseline Calculation
  ## (mirrors Python: create_site_or_lab_method -> apply_baseline_updates -> finish_baseline)
  ##################################################
  
  initial_days_val <- if (target == "sars-cov-2") 182 else 365
  
  ## ---- PHASE 1 (create_site_or_lab_method): unconditional rolling stats, every row ----
  
  # site-specific rolling p10 / std over GROUP_KEY - computed for EVERY row,
  # NOT gated by update_cond (that gating happens later, in Phase 2)
  all_final4 <- all_final3 %>%
    group_by(GROUP_KEY) %>%
    arrange(sample_collect_date, .by_group = TRUE) %>%
    mutate(
      site_roll_p10_24mo = slide_index_dbl(
        pcr_target_avg_conc_ln_winsor, sample_collect_date,
        .before = 730, .after = 0, .complete = FALSE,
        .f = ~ polars_quantile(.x, 0.10)
      ),
      site_roll_std_24mo = slide_index_dbl(
        pcr_target_avg_conc_ln_winsor, sample_collect_date,
        .before = 730, .after = 0, .complete = FALSE,
        .f = ~ sd(.x, na.rm = TRUE)
      )
    ) %>%
    ungroup()
  
  # collapse same-date duplicate rows within GROUP_KEY to one shared value
  # (matches Python's `.last().over([GROUP_KEY, "sample_collect_date"])`)
  all_final4 <- all_final4 %>%
    group_by(GROUP_KEY, sample_collect_date) %>%
    mutate(
      site_roll_p10_24mo = dplyr::last(site_roll_p10_24mo),
      site_roll_std_24mo = dplyr::last(site_roll_std_24mo)
    ) %>%
    ungroup()
  
  # lab-based rolling p10 / std over LAB_BASED_KEY - same pattern, independent sort
  all_final4 <- all_final4 %>%
    group_by(LAB_BASED_KEY) %>%
    mutate(
      first_lab_collection_date = min(sample_collect_date),
      days_since_first_lab_sample = as.integer(sample_collect_date - first_lab_collection_date)
    ) %>%
    arrange(sample_collect_date, .by_group = TRUE) %>%
    mutate(
      lab_roll_p10_24mo = slide_index_dbl(
        pcr_target_avg_conc_ln_winsor, sample_collect_date,
        .before = 730, .after = 0, .complete = FALSE,
        .f = ~ polars_quantile(.x, 0.10)
      ),
      lab_roll_std_24mo = slide_index_dbl(
        pcr_target_avg_conc_ln_winsor, sample_collect_date,
        .before = 730, .after = 0, .complete = FALSE,
        .f = ~ sd(.x, na.rm = TRUE)
      )
    ) %>%
    ungroup()
  
  all_final4 <- all_final4 %>%
    group_by(LAB_BASED_KEY, sample_collect_date) %>%
    mutate(
      lab_roll_p10_24mo = dplyr::last(lab_roll_p10_24mo),
      lab_roll_std_24mo = dplyr::last(lab_roll_std_24mo)
    ) %>%
    ungroup()
  
  # use_lab_based_method / method_source - computed over GROUP_KEY, after both streams exist
  all_final4 <- all_final4 %>%
    group_by(GROUP_KEY) %>%
    mutate(
      use_lab_based_method = max(days_since_first_sample, na.rm = TRUE) <= initial_days_val,
      method_source = if_else(use_lab_based_method, "Lab-based method", "Site-specific method")
    ) %>%
    ungroup()
  
  ## ---- PHASE 2 (apply_baseline_updates): mask by update_cond, forward-fill, pick stream ----
  
  if (target == "sars-cov-2") {
    all_final4 <- all_final4 %>%
      mutate(
        site_update_cond = coalesce(
          (days_since_first_sample <= 182) |
            (target_candidate & (days_to_nearest_target == site_min_days_to_target)),
          FALSE
        ),
        lab_update_cond = coalesce(
          (days_since_first_lab_sample <= 182) |
            (target_candidate & (days_to_nearest_target == lab_min_days_to_target)),
          FALSE
        )
      )
  } else {
    all_final4 <- all_final4 %>%
      mutate(
        site_update_cond = coalesce(
          (days_since_first_sample <= 365) |
            (target_candidate & (days_to_nearest_target == site_min_days_to_target)),
          FALSE
        ),
        lab_update_cond = coalesce(
          (days_since_first_lab_sample <= 365) |
            (target_candidate & (days_to_nearest_target == lab_min_days_to_target)),
          FALSE
        )
      )
  }
  
  # site stream: mask, then forward-fill within GROUP_KEY (needs date-sorted rows)
  all_final4 <- all_final4 %>%
    group_by(GROUP_KEY) %>%
    arrange(sample_collect_date, .by_group = TRUE) %>%
    mutate(
      site_baseline_ln = if_else(site_update_cond, site_roll_p10_24mo, NA_real_),
      site_std_ln = if_else(site_update_cond, site_roll_std_24mo, NA_real_),
      site_baseline_ln = na.locf(site_baseline_ln, na.rm = FALSE),
      site_std_ln = na.locf(site_std_ln, na.rm = FALSE)
    ) %>%
    ungroup()
  
  # lab stream: mask, then forward-fill within LAB_BASED_KEY (needs date-sorted rows)
  all_final4 <- all_final4 %>%
    group_by(LAB_BASED_KEY) %>%
    arrange(sample_collect_date, .by_group = TRUE) %>%
    mutate(
      lab_baseline_ln = if_else(lab_update_cond, lab_roll_p10_24mo, NA_real_),
      lab_std_ln = if_else(lab_update_cond, lab_roll_std_24mo, NA_real_),
      lab_baseline_ln = na.locf(lab_baseline_ln, na.rm = FALSE),
      lab_std_ln = na.locf(lab_std_ln, na.rm = FALSE)
    ) %>%
    ungroup()
  
  # pick the active stream per row - re-sorted by [GROUP_KEY, date], matching Python's final sort
  all_final4 <- all_final4 %>%
    group_by(GROUP_KEY) %>%
    arrange(sample_collect_date, .by_group = TRUE) %>%
    mutate(
      baseline_avg_conc_ln = if_else(use_lab_based_method, lab_baseline_ln, site_baseline_ln),
      std_avg_conc_ln = if_else(use_lab_based_method, lab_std_ln, site_std_ln),
      baseline_avg_conc_ln = if_else(is.na(baseline_avg_conc_ln), NaN, baseline_avg_conc_ln),
      std_avg_conc_ln = if_else(is.na(std_avg_conc_ln), NaN, std_avg_conc_ln)
    ) %>%
    ungroup() %>%
    select(-site_update_cond, -lab_update_cond) %>%
    drop_na(pcr_target)
  
  ## ---- PHASE 3 (finish_baseline): current baseline/std + wval, computed in-place ----
  ## Replaces the old last_vals dataframe + left_join entirely.
  
  result <- all_final4 %>%
    group_by(GROUP_KEY) %>%
    arrange(sample_collect_date, .by_group = TRUE) %>%
    mutate(
      current_baseline_ln = dplyr::last(baseline_avg_conc_ln),
      current_std_ln = dplyr::last(std_avg_conc_ln),
      current_baseline_ln = if_else(is.na(current_baseline_ln), NaN, current_baseline_ln),
      current_std_ln = if_else(is.na(current_std_ln), NaN, current_std_ln),
      days_site_online = max(days_since_first_sample, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      ww_index_ln_raw = if_else(
        !is.na(pcr_target_avg_conc_ln_winsor),
        (pcr_target_avg_conc_ln_winsor - current_baseline_ln) / current_std_ln,
        NaN
      ),
      .wval_bad = is.na(pcr_target_avg_conc_ln_winsor) | (current_std_ln == 0),
      ww_index_ln = if_else(.wval_bad, NA_real_, ww_index_ln_raw),
      ww_index_ln_lin = if_else(.wval_bad, NA_real_, exp(ww_index_ln_raw)),
      # bridge names so downstream site_weekly_average() doesn't need to change
      last_baseline_ln = current_baseline_ln,
      last_std_ln = current_std_ln
    ) %>%
    select(-ww_index_ln_raw, -.wval_bad) %>%
    group_by(GROUP_KEY, week_end) %>%
    arrange(sample_collect_date, .by_group = TRUE) %>%
    ungroup() %>%
    distinct()
  
  # create intermediate dataset
  intermediate_result <- result
  
  ##################################################
  ## MLM Overlap Handling
  ##################################################
  # if a site_id_with_pcr has more than one MLM reporting in a week,
  # keep data only if ww_index_ln_lin is non-null and it is from a more recent MLM
  result <- result %>%
    group_by(GROUP_KEY) %>%
    mutate(min_sample_date_for_mlm = min(sample_collect_date, na.rm = TRUE)) %>%
    ungroup() %>%
    # determines latest min_sample_date_for_mlm across all GROUP_KEY
    # that share the same site_id_with_pcr_source
    # using rank descending in case of a situation where these min dates happen to be the same (in which case they will both be 1)
    group_by(site_id_with_pcr, week_end) %>%
    arrange(desc(min_sample_date_for_mlm), desc(major_lab_method)) %>%
    mutate(latest_min_date_rank = rank(-as.numeric(min_sample_date_for_mlm), ties.method = "min")) %>%
    ungroup() %>%
    
    # filter to retain only rows from the latest GROUP_KEY with non-null ww_index_ln_lin
    filter(latest_min_date_rank == 1, !is.na(ww_index_ln_lin))
  
  
  #################################################
  # Categorize WVAL
  #################################################
  # at the result level
  unaggregated_result <- result %>%
    mutate(wval_result = categorize_wval(ww_index_ln_lin, pcr_target)) %>% distinct()
  
  unaggregated_result <- unaggregated_result
  
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
        #all_sample_ids = paste(sample_id, collapse = ", "),
        source = source,
        county_names = county_names, 
        all_method_source = paste(method_source, collapse = ", "), 
        all_lab_ids = paste(lab_id, collapse = ", ")
        
      ) %>% 
      # categorize wval
      mutate(wval_site_wk_level = categorize_wval(wval_site_wk, pcr_target)
      )%>% distinct()
  }
  
  
  
  ############################
  ## Return aggregated results 
  ############################            
  
  site_ag <- site_weekly_average(unaggregated_result)
  
  return(list(site_ag = site_ag, unaggregated_result = unaggregated_result))
  
  
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