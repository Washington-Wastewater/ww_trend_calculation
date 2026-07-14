## Purpose: provide functions for transforming redcap site data into a site key

## Load Libraries

library(tidyverse)

## REDCap sample instrument to df

# redcap_samples_to_df() transforms REDCap instrument format with repeating instruments
# to a dafaframe with repeating sample id rows for each of the pcr test results
# redcap_ww_samples - dataframe of redcap samples project output
# the output is two dataframes:
# redcap_all_samples - contains complete sample metadata and test results in each row 
# seq_results - contains sequencing results with identifying metadata

redcap_samples_to_df <- function(redcap_ww_samples){
  
  # Subset wastewater sample data to rows that are not part of repeating 
  # instruments. Remove columns that contain data for repeating instruments 
  # ww sample data from non-repeating instruments
  redcap_samples<- redcap_ww_samples %>% 
    filter(is.na(redcap_repeat_instrument)) %>%
    # remove columns
    select(-c(redcap_repeat_instrument,
              test_result_date,
              redcap_repeat_instance,
              starts_with("pcr"),
              lod_sewage,
              quality_flag,
              quality_flag_notes,
              nwss_submit,
              nwss_submit_notes,
              contains(c("variant", "complete", "timestamp", 
                         "received", "printed",
                         "orig")))) 
  
  # Subset the wastewater samples dataset to rows with the repeating instrument
  # data and only keep columns associated with the repeating instrument form
  # (c3_pcr_target_report_form)
  redcap_targets <- redcap_ww_samples %>% 
    filter(redcap_repeat_instrument =="c3_pcr_target_report_form") %>%
    select(sample_id, 
           test_result_date,
           starts_with("pcr"),
           lod_sewage,
           nwss_submit,
           nwss_submit_notes,
           quality_flag,
           quality_flag_notes)
  
  ## Create reformatted dataset 
  redcap_all_samples <- redcap_targets %>% 
    merge(redcap_samples, by = "sample_id") 

  # Seq data from sequence_proportions
  seq_results <- redcap_ww_samples %>% 
    filter(redcap_repeat_instrument =="sequence_proportions") %>%
    select(sample_id,
           variant_name,
           variant_prop)
  
  return(list(redcap_all_samples = redcap_all_samples,
              seq_results = seq_results))
}


## Map the REDCap samples data dictionary values

# map_sample_dictionary_values() maps REDCap pcr_target and pcr_gene_target data 
# dictionary values and sets data types
# redcap_all_samples - dataframe of redcap samples project after instrument transformation
# Output is redcap_targets_mapped - contains mapped values and desired data types
map_sample_dictionary_values <- function(redcap_all_samples){
  redcap_targets_mapped <- redcap_all_samples %>%
    mutate(pcr_target = case_match(pcr_target, 
                                   1 ~ "sars-cov-2",
                                   2 ~ "delta",
                                   3 ~ "omicron",
                                   4 ~ "hMPXV",
                                   5 ~ "FLUAV",
                                   6 ~ "FLUBV",
                                   7 ~ "RSV",
                                   13 ~ "FLUAV A H1",
                                   14 ~ "FLUAV A H3",
                                   15 ~ "FLUAV A H5",
                                   16 ~ "MeV_WT",
                                   17 ~ "HAV", 
                                   18 ~ "NoV GII"), 
         pcr_gene_target = case_match(pcr_gene_target,
                                      1 ~ "n1",
                                      2 ~ "n2",
                                      3 ~ "p681r",
                                      4 ~ "del156-157",
                                      5 ~ "wt214",
                                      6 ~ "del69/70",
                                      7 ~ "ins214epe",
                                      8 ~ "del142-144",
                                      9 ~ "wt156-157",
                                      10 ~ "e9l-nvar",
                                      11 ~ "G2R_G",
                                      12 ~ "InfA1",
                                      13 ~ "InfA2",
                                      14 ~ "InfB",
                                      15 ~ "n",
                                      16 ~ "InfA1 and InfA2 combined",
                                      17 ~ "RSV-A and RSV-B combined",
                                      23 ~ "INFA_H1pdm09 (CDC)",
                                      24 ~ "INFA_H3 (CDC)",
                                      25 ~ "INFA_H5 (CDC)",
                                      26 ~ "MeV_M_Gene_WT1",
                                      27 ~ "HAV_5' UTR", 
                                      28 ~ "NoV GII_ORF1-2 Junction"))%>% 
  # set data types
  mutate(sample_id = as.character(sample_id),
         sample_collect_date = as.Date(sample_collect_date)) 
  
  return(redcap_targets_mapped)
}


## Clean REDCap df

# redcap_df_clean() joins and filters REDCap results
# input_data - contains list of two objects
#     samples - contains complete sample metadata and test results in each row 
#     seq - contains sequencing results with identifying metadata
# the output is a list two dataframes:
#     redcap_all_samples_clean - contains complete sample metadata and test results in each row 
#     seq_results_clean - contains sequencing results with identifying metadata

redcap_df_clean <- function(input_data, population){
  
  samples <- input_data[[1]]
  seq <- input_data[[2]]
  
  
  redcap_all_samples_clean <- samples %>%
    # map the dictionary values
    map_sample_dictionary_values() %>%
    # rename the sample site name to match naming in the population df
    rename(site_id = sample_site_name) %>% 
    # join to add the population data
    left_join(population, by = "site_id")  %>% 
    ## Filter and keep data desired for data processing analysis
    # remove data missing sample_site_ids
    filter(!is.na(site)) %>% 
    # remove invalid sample data
    filter(pcr_target_avg_conc >= 0) %>%
    # add county columns
    mutate(county = COUNTY,
           county_main = COUNTY) %>%
    # group by site id and gene target
    group_by(site_id, pcr_gene_target, .drop = F) %>% 
    # keep only first ww record if there are two for a site on the same day
    distinct(sample_collect_date, .keep_all = TRUE) %>%
    # targets tested, with no positive detections will
    # have a missing conc value. To show no detection, change NA to 0
    mutate(pcr_target_avg_conc = if_else(is.na(pcr_target_avg_conc) & pcr_target_detection == "Not Detected {52}",
                                         0, pcr_target_avg_conc)) %>%
    ungroup()
  
  
  ## Clean sequencing dataset that has variables necessary from sequencing instrument
  missing_seq_vars <- redcap_all_samples_clean %>% 
    select(sample_id, 
           sample_collect_date,
           site,
           county) %>%
    unique()
  
  # Join metadata with seq results
  seq_results_clean <- seq %>% 
    left_join(missing_seq_vars, by = "sample_id")
  
  return(list(redcap_all_samples_clean = redcap_all_samples_clean, 
              seq_results_clean = seq_results_clean))
}
