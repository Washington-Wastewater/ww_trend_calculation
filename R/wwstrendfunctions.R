# WWS TREND ANALYSIS  -----------------------------------------------------






# PURPOSE -----------------------------------------------------------------

# 1. Functions to analyze trends in wastwater surveillance data for weekly report. 


# CONTRIBUTIONS -----------------------------------------------------------

# Wiley Jennings CDC/NCEZID/DFWED/WDPB
# Zach Marsh CDC/NCEZID/DFWED/WDPB
# Rasha Elnimeiry - Special Projects Data Requests, WA DOH
# Breanna McArdle - Environmental Public Health Science

# Dependencies ------------------------------------------------------------

# Depends on the data frame 'nwss' produced by calc-sars-cov2.R, which has been
# processed to contain not more than 1 observation per day per site (key_plot)


# Functions ---------------------------------------------------------------

make_df_for_trends <- function(.data, .date, .trend_var, .window_type ,
# make_df_for_trends <- function(.data, .date, .trend_var, .window_type = "time",
                               .window_size) {
  # Given:
  # .data (df): data frame of wastewater measurements for one group (key_plot)
  # .date (string): the column in .data giving dates to use for trends
  # .trend_var (string): the column in .data on which to calculate trends
  # .window_type (string): whether regression window should be defined in terms of days
  # - or measurements. Takes the value "time" for days or "measurements" for consecutive measurements.
  # Returns:
  # A data frame containing the data needed for trend regressions 
  
  # If calculating trends on measurements, can simply remove days without measurements
  if (.window_type == "measurements") {
    .data <- .data[!is.na(.data[[.trend_var]]), ]
  } else if (.window_type == "time") {
    # If calculating trends on days, need make .data contain a complete set of dates
    # Generate a complete set of dates. Also need to ensure that complete set of dates
    # is at least as long as the window_size; if not, need to add padding
    if( (Sys.Date() - min(.data[[.date]], na.rm = T) + 1) < .window_size) {
      dates_complete <- data.frame(
        date = seq(Sys.Date() - .window_size + 1, Sys.Date(), by = "1 day")
      )
    } else {
      dates_complete <- data.frame(
        date = seq(min(.data[[.date]], na.rm = T), Sys.Date(), by = "1 day")
      )
    }
    
    names(dates_complete) <- .date
    
    # Join .data to complete dates
    .data <- merge(dates_complete, .data, by = .date, all.x = TRUE)
  }
  
  # Return data frame
  .data
}



calc_trends <- function(.data, .date, .trend_var, .trend_var_variance, 
                        .window_size = 15, .window_type ) {
  
  # Given:
  # .data (df): data frame of wastewater measurements for one group (key_plot)
  # .date (string): the column in .data giving dates to use for trends
  # .trend_var (string): the column in .data on which to calculate trends
  # -Currently expects only a log10 scale variable with accompanying log10 scale variance
  # .trend_var_variance (string): the column in .date containing (log10 scale) variance estimates for .trend_var
  # .window_size (numeric): size of the window used for regression, specified in days (if
  # - if .window_type = "time") or number of measurements (if .window_type = "measurements")
  # .window_type (string): whether regression window should be defined in terms of days
  # - or measurements. Takes the value "time" for days or "measurements" for consecutive measurements.
  # Returns:
  # A data frame summarizing trend statistics (slope, p-values, percent daily 
  # change, total change, prediction upper bound for alerts, and whether data 
  # point constitutes alert) for each measurement 
  
  # Make data frame for trend calculations
  .data <- make_df_for_trends(.data = .data, .date = .date, .trend_var = .trend_var, 
                              .window_type = .window_type, .window_size = .window_size)
  
  # Set up matrix for storing trend summaries
  # -slopes: regression slope (equal to average daily concentration change); 
  # -pvals: p-value of slope; 
  # -pred_upr_exc: upper 90% prediction interval bound (i.e., 5% alpha, since
  # -- considering only 1-sided test)
  # -last_<trend_var>: the last concentration (or other trend var) measured
  last_meas <- paste0("last_", .trend_var)
  trend_sum_names <- c(
    "slopes", "pvals", "pdc", "ptc", "days_spanned", "preds_upr_exc", 
    "trend_num_meas", last_meas)
  trend_sum <- matrix(nrow = nrow(.data), ncol = length(trend_sum_names)) 
  colnames(trend_sum) <- trend_sum_names
  
  # Check for sufficient values present to compute trends; if not return NAs
  # For .window_type = "measurements", can simply make sure num measurements in 
  # entire data set is at least as large as .window_size 
  # For .window_type = "time", must check that there are sufficient measurements
  # in each window.
  if (.window_type == "measurements") {
    if (nrow(.data) < .window_size) {
      trend_sum <- as.data.frame(trend_sum)
      dates <- .data[[.date]]
      trend_sum$date <- dates
      trend_sum$alert <- NA_character_
      trend_sum <- trend_sum[, !(names(trend_sum) %in% last_meas)]  # remove last_meas var
      return(trend_sum)
    }
  }
  
  # Order data and create day number variable (since arbitrary date)
  .data <- .data[order(.data[[.date]]), ]
  .data[["day_num"]] <- as.numeric(.data[[.date]] - as.Date("2020-01-01"))
  
  # Assign weights of measurements, regression formula, summary variables to store
  # For now, only log10 values accepted, and trends calculated on log10 scale
  # KLJ - 
  .data[["weights"]] <- 1 #/.data[[.trend_var_variance]] 
  lmod_formula <- as.formula(paste(.trend_var, "~ day_num"))
  
  # Calculate num days spanned by trend (if .window_type = "time", it's constant)
  if (.window_type == "time") {
    days_spanned <- .window_size  # COULD MOVE THIS BEFORE FOR LOOP
  }
  
  # Compute and store trend variables
  lmod <- NULL
  # Define index for the beginning of the final window to be used
  ind_max <- nrow(.data) - (.window_size - 1) # Note: row number index can be used to index
  # by time as well as measurements because day_number index dates were forced to be complete for 
  # .window_type = "time"
  
  for (i in 1:ind_max) {
    # Define index of last measurement in measurement window for trend summary calculations
    ind_last <- i + .window_size - 1 
    
    lmod_prev <- lmod  # for prediction: meas period prior to most recent one
    
    #KLJ removed weights argument
    # Prep for regression
    lmod_data <- .data[i:ind_last, ]
    lmod_data <- lmod_data[!is.na(lmod_data[[.trend_var]]), ]
    # Calculate num days spanned (if .window_type = "time", it's already explicit)
    if (.window_type == "measurements") {
      days_spanned <- .data[ind_last, "day_num"] - .data[i, "day_num"] + 1  # inclusive interval
    }
    
    # Describe the trend window, regardless of whether calculate those trends.
    trend_sum[ind_last, "trend_num_meas"] <- nrow(lmod_data)
    trend_sum[ind_last, "days_spanned"] <- days_spanned  # inclusive interval
    
    # Compute regression. Notes on lm:
    # lm will silently toss NAs, so don't need to explicitly remove dates without measurements.
    # lm will calculate regression on 2 data points without complaining, so need
    # to check that at least 2 data points present. Easiest way
    # is just to filter out NAs (after lmod data frame is defined), and only run lm if
    # data are sufficient.
    if (nrow(lmod_data) > 1) {
      lmod <- lm(lmod_formula, data = lmod_data)  
      
      # Suppress warning that 'essentially perfect fit: 
      # summary may be unreliable', which occurs when all are non-detects
      lmod_sum <- withCallingHandlers(
        {summary(lmod)},
        suppress_this_warning = function(w) {
          if (grepl("essentially perfect fit", conditionMessage(w))) {
            invokeRestart("muffleWarning")}
        }
      )
      # Store trend vars from regression: vars that can be estimated from only 2 meas
      slope <- coef(lmod_sum)[2, 1]
      trend_sum[ind_last, "slopes"] <- slope
      trend_sum[ind_last, "pdc"] <- (10^slope - 1)*100  # percent daily change
      # Calculate Percent total change (ptc) from percent daily change
      trend_sum[ind_last, "ptc"] <- ((10^slope)^(days_spanned - 1) - 1) * 100
      
      # Store trend vars from regression: vars that need more than 2 meas
      if (nrow(lmod_data) > 2) {
        trend_sum[ind_last, "pvals"] <- coef(lmod_sum)[2, 4]
      }
      
      # Prediction interval is based on the measurement period preceding the most 
      # recent measurement period. Only execute if lmod has been estimated for previous 
      # window, AND a standard error was estimated (which requires >2 measurements).
      # Also only predict if slope is greater than ~machine level precision.
      trend_num_meas_prev <- trend_sum[ind_last - 1, "trend_num_meas"]
      
      if (!is.null(lmod_prev)) {
        if (abs(coef(lmod_prev)[2]) < 1e-10) {  #CHANGE
          trend_sum[ind_last, "preds_upr_exc"] <- NA_real_
        } else if (!is.na(trend_num_meas_prev) & trend_num_meas_prev > 2) {
          trend_sum[ind_last, "preds_upr_exc"] <- predict(
            lmod_prev, newdata = .data[ind_last, ], interval = "prediction",
            level = 0.90)[, "upr"]}
      }
    }    
  }
  
  trend_sum <- as.data.frame(trend_sum)
  
  # Add metadata to trend summaries; trend_sum has same num rows as .data, though
  # num rows differs by .window_type
  trend_sum[[last_meas]] <- as.numeric(.data[[.trend_var]]) 
  dates <- .data[[.date]]
  trend_sum$date <- dates
  
  # Add alert
  trend_sum$alert <- ifelse(
    is.na(trend_sum[["preds_upr_exc"]]), "no", ifelse(
      is.na(trend_sum[[last_meas]]), "no", ifelse(
        trend_sum[[last_meas]] > trend_sum[["preds_upr_exc"]], "yes", "no"
      )
    )
  )
  trend_sum <- trend_sum[, !(names(trend_sum) %in% last_meas)]  # no longer need last_meas var
  
  # Return
  trend_sum
}


# Function to classify trends calculated by calc_trends() 
classify_trends <- function(.trend_summary) {
  # Given:
  # .trend_summary: data frame of trend calculation summary vars. Currently,
  # - must include vars: slopes_sus (sustained slopes), pvals_sus (sustained 
  # - slope p values), slopes_short, pvals_short
  # Returns:
  # The .trend_summary data frame, with trend classified by the most recent
  # trend regression
  
  if (nrow(.trend_summary) == 0) {
    .trend_summary$trend <- NA_character_
    return(.trend_summary)
  }
  
  # Add 'trend' var, indicating trend for that measurement (not back-filled)
  mchn_zero <- 1e-10  # ~machine precision zero
  # .trend_summary$trend <- with(.trend_summary, ifelse(
  #   pvals_sus < 0.05 & slopes_sus < -mchn_zero & !is.na(slopes_sus), "sustained decrease", ifelse(
  #     pvals_sus < 0.05 & slopes_sus > mchn_zero & !is.na(slopes_sus), "sustained increase", ifelse(
  #       pvals_short < 0.05 & slopes_short < -mchn_zero & !is.na(slopes_short), "decrease", ifelse(
  #         pvals_short < 0.05 & slopes_short > mchn_zero & !is.na(slopes_short), "increase", ifelse(
  #           (abs(pvals_short) >= 0.05 | abs(slopes_short) <= mchn_zero) & !is.na(slopes_short), "plateau", NA_character_))))
  .trend_summary$trend <- with(.trend_summary, ifelse(
    pvals_sus < 0.1 & slopes_sus < -mchn_zero & !is.na(slopes_sus), "sustained decrease", ifelse(
      pvals_sus < 0.1 & slopes_sus > mchn_zero & !is.na(slopes_sus), "sustained increase", ifelse(
        pvals_short < 0.1 & slopes_short < -mchn_zero & !is.na(slopes_short), "decrease", ifelse(
          pvals_short < 0.1 & slopes_short > mchn_zero & !is.na(slopes_short), "increase", ifelse(
            (abs(pvals_short) >= 0.1 | abs(slopes_short) <= mchn_zero) & !is.na(slopes_short), "plateau", NA_character_))))
  ))
  
  # Return
  .trend_summary
}


# Function to both calculate and classify trends for one group for both short
# and sustained length trends, wrapping the calc_trends() and classify_trends() 
# functions
calc_classify_trends <- function(
  .data, .date, .trend_var, .trend_var_variance, .window_size_sus, 
  .window_size_short, .window_type) { 
  # Given:
  # -.data (df): data frame containing wastewater measurements for assessing 
  # --trends for one group (key_plot)
  # -.date (string): the column in .data giving dates to use for trends
  # -.trend_var (string): the column in .data on which to calculate trends
  # --Currently expects only a log10 scale variable with accompanying log10 scale variance
  # -.trend_var_variance (string): the column in .date containing (log10 scale) variance estimates for .trend_var
  # -.window_size_sus (numeric): size of window used to calculate sustained trends 
  # -.window_size_short (numeric): size of window used to calculated short trends
  # -.window_type (string): whether regression window should be defined in terms of days
  # --or measurements. Takes the value "time" for days or "measurements" for consecutive measurements.  
  # Returns:
  # A data frame, with trends calculated and classified for one group
  
  trend_sum_sus <- calc_trends(
    .data = .data, .date = .date, .trend_var = .trend_var,
    .trend_var_variance = .trend_var_variance,
    .window_type = .window_type, .window_size = .window_size_sus)
  trend_sum_short <- calc_trends(
    .data = .data, .date = .date, .trend_var = .trend_var, 
    .trend_var_variance = .trend_var_variance,
    .window_type = .window_type, .window_size = .window_size_short)
  
  # Join short and sustained trends
  trend_sum_short <- 
    trend_sum_short[, names(trend_sum_short) %in%
                      c("slopes", "pvals", "date", "pdc", "ptc", "days_spanned", 
                        "trend_num_meas", "alert", "preds_upr_exc")]
  trend_sum <- merge(trend_sum_sus, trend_sum_short, 
                     by = "date", all.x = TRUE, 
                     suffixes = c("_sus", "_short"))
  
  # Classify trends and return data frame
  classify_trends(.trend_summary = trend_sum)
}  




