# Wastewater Trend Calculation Code

The purpose of this code is to process wastewater testing results and calculate trends in data.


# R Files

**WAWBE_trend_calculation_share_code.Rmd**

- This file contains the code to pull data from redcap and calculate wastewater trends. 
- Inputs: sewershed population data, REDCap wastewater sample data
- Outputs: core dataset (analyzed wastewater results)


**wwstrendfunctions.R**

- This file contains trend calculation functions provided by NWSS that are used in WAWBE_trend_calculation_share_code.Rmd

**credentials_ref.R**

- This is a reference file for the credentials.R file that is necessary to run the script. The purpose of this file is to store the REDCap API token (password)


# Usage

This code should be re-run every time new wastewater results are available.

# General Notes

Sewershed population data and REDCap access is necessary to run this script. Please contact the WAWBE for access.

# Contact Information
To provide feedback, ask questions, or inquire about contributing to this repository, please contact wastewaterbasedepi@doh.wa.gov.

