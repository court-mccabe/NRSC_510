# Update file path
file_path <- "C:/NRSC_510/IPV_cohort_simulated_v5.xlsx"
#file path sets excel so R can read from my computer

#combine data
acute_raw   <- read_excel(file_path, sheet = "Acute_Subacute")
chronic_raw <- read_excel(file_path, sheet = "Chronic")
cohort_raw <- bind_rows(acute_raw, chronic_raw)
# excel reads each sheet separately. use read_excel to read only the sheet
# then create new variable to combine individuals rows to complete dataset
#capturing both data sets. bind_rows stacks vertically, bind_cols stacks
# horizontally (makes wider dataset)


#neuropsych data
neuro_test_cols <- c(
  "RAVLT_total", "RAVLT_delayed_recall",
  "TMT_A_Tscore", "TMT_B_Tscore",
  "TOPF_standard_score", "COWA_total_words",
  "SDMT_oral_z", "SDMT_written_z",
  "WAIS_PSI", "WAIS_WMI",
  "TOMM_trial1", "TOMM_trial2",
  "PCL5_total", "PHQ9_total", "NSI_total",
  "BISA_BI_count", "HELPS_total"
)
#manually written character vector that match column names -- will reference
#later

#identify characters
char_missing_cols <- neuro_test_cols[
  sapply(cohort_raw[neuro_test_cols], is.character)
]
# sapply runs is.character on each column in neuro_test_cols within cohort_raw, 
#returning a named logical vector (TRUE/FALSE per column). The outer 
#neuro_test_cols[...] then subsets only the names where the result was 
#TRUE — giving us only the columns that are stored as character type, 
#which indicates where “N/A” was manually entered as text instead of true missing data.


# Apply na.if
cohort <- cohort_raw %>%
  mutate(across(all_of(char_missing_cols), ~ na_if(.x, "N/A")))
#we manually input missing data as "N/A" but need to replace that we NA in R
# which is coded as missing data. apply mutate function across (this performs 
#function across multiple columns at once) and all of the character vectors
#all_of acts as a helper function inside of across to use the list of columns
#provided > .x is for each element in each column (so with all_of and .x
# we cover all columns and all elements to be checked for "N/A" and be 
#replaced with NA)


#across() is the engine — it tells mutate to apply a function to multiple 
#columns at once, and takes a column selection and a function as its two arguments

#all_of() is a selection helper used inside across() — it tells across 
#“use this external character vector as the list of columns.” It’s needed 
#because inside across(), R doesn’t automatically know that a variable like 
#yes_no_cols refers to column names rather than data

#Create clean group labels
cohort <- cohort %>%
  mutate(
    study_group = case_when(
      phase == "Acute/subacute" & group == "IPV-BI"     ~ "Acute IPV-BI",
      phase == "Acute/subacute" & group == "Non-IPV-BI" ~ "Acute Non-IPV-BI",
      phase == "Chronic"        & group == "IPV-BI"     ~ "Chronic IPV-BI"
    ),
    study_group = factor(study_group,
                         levels = c("Acute IPV-BI",
                                    "Acute Non-IPV-BI",
                                    "Chronic IPV-BI"))
  )

# within the mutate function, we assign a case_when where the phase is either 
# acute/subacute, "IPV-BI" or non IPV-BI. Both conditions have to be true for a 
# paticipants to be assigned to a group. If a participant does not meet the 
# conditions, they are assigned NA. Then we use factor to give rank to the 
# groups, this is not a weighted rank because we do not have ordered = TRUE 
#it is just so it displays in a way where the acute groups are together. 
# This could be rearranged to show IPV-BI groups together too
#without factor R will automatically assign them alphabetically

#identify binary variables
yes_no_cols <- cohort %>%
  select(
    starts_with("incident_"),
    starts_with("dx_")
  ) %>%
  names()

#we use a piping format into the select format (select returns a new dataframe
#with only the variables that match the criteria set) and starts_with will choose any 
#variables within cohort that start with this wording. then we call the 
#names of these variables. Names() extracts the column names of this dataframe
#as a character vector

cohort <- cohort %>%
  mutate(across(
    all_of(yes_no_cols),
    ~ case_when(
      .x == "Yes" ~ TRUE,
      .x == "No"  ~ FALSE,
      TRUE         ~ NA
    )
  ))

# we are assigning more meaning to the binary columns by using the mutate 
#function. We using across to work through every variable. then we use a 
# case_when function to use .x which goes through each element and applies 
# yes to true, no to false, and true to NA, so anything leftover would be NA
#this is important because it allows us to use T/F in statistical models
#which yes/no is limiting in