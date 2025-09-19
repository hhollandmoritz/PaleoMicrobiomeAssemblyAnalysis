#' ## Prepare and combine Assembly Analysis results
#' This step takes the output from 01_calc_betaNTI.R and 02_calc_rcbc.R and runs
#' transforms them into useable formats for downstream analysis

#+ include=FALSE
# some setup options for outputing markdown files; feel free to ignore these
knitr::opts_chunk$set(eval = FALSE, 
                      include = TRUE, 
                      warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      dpi = 300,
                      fig.dim = c(9, 9),
                      out.width = '98%',
                      out.height = '98%')
#'
#+ include=TRUE
# Loading necessary packages and data
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(here)

# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("outputs")
figures.fp <- here("figures")

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}


# Set dataset
data_set <- "cyanos" # default is "all"
#### ====================================================================== ####

# Helpful functions
#### ====================================================================== ####
transform_betaNTI <- function(betaNTI_fp = here("outputs", "weighted_bNTI.csv"),
                              sample_metadata = input$sample_metadata){
  # Read in betaNTI results
  betaNTI <- read_csv(file = betaNTI_fp) %>%
    rename(SampleID = 1) %>%
    column_to_rownames(var = "SampleID")
  
  # Convert to long format
  ind <- which(lower.tri(betaNTI, diag = FALSE), arr.ind = TRUE)
  nn <- dimnames(betaNTI)
  betaNTI_long <- data.frame(row = nn[[1]][ind[, 1]],
                             col = nn[[2]][ind[, 2]],
                             val = betaNTI[ind])
  
  betaNTI_metadata <- betaNTI_long %>%
    rename(Site1 = 1, Site2 = 2, BetaNTI = 3) %>%
    left_join(sample_metadata, by = c("Site1" = "SampleID")) %>%
    left_join(sample_metadata, by = c("Site2" = "SampleID"), suffix = c(".Site1", ".Site2")) %>%
    # Add column of BetaNTI interpretation
    mutate(Assembly_Process = ifelse(BetaNTI < -2, "Homogenous selection",
                                     ifelse(BetaNTI > 2, "Heterogenous selection", 
                                            "Stochastic")))
  
  return(betaNTI_metadata)  
  
}

#### ====================================================================== ####

#### Read in data 
#### ====================================================================== ####
# Select data:
input <- choose_data_set(data_set = data_set)

# Transform BetaNTI from different sources
betaNTI_all <- transform_betaNTI(betaNTI_fp = here("outputs", paste0("weighted_bNTI_", data_set, ".csv")),
                                 sample_metadata = input$sample_metadata)

betaNTI <- read_csv(file = here("outputs", paste0("weighted_bNTI_", data_set, ".csv"))) %>%
  rename(SampleID = 1) %>%
  column_to_rownames(var = "SampleID")
# # Create a long-form dummy data frame to combine results
# samp.df <- data.frame(Site1 = colnames(betaNTI_sing), Site2 = rownames(betaNTI_sing))
# 
# 
# betaNTI_sing[upper.tri(betaNTI_sing)] <- 999 
# 
# Read in RCBC results
RCBC.raw <- read_csv(file = here("outputs", paste0("/rcbc_matrix_",data_set, ".csv")))

RCBC <- RCBC.raw[c(as.character(1:ncol(input$otu_table[-1]))),
                 c(as.character(1:ncol(input$otu_table[-1])))]
RCBC <- as.matrix(RCBC)
colnames(RCBC) <- names(input$otu_table[-1])
rownames(RCBC) <- names(input$otu_table[-1])
# 
# # Set the upper corner of RCBC to 999 to distinguish self-comparison NAs from 
# # NAs generated during RCBC calculation; 
# # This is important if you set speedup = TRUE because there will be NAs in the 
# # matrix where no comparison was performed
# RCBC[upper.tri(RCBC)] <- 999 
# 
# # Read in long format of rcbc
#RCBC.lf <- read_csv(file = paste0(outputs.fp, "/rcbc_long_form.csv"))

#### ====================================================================== ####


# Filter out comparisons we won't use
# Placeholder - so far planning to use all comparisons
#### ====================================================================== ####
# Since we're working with depths, let's make a list of all consecutive comparisons and only filter those out
Ranked_depth.df <- input$sample_metadata %>%
  select(SampleID, AvgDepth) %>%
  mutate(DepthRank = rank(AvgDepth),
         SampleID2 = lead(SampleID)) %>%
  rename(Site2 = SampleID,
         Site1 = SampleID2) # this is reversed because of the way Site1 and 2 are in the betaNTI all




# All comparisons within site
betaNTI_all_filt <- betaNTI_all %>%
  mutate(TempConditionType = paste0(`Temperature condition.Site1`, "-", `Temperature condition.Site2`)) %>%
  mutate(TempConditionType = ifelse(TempConditionType == "Cold-Warm", "Warm-Cold", TempConditionType)) %>%
  mutate(EpochType = case_when(`Climate epoch.Site1` == `Climate epoch.Site2` ~ `Climate epoch.Site1`,
                               .default = paste0(`Climate epoch.Site1`, ":", `Climate epoch.Site2`))) %>%
  # Filter out Pre-LGS:Holocene comparisons because it violates rules of time
  filter(EpochType != "Pre-LGS:Holocene") %>%
  select(Site1, Site2, BetaNTI, Assembly_Process, ends_with("Type")) %>%
  mutate(Assembly_Process_expl = case_when(Assembly_Process == "Homogenous selection" ~ "Homogenous (abiotic) selection",
                                           Assembly_Process == "Heterogenous selection" ~ "Heterogenous (biotic) selection", 
                                           Assembly_Process == "Stochastic" ~ "Stochastic"))


# All comparisons within site
betaNTI_consecutive_only <- betaNTI_all %>%
  mutate(TempConditionType = paste0(`Temperature condition.Site1`, "-", `Temperature condition.Site2`)) %>%
  mutate(TempConditionType = ifelse(TempConditionType == "Cold-Warm", "Warm-Cold", TempConditionType)) %>%
  mutate(EpochType = case_when(`Climate epoch.Site1` == `Climate epoch.Site2` ~ `Climate epoch.Site1`,
                               .default = paste0(`Climate epoch.Site1`, ":", `Climate epoch.Site2`))) %>%
  mutate(Assembly_Process_expl = case_when(Assembly_Process == "Homogenous selection" ~ "Homogenous (abiotic) selection",
                                           Assembly_Process == "Heterogenous selection" ~ "Heterogenous (biotic) selection", 
                                           Assembly_Process == "Stochastic" ~ "Stochastic")) %>%
  right_join(Ranked_depth.df, by = c("Site1", "Site2")) %>%
  filter(!is.na(Site1))

#### ====================================================================== ####
# Plot comparisons
#### ====================================================================== ####
# # Compare the assembly processes for cold vs warm years
# ggplot(betaNTI_all_filt, aes(x = ComparisonType, y = BetaNTI)) +
#   geom_point(aes(color =Assembly_Process_expl), position = "jitter") +
#   theme_bw()



#### ====================================================================== ####
# Transform betaNTI and RCBC results to long format
#### ====================================================================== ####
betaNTI.lf <- betaNTI_all

# Sanity check, do we have the right number of rows?
# There are 528 combinations (without self-comparisons) of 33 samples
nrow(betaNTI.lf) == dim(combn(nrow(input$sample_metadata), 2))[2] # TRUE == good (combn(33, 2))

# Get list of "stochastic" comparisons, by setting any
# non-stochastic comparisons to NA
rcbc <- as.matrix(betaNTI)
rcbc[betaNTI < -2 | betaNTI > 2] <- NA
rcbc[which(!is.na(rcbc))]
image(as.matrix(rcbc)) # lots of stochasticity
image(as.matrix(RCBC))
image(as.matrix(betaNTI))

# RCBC
## -2 < betaNTI < 2 and RCBC < -0.95 indicates homogenizing dispersal
## -2 < betaNTI < 2 and RCBC > 0.95 indicates dispersal limitation and drift
## -2 < betaNTI < 2 and -0.95 < RCBC < 0.95 indicates drift
RCBC.lf <- RCBC %>%
  as.data.frame() %>%
  rownames_to_column(var = "Site1") %>% 
  pivot_longer(cols = !Site1, names_to = "Site2",
               values_to = "RCBC") %>%
  # remove duplicated comparisons
  filter(Site1 != Site2) %>% # self-comparisons
  filter(!is.na(RCBC)) %>% # duplicate comparisons from upper triangle
  # remove samples that are not in the metadata
  filter(Site1 %in% input$sample_metadata$SampleID) %>%
  filter(Site2 %in% input$sample_metadata$SampleID) %>%
  # Add column of BetaNTI interpretation
  mutate(Assembly_Process = ifelse(RCBC < -0.95, "Homogenizing dispersal",
                                   ifelse(RCBC <= 0.95, "Drift",
                                          ifelse(RCBC > 0.95, "Dispersal limitation and drift",
                                                 NA))))

# Sanity check, do we have the right number of rows?
# ncol(combn(xxxsamples, 2))
nrow(RCBC.lf) == ncol(combn(nrow(input$sample_metadata), 2))
#### ====================================================================== ####

#' Create a table with betaNTI and rcbc results combined
#### ====================================================================== ####
betanull_consecutive.lf <- left_join(betaNTI_consecutive_only %>%
                           mutate(Assembly_Process = ifelse(Assembly_Process == 
                                                              "Stochastic",
                                                            NA, Assembly_Process)),
                         RCBC.lf, by = c("Site1", "Site2")) %>% 
  mutate(Assembly_Process = coalesce(Assembly_Process.x, Assembly_Process.y),
         StochasticDeterministic = ifelse(grepl("selection", Assembly_Process), 
                                          "Deterministic",
                                          "Stochastic"),
         Assembly_Process = factor(Assembly_Process, 
                                   levels = c("Homogenous selection",
                                              "Heterogenous selection",
                                              "Homogenizing dispersal", 
                                              "Dispersal limitation and drift", 
                                              "Drift"))) %>%
  rename(Assembly_Process.RCBC = Assembly_Process.y,
         Assembly_Process.BetaNTI = Assembly_Process.x) %>%
  mutate(Assembly_Process.RCBC = ifelse(grepl("selection", 
                                              Assembly_Process.BetaNTI), 
                                        NA, Assembly_Process.RCBC),
         RCBC.nona = RCBC,
         RCBC = ifelse(grepl("selection", Assembly_Process.BetaNTI), 
                       NA, RCBC)) %>%
  dplyr::select(Site1, Site2, BetaNTI, RCBC, Assembly_Process.BetaNTI, 
         Assembly_Process.RCBC, Assembly_Process, everything())

betanull.lf <- left_join(betaNTI_all_filt %>%
                                       mutate(Assembly_Process = ifelse(Assembly_Process == 
                                                                          "Stochastic",
                                                                        NA, Assembly_Process)),
                                     RCBC.lf, by = c("Site1", "Site2")) %>% 
  mutate(Assembly_Process = coalesce(Assembly_Process.x, Assembly_Process.y),
         StochasticDeterministic = ifelse(grepl("selection", Assembly_Process), 
                                          "Deterministic",
                                          "Stochastic"),
         Assembly_Process = factor(Assembly_Process, 
                                   levels = c("Homogenous selection",
                                              "Heterogenous selection",
                                              "Homogenizing dispersal", 
                                              "Dispersal limitation and drift", 
                                              "Drift"))) %>%
  rename(Assembly_Process.RCBC = Assembly_Process.y,
         Assembly_Process.BetaNTI = Assembly_Process.x) %>%
  mutate(Assembly_Process.RCBC = ifelse(grepl("selection", 
                                              Assembly_Process.BetaNTI), 
                                        NA, Assembly_Process.RCBC),
         RCBC.nona = RCBC,
         RCBC = ifelse(grepl("selection", Assembly_Process.BetaNTI), 
                       NA, RCBC)) %>%
  dplyr::select(Site1, Site2, BetaNTI, RCBC, Assembly_Process.BetaNTI, 
                Assembly_Process.RCBC, Assembly_Process, everything())

# Sanity check, do we have the right number of rows?
nrow(betanull.lf) == ncol(combn(nrow(input$sample_metadata), 2))


# Transform combined dataframe to wide format
betanull.wf <- betanull.lf %>%
  select(Site1, Site2, Assembly_Process, BetaNTI, RCBC) %>%
  # recode Assembly processes
  mutate(Assembly_Process_Code = factor(Assembly_Process, 
                                        levels = c("Heterogenous selection",
                                                   "Dispersal limitation and drift",
                                                   "Drift",
                                                   "Homogenizing dispersal",
                                                   "Homogenous selection")),
         Assembly_Process_level = as.numeric(Assembly_Process_Code)) %>%
  dplyr::select(Site1, Site2, Assembly_Process_level) %>%
  pivot_wider(names_from = Site1, values_from = Assembly_Process_level) %>% 
  column_to_rownames(var = "Site2")

# Fill in upper diagonal
diag(betanull.wf) <- NA
betanull.mat <- as.matrix(betanull.wf)

#### ====================================================================== ####

# Save outputs
#### ====================================================================== ####
# # Data
# assembly results in long format
saveRDS(betanull_consecutive.lf, paste0(outputs.fp, "/betanull.consecutive.lf_", data_set, ".RDS"))
write.csv(betanull_consecutive.lf,
          paste0(outputs.fp, "/betanull.consecutive.lf_", data_set, ".csv"),
          quote=TRUE, row.names = FALSE)

saveRDS(betanull.lf, paste0(outputs.fp, "/betanull.lf_", data_set, ".RDS"))
write.csv(betanull.lf,
          paste0(outputs.fp, "/betanull.lf_", data_set, ".csv"),
          quote=TRUE, row.names = FALSE)

saveRDS(betanull.wf, paste0(outputs.fp, "/betanull.wf_", data_set, ".RDS"))
write.csv(betanull.wf,
          paste0(outputs.fp, "/betanull.wf_", data_set, ".csv"),
          quote=FALSE, row.names = FALSE)

