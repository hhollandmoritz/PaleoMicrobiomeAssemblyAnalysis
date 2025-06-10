#' ## Calculate betaNTI
#' This step calcluates the betaNTI or phylogenetic turnover among samples
#' It requires a phylogenetic tree and a species-by-sample table
#' This step takes a long time and would be much better run on a super computer

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
library(picante); packageVersion("picante") # handles tree manipulations
library(vegan); packageVersion("vegan") # for ecological applications
library(here)
library(parallel)
library(foreach)
library(doMPI) # For running on super computer

# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("outputs")
figures.fp <- here("figures")

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

#### =================================================================================================== ####
#' First we'll set up user settings and command line arguments: 
#' Users should modify these settings for their own environment.  
#'  - Set the number of replicate randomizations for the BetaMNTD calculations `beta.reps`
#'  - Set the relative abundance flag if your OTU table is not integers but real numbers instead
#'  - If you're working in an environment with parallel computing set `paral <- TRUE`
# User settings and command line arguments
#### =================================================================================================== ####
beta.reps <- 1000; # number of randomizations for beta MNTD

rel_abund <- FALSE # is the otu table in relative abundances such as read mappings?

paral <- TRUE

# Tell the user their settings:
writeLines(paste("User settings are as follows: \n",
                 "Parallel:", paral, "\n",
                 "Replicates:", beta.reps, "\n",
                 "Relative Abundances:", rel_abund))

#Install and load all required packages

# If parallel is true, detect available cores
if(paral==T) {
  # Figure out parallelization engine; ideally MPI, if not mpi, parallel, if not parallel, none
  if(require(Rmpi) & require(doMPI)) {
    n.cores <- max(1, mpi.universe.size()-1) # set n.cores to 1 or mpi allocation size, whichever is bigger
  } else if(require(parallel)) { # else if parallel is available use it to calculate cores
    n.cores <- max(1, detectCores() - 1)
  } else {
    writeLines("No parallization pacakge detected, setting n.cores to 1")
    n.cores <- 1
  }
  
  #print(n.cores)
  print(paste(n.cores, "worker nodes are available."))
  
  if (n.cores > 24) {
    # we are using more than 1 node, so really run in parallel mode
    cl <- startMPIcluster(count = n.cores)
    registerDoMPI(cl) # tell foreach about the cluster
    print(paste("Running in parallel mode on",n.cores,"worker nodes."))
    
  } else if (n.cores > 1) {
    library(doParallel)
    cl <- parallel::makeForkCluster(n.cores)
    registerDoParallel(cl)
    print(paste("Running in parallel mode without mpi on", n.cores, "worker nodes."))
  } else {
    foreach::registerDoSEQ() # tells foreach to use sequential mode
    paral <- FALSE
    print(paste("Only", n.cores,
                "available. Not enough to run in parallel, running in sequential mode."))
  }
} else {
  foreach::registerDoSEQ()
  print("Running in sequential mode.")
} 

#### =================================================================================================== ####


#' ## Read in data 
#### =================================================================================================== ####
# Read in OTU table and env data
writeLines("reading in input data")

# Read in tree
otutree <- input_all$tree

# Verify the everything is matching
match.phylo.otu <- match.phylo.data(otutree, input_all$otu_table[-1])

#### =================================================================================================== ####

#' ## Step 1: First we calculate the empirical betaMNTD
#' This is the observed betaMNTD which we will compare to our null distribution that we will generate below
#### =================================================================================================== ####
tip_dist <- cophenetic(match.phylo.otu$phy) # calculate this separately so that if you want, you can save it

beta.mntd.weighted <- as.matrix(comdistnt(t(match.phylo.otu$data),
                                          tip_dist,
                                          abundance.weighted=T))
dim(beta.mntd.weighted)


#### =================================================================================================== ####

#' ## Step 2: Calculate null expectation (randomized) betaMNTD
#' This is our null distribution. Calculating this will take a while - best to run on server
#### =================================================================================================== ####
null_bMNTD <- NULL # make a null object to hold list of results
null_bMNTD <- foreach (i=1:beta.reps, .packages = c("picante")) %dopar% {
  as.matrix(comdistnt(t(match.phylo.otu$data),
                      taxaShuffle(cophenetic(match.phylo.otu$phy)),
                      abundance.weighted=T,
                      exclude.conspecifics = F))
  } # end of parallel computation

# reformat data output from foreach
# Unlisting the list of null beta.reps; This creates an 3d array
# of size: # samples x # samples x # beta.reps
rand.weighted.bMNTD.comp <- array(as.numeric(unlist(null_bMNTD)), 
                                  dim=c(ncol(match.phylo.otu$data),
                                        ncol(match.phylo.otu$data), 
                                        beta.reps))

dim(rand.weighted.bMNTD.comp) # should be # samples x # samples x # beta.reps

#### =================================================================================================== ####

#' ## Step 3: Now calculate betaNTI by normalizing observed against the null expectation
#' Reminder: 
#' -2 < bNTI < 2 indicates stochastic process which can be further refined using values from RCbray in next step.
#' bNTI < -2 indicates homogeneous selection
#' bNTI > 2 indicates heterogeneous selection
#### =================================================================================================== ####
# Save time by setting up matrix first
weighted.bNTI <- matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));

dim(weighted.bNTI); # should be # samples x # samples

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) { # columns + 1 so samples are not compared to themselves
    
    # Pull out the randomized betaMNTD values
    rand.vals <- rand.weighted.bMNTD.comp[rows,columns,]; # length of this object = beta.reps
    # Normalize the observed against the randomized values
    weighted.bNTI[rows,columns] <- (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals); # normalize
    rm("rand.vals");
    
  };
};

# Rename columns and rows to get a samplexsample table
rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
# weighted.bNTI;
#### =================================================================================================== ####

#' Step 4: Verification - Plot Null and empirical distributions
#### =================================================================================================== ####
dist_plot.df <- map_df(null_bMNTD, ~{as.data.frame(.x)}, .id = "null_id") %>% 
  mutate(null_id = paste0("null_", null_id)) %>% 
  bind_rows(.,beta.mntd.weighted %>%
              as_tibble() %>%
              mutate(null_id = "empirical")) %>%
  #group_by(null_id) %>%
  #nest()
  mutate(null_id = as.factor(null_id)) %>%
  pivot_longer(cols = !matches("null_id"),
               names_to = "rep", values_to = "betaMNTD")

dist_plot <- dist_plot.df %>%
  filter(null_id != "empirical") %>%
  ggplot(aes(x = betaMNTD, fill = null_id)) +
  geom_density(alpha = 0.3) +
  scale_fill_grey() +
  geom_density(data = dist_plot.df %>%
                 filter(null_id == "empirical"),
               aes(x = betaMNTD),
               fill = "red", alpha = 0.3) +
  theme(legend.position = "none")

#### =================================================================================================== ####

#### Save Data and Figures
#### =================================================================================================== ####

# Data
# BetaNTI Matrix
write.csv(weighted.bNTI,
          paste0(outputs.fp, "/weighted_bNTI.csv"),
          quote=FALSE)

# Figures
ggsave(dist_plot, device = "png",
       filename = paste0(figures.fp, "/weighted_bNTI_dist.png"))
