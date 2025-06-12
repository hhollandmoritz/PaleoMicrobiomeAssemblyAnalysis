#' ## Interpret Assembly Analysis results for Paleo Microbiome Project
#' ***This step takes the output from 
#' 03_prepare_and_combine_assembly_analysis_results.R and runs initial figures 
#' and analyses on them. It has been knit into a notebook for ease of interpretation***

#' # Driving questions:
#'  1. What assembly processes predominate in these paleo-microbiome samples?
#'  2. What assembly processes characterize each Epoch, and the Epoch transitions?
#'  3. What assembly processes drive cold to warm transitions?
#'  4. Do those processes differ in different Climate Epochs?
#'  Bonus: Are any other factors related to assembly processes?

#' #### How does assembly analysis work?
#' Find out more at this link: [SlideShow](https://docs.google.com/presentation/d/1BSLtMNZZrXxR9Nk0GER2RBBGPU0I1v4zAgnr9jylAQw/edit?usp=sharing)


#+ include=FALSE
# some setup options for outputing markdown files; feel free to ignore these
knitr::opts_chunk$set(eval = TRUE, 
                      echo = TRUE,
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
library(vegan); packageVersion("vegan") # for ecological applications
library(viridis)
library(cowplot) # Pretty plotting
library(here)

# Load required data
source(here("setup.R"))

## Set up input and output directories
outputs.fp <- here("outputs")
figures.fp <- here("figures")

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

#+ echo=F

# Setup plotting parameters (this is hidden for ease of interpretability)
#### ====================================================================== ####
colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), c("#737373", "#FFFFFF", "#000000")), c("blue", "darkblue", "green", "darkgreen", "red", "darkred", "orange", "darkorange", "purple", "darkpurple", "yellow", "brown", "grey", "white", "black"))
type_levels <- c("BL", "RF")
colour_type <- c("#703C1B","brown", "#058000")
fill_type   <- colour_type
assembly_levels <- c("Homogenous selection", "Heterogenous selection",
                     "Homogenizing dispersal", "Dispersal limitation and drift",
                     "Drift")
assembly_labels <- assembly_levels
#colour_assembly <- c("#133253", "#738CA6","#1C4A00","#80BA5D", "#B1DF95") # blue/green
#colour_assembly <- c("#3F002E", "#BE7FAD","#871200","#DA3015", "#FF816D") # purple/red
#colour_assembly <- c("#6C358D", "#AB81C4","#085A4F","#2D8478", "#79BDB4") # purple/teal
colour_assembly <- c("#521168", "#8e318f","#BD4D0C","#FF8000", "#FADBAC") # purple/orange
fill_assembly <- colour_assembly
#### ====================================================================== ####



#'### Read in data not included in `setup.R`

#' Some important data processing notes: 
#' - Comparisons between Pre-LGS and Holocene samples have been filtered from the dataset. The logic behind this is that it doesn't make much sense to compare samples that are not consecutive in time.
#' - Epochs are ordered by time
#' - Temperature comparisons are ordered by "same-same", "different"
#' - After filtering out the Pre-LGS:Holocene comparisons, we are left with 468 comparisons

# Read in assembly analysis results
betanull.lf <- read_csv(file = paste0(outputs.fp, "/betanull.lf.csv"))
betanull_consecutive.lf <- read_csv(file = paste0(outputs.fp, "/betanull.consecutive.lf.csv"))


# Change variables into a factors
betanull.lf <- betanull.lf %>%
  mutate(across(starts_with("TempCondition"), ~ factor(.x, 
                                               levels = c("Cold-Cold", 
                                                          "Warm-Warm", 
                                                          "Warm-Cold")))) %>%
  mutate(EpochType = factor(EpochType, 
                            levels = c("Holocene", 
                                       "LGS:Holocene", 
                                       "LGS", 
                                       "Pre-LGS:LGS", 
                                       "Pre-LGS"))) %>%
  mutate(Assembly_Process = factor(Assembly_Process, 
                                   levels = c("Homogenous selection",
                                              "Heterogenous selection",
                                              "Homogenizing dispersal",
                                              "Dispersal limitation and drift",
                                              "Drift")))

# Change variables into a factors
betanull_consecutive.lf <- betanull_consecutive.lf %>%
  mutate(across(starts_with("TempCondition"), ~ factor(.x, 
                                                       levels = c("Cold-Cold", 
                                                                  "Warm-Warm", 
                                                                  "Warm-Cold")))) %>%
  mutate(EpochType = factor(EpochType, 
                            levels = c("Holocene", 
                                       "LGS:Holocene", 
                                       "LGS", 
                                       "Pre-LGS:LGS", 
                                       "Pre-LGS"))) %>%
  mutate(Assembly_Process = factor(Assembly_Process, 
                                   levels = c("Homogenous selection",
                                              "Heterogenous selection",
                                              "Homogenizing dispersal",
                                              "Dispersal limitation and drift",
                                              "Drift")))

#'  ### Question 1: What assembly processes predominate in these paleo-microbiome samples?

#+ eval=TRUE
# Calculate proportions of pairwise comparisons overall
#### ====================================================================== ####
# Proportions across all samples
betanull.lf %>%
  select(Site1, Site2, BetaNTI, RCBC, Assembly_Process) %>%
  group_by(Assembly_Process) %>%
  tally() %>%
  mutate(Total = sum(n),
         Percent = round(100*n/Total, digits = 2)) %>%
  arrange(desc(Percent)) %>% knitr::kable()

#' - Most of the pairwise comparisons are characterized by "Homogenous selection". This is typically considered to be the result of abiotic selection that favors some clades over others, rather than filtering sister taxa. 
#' - The second most common assembly process is ecological drift. This is a purely stochastic process. Assembly cannot be attributed to selection or dispersal. 
#' - We see essentially no evidence (only 2/468 comparisons) of heterogeneous (biotic) selection
#'   
#'   
#' ### Question 2: What assembly processes characterize each Epoch, and the Epoch transitions?
#' 
#' 
#+ eval=TRUE
# Summary plot by type
epochtype.plot.prop.df <- betanull.lf %>%
  select(Site1, Site2, Assembly_Process, EpochType) %>%
  group_by(EpochType, Assembly_Process) %>%
  tally() %>% 
  ungroup() %>% group_by(EpochType) %>%
  mutate(Total = sum(n),
         Percent = 100*n/Total)
  
epochtype.plot.prop <- epochtype.plot.prop.df %>%
  ggplot(aes(x= EpochType, y = Percent, 
             fill = Assembly_Process)) +
  facet_wrap(~EpochType, shrink = TRUE, drop = TRUE, nrow = 1,
             scales = "free_x") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Assembly Process", values = fill_assembly, 
                     breaks = assembly_levels, 
                     labels = assembly_labels) +
  xlab("") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  #coord_cartesian(expand = c(0,0)) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(1.5)),
        panel.spacing.x = unit(0, "lines"),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = rel(1)),
        strip.placement = "outside",
        legend.position = "bottom")

epochtype.plot.prop.legend <- get_legend(epochtype.plot.prop)
#type.plot.prop <- type.plot.prop + theme(legend.position = "none")

ggsave(epochtype.plot.prop, 
       filename = paste0(figures.fp, "/epochType_prop.png"),
       width = 10, height = 10, dpi = 400)
epochtype.plot.prop
#' - As time proceeds, we see an increase in the proportion of pariwise comparisons characterized by drift in our dataset.
#' - This pattern is often observed with depth in soils, more generally. It could be due to relic DNA (which experiences no selection), or it could be indicative that in previous epochs there were less selective pressures on the organisms forming these communities. 
#' - Interestingly, we don't see major increases in selctive processes at the Epoch-Epoch tansistions. 
#' - But, dipsersal limitation seems to have played a fairly strong role in the Pre-LGS:LGS transition. 
#' - The modern era communities are not at all characterized by dispersal, but there appears to be more of it in the LGS, and Pre-LGS periods. 
#' - This may suggest greater wind or water-mediated dispersal during these epochs.
#' 
#' 
#' ### Question 3: What assembly processes characterize temperature transitions
#+ eval=TRUE
# Temperature type
temptype.plot.prop.df <- betanull.lf %>%
  select(Site1, Site2, Assembly_Process, TempConditionType) %>%
  group_by(TempConditionType, Assembly_Process) %>%
  tally() %>% 
  ungroup() %>% group_by(TempConditionType) %>%
  mutate(Total = sum(n),
         Percent = 100*n/Total)

temptype.plot.prop <- temptype.plot.prop.df %>%
  ggplot(aes(x= TempConditionType, y = Percent, 
             fill = Assembly_Process)) +
  facet_wrap(~TempConditionType, shrink = TRUE, drop = TRUE, scales = "free_x") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Assembly Process", values = fill_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  xlab("") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  #coord_cartesian(expand = c(0,0)) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = rel(1)),
        strip.placement = "outside",
        legend.position = "bottom")

temptype.plot.prop.legend <- get_legend(temptype.plot.prop)
#type.plot.prop <- type.plot.prop + theme(legend.position = "none")

ggsave(temptype.plot.prop, 
       filename = paste0(figures.fp, "/tempType_prop.png"),
       width = 10, height = 10, dpi = 400)
temptype.plot.prop
#' I plotted the proportion of different types of assembly processes in Warm-Cold transitions as well as pairwise comparisons where both samples were considered either warm or cold
#' - I expected that temperature changes would have a stronger homogeneous selective effect than those between the same temperature, but that was not the case.
#' - In fact, Cold-Cold comparisons seemed to have the strongest selective pressure, while warm-warm and warm-cold had about the same amount of homogeneous (abiotic) selection. 
#' - This may mean that cold exerts a greater selective pressure than warm, or possibly that cold conditions lead to other selective processes that are not present when warm conditions are present.
#' - The cold-warm transition is one of the few places we see heterogeneous selection which could be consistent with the idea that once "woken up" competition forces play a larger role in some communities.
#' - The major difference between same-same comparisons and "different" comparisions is in the dispersal processes. These appear to play a greater role, but both high dispersal, and dispersal limitation are present. 
#'
#' ### Question 4: Do those temperature processes differ in different Climate Epochs?
#' Given that we know that drift increases with depth, it might be a good idea to understand how cold-warm transitions play out within a particular climate Epoch, in case this is skewing our results
#+ eval=TRUE
# Temperature comparisons within epoch type
epochtemptype.plot.prop.df <- betanull.lf %>%
  select(Site1, Site2, Assembly_Process, EpochType, TempConditionType) %>%
  group_by(EpochType, TempConditionType, Assembly_Process) %>%
  tally() %>% 
  ungroup() %>% group_by(EpochType, TempConditionType) %>%
  mutate(Total = sum(n),
         Percent = 100*n/Total)

epochtemptype.plot.prop <- epochtemptype.plot.prop.df %>%
  ggplot(aes(x= EpochType, y = Percent, 
             fill = Assembly_Process)) +
  facet_grid(TempConditionType~EpochType, drop = TRUE, scales = "free_x") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Assembly Process", values = fill_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  xlab("") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  #coord_cartesian(expand = c(0,0)) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(1.5)),
        panel.spacing.y =  unit(1, "lines"),
        panel.spacing.x =  unit(0, "lines"),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = rel(1)),
        strip.placement = "outside",
        legend.position = "bottom")


epochtemptype.plot.prop.legend <- get_legend(epochtemptype.plot.prop)
#type.plot.prop <- type.plot.prop + theme(legend.position = "none")

ggsave(epochtemptype.plot.prop, 
       filename = paste0(figures.fp, "/epochtempType_prop.png"),
       width = 10, height = 11, dpi = 400)
epochtemptype.plot.prop

#' - Here we see that there seems to be a greater increase in stochastic processes, as we saw before
#' - But one important difference is that the presence of dispersal processes does still seem to correspond to warm-cold transitions in the oldest climate epochs
#' - There are not clear trends within a given climate epoch of the differences in a cold-cold, cold-warm, or warm-warm trnaisition. Although Warm-warm and cold-warm transitions do seem to have more similarity than cold-cold. 
#' - With the exception of the Pre-LGS:LGS transition, cold-cold transitions are very strongly shaped by homogenizing selection, much more so than the warm-warm and cold-warm transitions in the same epoch.
#'
#'

#' ## How do results change when we restrict to only consecutive comparisons?
#' 
#+ eval=TRUE
# Proportions across all samples
betanull_consecutive.lf %>%
  select(Site1, Site2, BetaNTI, RCBC, Assembly_Process) %>%
  group_by(Assembly_Process) %>%
  tally() %>%
  mutate(Total = sum(n),
         Percent = round(100*n/Total, digits = 2)) %>%
  arrange(desc(Percent)) %>% knitr::kable()

#' - largely the same
#' ### How do the assembly processes change through time when only consecutive times are considered?

epochtemptype.plot.prop.df <- betanull_consecutive.lf %>%
  select(Site1, Site2, Assembly_Process, EpochType, TempConditionType, AvgDepth, AvgDepth.Site1) %>%
  mutate(DepthRange = paste0(round(AvgDepth, 1), "-", round(AvgDepth.Site1, 1)),
         DepthRange = fct_reorder(DepthRange, AvgDepth)) %>%
  group_by(DepthRange, EpochType, TempConditionType, Assembly_Process) %>%
  tally() %>% 
  ungroup() %>%
  mutate(DescreteDepth = fct_rev(DepthRange))

temperatureinfo <- epochtemptype.plot.prop.df %>% 
  select(DescreteDepth, EpochType, TempConditionType, n) %>%
  distinct()

#+ eval =TRUE
epochsummary.prop.plot.consec <- ggplot(epochtemptype.plot.prop.df, aes(x = EpochType)) +
  geom_bar( aes(fill = Assembly_Process), position = "fill") +
  guides(fill = guide_legend(nrow = 2)) + 
  scale_fill_manual(name = "Assembly Process", values = fill_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_x_discrete(expand = c(0,0)) +
  #coord_cartesian(expand = c(0,0)) +
  ylab("") +
  theme_bw() + 
  theme(axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title.x = element_blank(),
        panel.spacing.y =  unit(1, "lines"),
        panel.spacing.x =  unit(0, "lines"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = rel(1)),
        strip.placement = "outside",
        legend.position = "bottom")

ggsave(epochsummary.prop.plot.consec, 
       filename = paste0(figures.fp, "/consecutive_epoch_prop.png"),
       width = 5, height = 12, dpi = 400)
epochsummary.prop.plot.consec

#' - Now that we have restricted to only consecutive comparisons, the increasing influence of drift through time has become much reduced, and less linear. This is perhaps further evidence that the that signal was due to comparing communities that were very old in time to those that were very new. 
#' - Keep in mind that although these are reported as percentages, the total number of comparisons is very different among time periods. The Epoch transitions now represent only 1 sample (see below for an unscaled plot).  
#+ eval =TRUE, echo=FALSE
epochtemptype.plot.prop.df %>% 
  group_by(EpochType) %>% tally() %>%
  rename(`Number of Comparisons` = n) %>%
  knitr::kable()
#' - What we see here is that selection is still prominent in the dataset but that the Holocene is particularly marked by selection, while the LGS and pre-LGS are much more variable.
#+ eval=TRUE,
epochsummary.prop.plot.consec  +
  geom_bar( aes(fill = Assembly_Process)) +
  scale_y_continuous(expand = c(0,0)) + ylab("Count")


#' - If we plot the same data along the depth of the core with the temperature transitions, we can see that during the LGS these sometimes, but not always correspond to temperature changes. 
#' 
#+ eval =TRUE, fig.dim = c(6, 12)
epochtemptype.plot.consec <- ggplot(epochtemptype.plot.prop.df, aes(x = DescreteDepth, y = n)) +
  geom_bar(stat = "identity", aes(fill = Assembly_Process),
           linewidth = 5, # temperature border thickness
           width = 1) + # bar thickness
  geom_errorbar(data = temperatureinfo, 
            aes(y = n, ymax = n*0.02, ymin = n*0.02, color = TempConditionType),
                linewidth = 2, width = 1) +
  # geom_text(aes(label = TempConditionType, 
  #               y = n*0.2, color = TempConditionType),
  #           vjust = 1) +
  coord_flip() +
  facet_grid(EpochType~., drop = TRUE, scales = "free", space = "free_y") +
  scale_fill_manual(name = "Assembly Process", values = fill_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  guides(fill = guide_legend(nrow = 5)) + 
  scale_color_manual(name = "Temperature Transition", values = c("blue", "red", "black"),
                    breaks = c("Cold-Cold", "Warm-Warm", "Warm-Cold"),
                    labels = c("Cold-Cold", "Warm-Warm", "Warm-Cold")) +
  ylab("") + 
  xlab("Community Transitions (by increasing depth)") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() + 
  theme(axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(1)),
        axis.text.x = element_blank(),
        panel.spacing.y =  unit(0, "lines"),
        panel.spacing.x =  unit(0, "lines"),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        legend.title.position = "top",
        legend.position = "right")

ggsave(epochtemptype.plot.consec, 
       filename = paste0(figures.fp, "/consecutive_epoch.png"),
       width = 5, height = 12, dpi = 400)
epochtemptype.plot.consec

#' ### How do the assembly processes with temperature when only consecutive times are considered?

#+ eval = TRUE
temptype.plot.prop.consec <- ggplot(epochtemptype.plot.prop.df, aes(x = TempConditionType)) +
  geom_bar( aes(fill = Assembly_Process), position = "fill") +
  guides(fill = guide_legend(nrow = 2)) + 
  scale_fill_manual(name = "Assembly Process", values = fill_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_x_discrete(expand = c(0,0)) +
  #coord_cartesian(expand = c(0,0)) +
  ylab("") +
  theme_bw() + 
  theme(axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title.x = element_blank(),
        panel.spacing.y =  unit(1, "lines"),
        panel.spacing.x =  unit(0, "lines"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = rel(1)),
        strip.placement = "outside",
        legend.position = "bottom")

ggsave(temptype.plot.prop.consec, 
       filename = paste0(figures.fp, "/consecutive_temptype.png"),
       width = 4, height = 12, dpi = 400)
temptype.plot.prop.consec
#' 
#' 
#' 
temptype.epoch.plot.prop.consec <- temptype.plot.prop.consec + 
  facet_grid(~EpochType, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 10),
        panel.spacing.x =  unit(0.5, "lines"))
ggsave(temptype.epoch.plot.prop.consec, 
       filename = paste0(figures.fp, "/consecutive_temptype_splitbyepoch.png"),
       width = 4, height = 12, dpi = 400)
temptype.epoch.plot.prop.consec
#' - *Cold-cold*: We see that selection is more prominent in cold-cold transitions generally, although in the Holocene, nearly all transitions are dominated by selection regardless of temperature transition
#' - *Warm-Warm* and *Warm-Cold*: In the epoch's prior to the Holocene there are greater stochastic processes when a warm transition is involved (either warm-warm or warm-cold transition). More dispersal, and drift dominate. However there is still an appreciable selective signal even in these transitions, although the selection is not consistent between periods (Warm-Warm has selection in the pre-LGS time, Warm-Cold has selection in the LGS). 
#' 
#' 
#' 
#' ### Bonus: Are any other factors related to assembly processes?
#' #### Dust differences and Assembly process
#+ eval=TRUE
#### ====================================================================== ####
DustDiff.betanull <- betanull.lf %>%
  left_join(input_all$sample_metadata %>% select(SampleID, `Dust count per ml ice (diameter >0.63 μm)`), 
            by = c("Site1" = "SampleID")) %>%
  rename(Dust1 = `Dust count per ml ice (diameter >0.63 μm)`) %>%
  left_join(input_all$sample_metadata %>% select(SampleID, `Dust count per ml ice (diameter >0.63 μm)`), 
            by = c("Site2" = "SampleID")) %>%
  rename(Dust2 = `Dust count per ml ice (diameter >0.63 μm)`) %>%
  filter(Dust1 < 60, Dust2 < 60) %>% # filter out the really dusty sample
  mutate(DustDiff = abs(Dust1-Dust2)) %>%
  pivot_longer(cols = all_of(c("BetaNTI", "RCBC")), 
               names_to = "AssemblyMetric", values_to = "AssemblyValue")

corr.dustdiff.bnti <- DustDiff.betanull %>%
  filter(AssemblyMetric == "BetaNTI") %>%
  select(AssemblyValue, DustDiff) %>%
  cor.test(~ AssemblyValue + DustDiff, data = ., method = "pearson")

corr.dustdiff.RCBC <- DustDiff.betanull %>%
  filter(AssemblyMetric == "RCBC") %>%
  select(AssemblyValue, DustDiff) %>%
  cor.test(~ AssemblyValue + DustDiff, data = ., method = "pearson")


ggplot(DustDiff.betanull, 
       aes(x = DustDiff, y = AssemblyValue)) +
  geom_point(aes(color = Assembly_Process)) +
  scale_color_manual(name = "Assembly Process", values = fill_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels) +
  facet_wrap(~AssemblyMetric, scales = "free_y") +
  theme_bw()

#' As dust amounds increase in difference between samples, dispersal processes 
#' do not appear to increase
#' 
#' 

#' #### Ice Volume differences and Assembly process
#+ eval=TRUE
IceDiff.betanull <- betanull.lf %>%
  left_join(input_all$sample_metadata %>% select(SampleID, `Ice volume (ml)`), 
            by = c("Site1" = "SampleID")) %>%
  rename(Ice1 = `Ice volume (ml)`) %>%
  left_join(input_all$sample_metadata %>% select(SampleID, `Ice volume (ml)`), 
            by = c("Site2" = "SampleID")) %>%
  rename(Ice2 = `Ice volume (ml)`) %>%
  mutate(IceDiff = abs(Ice1-Ice2)) %>%
  pivot_longer(cols = all_of(c("BetaNTI", "RCBC")), 
               names_to = "AssemblyMetric", values_to = "AssemblyValue")

corr.icediff.bnti <- IceDiff.betanull %>%
  filter(AssemblyMetric == "BetaNTI") %>%
  select(AssemblyValue, IceDiff) %>%
  cor.test(~ AssemblyValue + IceDiff, data = ., method = "pearson")

corr.icediff.RCBC <- IceDiff.betanull %>%
  filter(AssemblyMetric == "RCBC") %>%
  select(AssemblyValue, IceDiff) %>%
  cor.test(~ AssemblyValue + IceDiff, data = ., method = "pearson")

ggplot(IceDiff.betanull, 
       aes(x = IceDiff, y = AssemblyValue)) +
  geom_point(aes(color = Assembly_Process)) +
  scale_color_manual(name = "Assembly Process", values = fill_assembly, 
                     breaks = assembly_levels, 
                     labels = assembly_labels) +
  facet_wrap(~AssemblyMetric, scales = "free_y") +
  theme_bw()

#' #### Cryobacterium differences and Assembly process
#' As ice amounds increase in difference between samples, dispersal processes 
#' do not appear to increase
#+ eval=TRUE
# Cryobacterium names
cryobacterium <- input_all$taxonomy %>% filter(Genus == "D_5__Cryobacterium")
Cryo_abund <- input_all$otu_table %>%
  mutate(Cryobacterium = ifelse(OTU_ID %in% cryobacterium$OTU_ID, "Cryobacterium", "Other")) %>%
  pivot_longer(starts_with("D"), names_to = "SampleID", values_to = "Counts") %>%
  group_by(SampleID) %>%
  mutate(TotalCount = sum(Counts)) %>% ungroup() %>%
  group_by(SampleID, Cryobacterium) %>%
  mutate(CryoCount = sum(Counts)) %>%
  ungroup() %>% group_by(SampleID) %>%
  mutate(RelAbundCryo =  CryoCount/TotalCount) %>%
  select(SampleID, Cryobacterium, RelAbundCryo, TotalCount, CryoCount) %>%
  filter(Cryobacterium!="Other") %>%
  distinct()


#' Does the relative abundance of Cryobacterium appear to influence assemlby processes?
#+ eval=TRUE
CryoDiff.betanull <- betanull.lf %>%
  left_join(Cryo_abund %>% select(SampleID, RelAbundCryo), 
            by = c("Site1" = "SampleID")) %>%
  rename(Cryo1 = RelAbundCryo) %>%
  left_join(Cryo_abund %>% select(SampleID, RelAbundCryo), 
            by = c("Site2" = "SampleID")) %>%
  rename(Cryo2 = RelAbundCryo) %>%
  mutate(CryoDiff = abs(Cryo1-Cryo2)) %>%
  pivot_longer(cols = all_of(c("BetaNTI", "RCBC")), 
               names_to = "AssemblyMetric", values_to = "AssemblyValue")

corr.cryodiff.bnti <- CryoDiff.betanull %>%
  filter(AssemblyMetric == "BetaNTI") %>%
  select(AssemblyValue, CryoDiff) %>%
  cor.test(~ AssemblyValue + CryoDiff, data = ., method = "pearson")

corr.cryodiff.RCBC <- CryoDiff.betanull %>%
  filter(AssemblyMetric == "RCBC") %>%
  select(AssemblyValue, CryoDiff) %>%
  cor.test(~ AssemblyValue + CryoDiff, data = ., method = "pearson")



ggplot(CryoDiff.betanull, 
       aes(x = CryoDiff, y = AssemblyValue)) +
  geom_point(aes(color = Assembly_Process)) +
  scale_color_manual(name = "Assembly Process", values = fill_assembly, 
                     breaks = assembly_levels, 
                     labels = assembly_labels) +
  facet_wrap(~AssemblyMetric, scales = "free_y") +
  theme_bw()

#' - The relative abundance of Cryobacterium does appear to influence assembly processes (pearson correlation ~ 0.4 in each case, significant at the 0.05 level). 
#' - Samples with more similar Cryobacterium communities have stronger homogenizing selective processes than those with very different cyobacterium relative abundances. Similarly, for stochastic processes, very different cryobacterium relative abundances seem related to dispersal limitation, while similar abundances have more homogenizing selection. 
#' 
#' 
#' ## Concluding Thoughts
#' - We see strong evidence of abiotic selective processes driving the structure of these communities. Ecological drift is a second dominant assembly process. There is little evidence of biotic structuring within the community.
#' - This doesn't mean that there is no biotic structuring, just that it doesn't show up as a sister-taxa competition. Other forms of biotic structuring, such as food-web type structuring cannot be ruled out via this method. What this does imply is that the conditions at hand, strongly favor particular sets of clades. 
#' - Ecological drift increases with time, although when the analysis is restricted to consecutive comparisons only, this increase is less linear. This is hard to attribute. We cannot rule out the influence of relic DNA, or active entrained cells in this pattern.  However one piece of evidence against active entrained cells is a lack of "dispersal limitation with drift", which I would expect to see more of if cells were actively growing in these ice cores. What is clear is that there are stronger abiotic selective pressures nearer the top of the core.
#' - Climate epoch has a stronger influence on assembly process than temperature variation within an epoch. This may be due to different environmental/ecological forces being at play during each epoch. For example, climatic shifts, deposition shifts, or different sources of microbial communities. 
#' - In particular, the Holocene is a very selection-dominated time, despite the capture of multiple and Warm-Warm, Warm-Cold transitions. In the LGS, by contrast, Warm-Warm and Warm-Cold transitions are dominated by more stochastic processes, drift, and high dispersal. A similar pattern is seen in the Pre-LGS period, but the Warm-cold/Warm-Warm transitions are not consistent between them in amount of selection at play. 
#' - Microbial communities are significantly different in warm vs. cold periods. Interestingly cold seems to impart a greater selective pressure than warmth (though this is really only apparent prior to the Holocene). This may align with the role of photorophs in structuring the community, since colder years likely include less available light for phototrophs which could have downstream selective pressures on the community. 

