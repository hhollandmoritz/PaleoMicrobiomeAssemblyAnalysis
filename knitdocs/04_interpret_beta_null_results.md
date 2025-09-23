---
title: Interpret Assembly Analysis results for Paleo Microbiome Project
date: '2025-06-11'
date-modified: '`r Sys.Date()`'
author: Hannah Holland-Moritz
format:
  html:
    toc: true
    toc-depth: 4
    page-layout: full
    embed-resources: true
    self-contained-math: true
    html-math-method: katex
    theme: sandstone
    code-fold: true
    code-summary: Show the code
  pdf:
    toc: true
    number-sections: true
    colorlinks: true
  md:
    variant: gfm+yaml_metadata_block+definition_lists
    prefer-html: true
    fig-format: retina
    fig-width: 8
    fig-height: 5
    wrap: preserve
  docx:
    toc: true
    number-sections: true
    highlight-style: github
    prefer-html: true
editor: visual
editor_options:
  chunk_output_type: inline
---


<script src="04_interpret_beta_null_results_files/libs/kePrint-0.0.1/kePrint.js"></script>
<link href="04_interpret_beta_null_results_files/libs/lightable-0.0.1/lightable.css" rel="stylesheet" />


***This step takes the output from 03_prepare_and_combine_assembly_analysis_results.R and runs initial figures and analyses on them. It has been knit into a notebook for ease of interpretation***

# Driving questions:

1.  What assembly processes predominate in these paleo-microbiome samples?
2.  What assembly processes characterize each Epoch, and the Epoch transitions?
3.  What assembly processes drive cold to warm transitions?
4.  Do those processes differ in different Climate Epochs? Bonus: Are any other factors related to assembly processes?

#### How does assembly analysis work? Find out more at this link: [SlideShow](https://docs.google.com/presentation/d/1BSLtMNZZrXxR9Nk0GER2RBBGPU0I1v4zAgnr9jylAQw/edit?usp=sharing)

## Some important data processing notes:

-   Comparisons between Pre-LGS and Holocene samples have been filtered from the dataset. The logic behind this is that it doesn’t make much sense to compare samples that are not consecutive in time.

-   Epochs are ordered by time

-   Temperature comparisons are ordered by “same-same”, “different”

-   After filtering out the Pre-LGS:Holocene comparisons, we are left with 468 comparisons

### Question 1: What assembly processes predominate in these paleo-microbiome samples?

``` r
# Calculate proportions of pairwise comparisons overall
#### ====================================================================== ####

# extra data from other papers:
# from: https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.15651
#other_data <- data.frame(Assembly_Process = c("Homogenous selection", "Heterogenous selection", "homogenizing dispersal", "Drift"))

# another good paper: 
# https://academic.oup.com/ismej/article/16/12/2653/7474068


# arctic / antarctic ocean sediments https://doi.org/10.1016/j.marenvres.2025.107261



# Proportions across all samples
betanull.lf %>%
  select(Site1, Site2, BetaNTI, RCBC, Assembly_Process) %>%
  group_by(Assembly_Process) %>%
  tally() %>%
  mutate(Total = sum(n),
         Percent = round(100*n/Total, digits = 2)) %>%
  arrange(desc(Percent)) %>% 
  kbl() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| Assembly_Process               |   n | Total | Percent |
|:-------------------------------|----:|------:|--------:|
| Homogenous selection           | 282 |   468 |   60.26 |
| Drift                          | 150 |   468 |   32.05 |
| Dispersal limitation and drift |  24 |   468 |    5.13 |
| Homogenizing dispersal         |  10 |   468 |    2.14 |
| Heterogenous selection         |   2 |   468 |    0.43 |

-   Most of the pairwise comparisons are characterized by “Homogenous selection”. This is typically considered to be the result of abiotic selection that favors some clades over others, rather than filtering sister taxa.
-   The second most common assembly process is ecological drift. This is a purely stochastic process. Assembly cannot be attributed to selection or dispersal.
-   We see essentially no evidence (only 2/468 comparisons) of heterogeneous (biotic) selection

### Question 2: What assembly processes characterize each Epoch, and the Epoch transitions?

``` r
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

#epochtype.plot.prop.legend <- get_legend(epochtype.plot.prop)
#type.plot.prop <- type.plot.prop + theme(legend.position = "none")

ggsave(epochtype.plot.prop, 
       filename = paste0(figures.fp, "/epochType_prop.png"),
       width = 10, height = 10, dpi = 400)
epochtype.plot.prop
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-6-1.png" style="width:98.0%;height:98.0%" />

-   As time proceeds, we see an increase in the proportion of pairwise comparisons characterized by drift in our dataset.
-   This pattern is often observed with depth in soils, more generally. It could be due to relic DNA (which experiences no selection), or it could be indicative that in previous epochs there were less selective pressures on the organisms forming these communities.
-   Interestingly, we don’t see major increases in selective processes at the Epoch-Epoch transistions.
-   But, dispsersal limitation seems to have played a fairly strong role in the Pre-LGS:LGS transition.
-   The modern era communities are not at all characterized by dispersal, but there appears to be more of it in the LGS, and Pre-LGS periods.
-   This may suggest greater wind or water-mediated dispersal during these epochs.

### Question 3: What assembly processes characterize temperature transitions

``` r
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
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-7-1.png" style="width:98.0%;height:98.0%" />

I plotted the proportion of different types of assembly processes in Warm-Cold transitions as well as pairwise comparisons where both samples were considered either warm or cold - I expected that temperature changes would have a stronger homogeneous selective effect than those between the same temperature, but that was not the case. - In fact, Cold-Cold comparisons seemed to have the strongest selective pressure, while warm-warm and warm-cold had about the same amount of homogeneous (abiotic) selection. - This may mean that cold exerts a greater selective pressure than warm, or possibly that cold conditions lead to other selective processes that are not present when warm conditions are present. - The cold-warm transition is one of the few places we see heterogeneous selection which could be consistent with the idea that once “woken up” competition forces play a larger role in some communities. - The major difference between same-same comparisons and “different” comparisions is in the dispersal processes. These appear to play a greater role, but both high dispersal, and dispersal limitation are present.

### Question 4: Do those temperature processes differ in different Climate Epochs?

Given that we know that drift increases with depth, it might be a good idea to understand how cold-warm transitions play out within a particular climate Epoch, in case this is skewing our results

``` r
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
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-8-1.png" style="width:98.0%;height:98.0%" />

-   Here we see that there seems to be a greater increase in stochastic processes, as we saw before
-   But one important difference is that the presence of dispersal processes does still seem to correspond to warm-cold transitions in the oldest climate epochs
-   There are not clear trends within a given climate epoch of the differences in a cold-cold, cold-warm, or warm-warm transition. Although Warm-warm and cold-warm transitions do seem to have more similarity than cold-cold.
-   With the exception of the Pre-LGS:LGS transition, cold-cold transitions are very strongly shaped by homogenizing selection, much more so than the warm-warm and cold-warm transitions in the same epoch.

## How do results change when we restrict to only consecutive comparisons?

Because of the way the ice core forms, we might logically decide to restrict comparisons to only consecutive samples. These results show the answers for those comparisions.

### Question 1: What assembly processes predominate in these paleo-microbiome samples?

``` r
# Proportions across all samples
betanull_consecutive.lf %>%
  select(Site1, Site2, BetaNTI, RCBC, Assembly_Process) %>%
  group_by(Assembly_Process) %>%
  tally() %>%
  mutate(Total = sum(n),
         Percent = round(100*n/Total, digits = 2)) %>%
  arrange(desc(Percent)) %>% kbl(digits = 3) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| Assembly_Process               |   n | Total | Percent |
|:-------------------------------|----:|------:|--------:|
| Homogenous selection           |  19 |    32 |   59.38 |
| Drift                          |   9 |    32 |   28.12 |
| Homogenizing dispersal         |   3 |    32 |    9.38 |
| Dispersal limitation and drift |   1 |    32 |    3.12 |

-   The proportions are largely the same as when we had non-consecutive comparisons

### Question 2: How do the assembly processes change through time when only consecutive times are considered?

``` r
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
```

``` r
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
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-1.png" style="width:98.0%;height:98.0%" />

-   Now that we have restricted to only consecutive comparisons, the increasing influence of drift through time has become much reduced, and less linear. This is perhaps further evidence that the that signal was due to comparing communities that were very old in time to those that were very new.
-   Keep in mind that although these are reported as percentages, the total number of comparisons is very different among time periods. The Epoch transitions now represent only 1 sample (see below for an unscaled plot).

| EpochType    | Number of Comparisons |
|:-------------|----------------------:|
| Holocene     |                     9 |
| LGS          |                    16 |
| LGS:Holocene |                     1 |
| Pre-LGS      |                     5 |
| Pre-LGS:LGS  |                     1 |

-   What we see here is that selection is still prominent in the dataset but that the Holocene is particularly marked by selection, while the LGS and pre-LGS are much more variable.

``` r
epochsummary.prop.plot.consec  +
  geom_bar( aes(fill = Assembly_Process)) +
  scale_y_continuous(expand = c(0,0)) + ylab("Count")
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-13-1.png" style="width:98.0%;height:98.0%" />

-   If we plot the same data along the depth of the core with the temperature transitions, we can see that during the LGS these sometimes, but not always correspond to temperature changes.

``` r
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
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-14-1.png" style="width:98.0%;height:98.0%" />

### Question 3: How do the assembly processes with temperature when only consecutive times are considered?

``` r
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
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-15-1.png" style="width:98.0%;height:98.0%" />

### Question 4: Do these consecutive temperature transitions differ by different epochs?

``` r
temptype.epoch.plot.prop.consec <- temptype.plot.prop.consec + 
  facet_grid(~EpochType, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 10),
        panel.spacing.x =  unit(0.5, "lines"))
ggsave(temptype.epoch.plot.prop.consec, 
       filename = paste0(figures.fp, "/consecutive_temptype_splitbyepoch.png"),
       width = 4, height = 12, dpi = 400)
temptype.epoch.plot.prop.consec
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-16-1.png" style="width:98.0%;height:98.0%" />

-   *Cold-cold*: We see that selection is more prominent in cold-cold transitions generally, although in the Holocene, nearly all transitions are dominated by selection regardless of temperature transition
-   *Warm-Warm* and *Warm-Cold*: In the epoch’s prior to the Holocene there are greater stochastic processes when a warm transition is involved (either warm-warm or warm-cold transition). More dispersal, and drift dominate. However there is still an appreciable selective signal even in these transitions, although the selection is not consistent between periods (Warm-Warm has selection in the pre-LGS time, Warm-Cold has selection in the LGS).

## Bonus: Are any other factors related to assembly processes?

(Analyses below have been restricted to consecutive sample comparisions).

### Dust differences and Assembly process

``` r
#### ====================================================================== ####
dust_data <- input_all$sample_metadata %>% select(SampleID, `Dust count per ml ice (diameter >0.63 μm)`)
compute_dust_difference <- function(betanull = betanull_consecutive.lf, 
                                    sample_dust_data = dust_data,
                                    dust_column = "Dust count per ml ice (diameter >0.63 μm)") {
 DustDiff <- betanull %>%
  left_join(sample_dust_data, 
            by = c("Site1" = "SampleID")) %>%
  rename(Dust1 = all_of(dust_column)) %>%
  left_join(sample_dust_data, 
            by = c("Site2" = "SampleID")) %>%
  rename(Dust2 = all_of(dust_column)) %>%
  mutate(DustDiff = abs(Dust1-Dust2)) %>%
  pivot_longer(cols = all_of(c("BetaNTI", "RCBC")), 
               names_to = "AssemblyMetric", values_to = "AssemblyValue")
 
 return(DustDiff)
}

DustDiff.betanull <- compute_dust_difference() %>%
  # filter out comporisons that include that crazy high dust sample
  filter(DustDiff < 60)

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
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-17-1.png" style="width:98.0%;height:98.0%" />

As dust greater than 0.63 micrometers increases in difference between samples, dispersal processes do not appear to increase

What about other sizes of dust?

``` r
dust_data <- dust_size %>% select(Sample_name, contains("dust_")) %>%
  rename(SampleID = Sample_name) %>%
  pivot_longer(-SampleID, names_to = "dust_size_class", values_to = "dust_amount") %>%
  group_by(dust_size_class)

ggplot(dust_data, aes(x = dust_amount)) +
  geom_histogram() +
  facet_wrap(~dust_size_class, scale = "free")
```

<figure>
<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-18-1.png" style="width:98.0%;height:98.0%" alt="Here is the distribution of dust sizes." />
<figcaption aria-hidden="true">Here is the distribution of dust sizes.</figcaption>
</figure>

``` r
dust_data_betanull <- dust_data %>% nest() %>%
  # prepare data for plotting and testing
  mutate(dust_data_betanull = purrr::map(data, ~compute_dust_difference(betanull = betanull_consecutive.lf, # use consecutive comparisons only
                                                                     sample_dust_data = .x,
                                                                     dust_column = "dust_amount"))) %>%
  # compute correlation of betaNTI
  mutate(betanti_corr = purrr::map(dust_data_betanull, ~filter(.x, AssemblyMetric == "BetaNTI") %>%
                                     select(AssemblyValue, DustDiff) %>%
                                     cor.test(~ AssemblyValue + DustDiff, data = ., method = "pearson"))) %>%
  # compute correlation of RCBC
  mutate(RCBC_corr = purrr::map(dust_data_betanull, ~filter(.x, AssemblyMetric == "RCBC") %>%
                                     select(AssemblyValue, DustDiff) %>%
                                     cor.test(~ AssemblyValue + DustDiff, data = ., method = "pearson"))) %>%
  # Plot
  mutate(dust_size_plot = purrr::map2(dust_data_betanull, dust_size_class, ~ggplot(.x, aes(x = DustDiff, y = AssemblyValue)) +
                                       geom_point(aes(color = Assembly_Process)) +
                                       scale_color_manual(name = "Assembly Process", values = fill_assembly,
                                                          breaks = assembly_levels,
                                                          labels = assembly_labels) +
                                       facet_wrap(~AssemblyMetric, scales = "free_y") +
                                       theme_bw() + ggtitle(.y) + theme(legend.position = "none")))
```

``` r
# now pull out the correlations

#response_plots <- purrr::map(flatten(dust_data_betanull$dust_size_plot), ~cowplot::plot_grid(plotlist = .x))

cowplot::plot_grid(plotlist = dust_data_betanull$dust_size_plot, ncol = 2)
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-20-1.png" style="width:98.0%;height:98.0%" />

``` r
# now pull out the correlations
dust_levels <- dust_data_betanull$dust_size_class
dust_betanti_corr <- dust_data_betanull %>% 
  mutate(betanti_tidy = purrr::map(betanti_corr, ~broom::tidy(.x))) %>%
  unnest(betanti_tidy) %>%
  select(dust_size_class, estimate:conf.high) %>%
  mutate(dust_size_class = factor(dust_size_class, levels = dust_levels))
  
dust_betanti_corr %>% 
  kbl(digits = 3) %>%
  column_spec(4, bold = ifelse((dust_betanti_corr$p.value > 0.05 | is.na(dust_betanti_corr$p.value)), FALSE, TRUE)) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| dust_size_class   | estimate | statistic | p.value | parameter | conf.low | conf.high |
|:------------------|---------:|----------:|--------:|----------:|---------:|----------:|
| dust_0.63_0.8µm   |   -0.049 |    -0.230 |   0.820 |        22 |   -0.444 |     0.362 |
| dust_0.8_1µm      |   -0.046 |    -0.218 |   0.829 |        22 |   -0.442 |     0.364 |
| dust_1_1.26µm     |   -0.043 |    -0.200 |   0.843 |        22 |   -0.438 |     0.367 |
| dust_1.26_1.59µm  |   -0.026 |    -0.121 |   0.905 |        22 |   -0.425 |     0.382 |
| dust_1.59_2.02µm  |    0.002 |     0.008 |   0.994 |        22 |   -0.402 |     0.405 |
| dust_2.02_2.52µm  |    0.044 |     0.209 |   0.837 |        22 |   -0.365 |     0.440 |
| dust_2.52_3.14µm  |    0.119 |     0.563 |   0.579 |        22 |   -0.299 |     0.499 |
| dust_3.14_4µm     |    0.187 |     0.891 |   0.382 |        22 |   -0.234 |     0.549 |
| dust_4_5.04µm     |    0.240 |     1.160 |   0.259 |        22 |   -0.181 |     0.587 |
| dust_5.04_6.35µm  |    0.239 |     1.157 |   0.260 |        22 |   -0.181 |     0.586 |
| dust_6.35_8µm     |    0.232 |     1.121 |   0.274 |        22 |   -0.189 |     0.581 |
| dust_8_10.08µm    |    0.236 |     1.141 |   0.266 |        22 |   -0.185 |     0.584 |
| dust_10.08_12.7µm |    0.293 |     1.437 |   0.165 |        22 |   -0.125 |     0.623 |
| dust_12.7_16.0    |    0.335 |     1.667 |   0.110 |        22 |   -0.079 |     0.650 |
| dust_total        |   -0.033 |    -0.154 |   0.879 |        22 |   -0.431 |     0.376 |

Difference in dust amount for different fractions are *not* *significantly* correlated to deterministic assembly processes. (Note: when this is run on non-consecutive samples, there are significant relationships for fractions \<4 µm, but I have not shown that analysis here)

``` r
# now pull out the correlations

dust_RCBC_corr <- dust_data_betanull %>% 
  mutate(RCBC_tidy = purrr::map(RCBC_corr, ~broom::tidy(.x))) %>%
  unnest(RCBC_tidy) %>%
  select(dust_size_class, estimate:conf.high) %>%
  mutate(dust_size_class = factor(dust_size_class, levels = dust_levels))
  
dust_RCBC_corr %>% 
  kbl(digits = 3) %>%
  column_spec(4, bold = ifelse((dust_RCBC_corr$p.value > 0.05 | is.na(dust_RCBC_corr$p.value)), FALSE, TRUE)) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| dust_size_class   | estimate | statistic | p.value | parameter | conf.low | conf.high |
|:------------------|---------:|----------:|--------:|----------:|---------:|----------:|
| dust_0.63_0.8µm   |    0.137 |     0.367 |   0.724 |         7 |   -0.580 |     0.735 |
| dust_0.8_1µm      |    0.246 |     0.670 |   0.524 |         7 |   -0.500 |     0.782 |
| dust_1_1.26µm     |    0.244 |     0.666 |   0.527 |         7 |   -0.501 |     0.782 |
| dust_1.26_1.59µm  |    0.366 |     1.042 |   0.332 |         7 |   -0.394 |     0.829 |
| dust_1.59_2.02µm  |    0.429 |     1.258 |   0.249 |         7 |   -0.328 |     0.851 |
| dust_2.02_2.52µm  |    0.471 |     1.414 |   0.200 |         7 |   -0.280 |     0.865 |
| dust_2.52_3.14µm  |    0.501 |     1.530 |   0.170 |         7 |   -0.245 |     0.874 |
| dust_3.14_4µm     |    0.520 |     1.611 |   0.151 |         7 |   -0.220 |     0.880 |
| dust_4_5.04µm     |    0.521 |     1.614 |   0.151 |         7 |   -0.219 |     0.880 |
| dust_5.04_6.35µm  |    0.504 |     1.543 |   0.167 |         7 |   -0.241 |     0.875 |
| dust_6.35_8µm     |    0.500 |     1.529 |   0.170 |         7 |   -0.245 |     0.874 |
| dust_8_10.08µm    |    0.487 |     1.476 |   0.183 |         7 |   -0.261 |     0.870 |
| dust_10.08_12.7µm |    0.513 |     1.581 |   0.158 |         7 |   -0.229 |     0.878 |
| dust_12.7_16.0    |    0.553 |     1.756 |   0.122 |         7 |   -0.176 |     0.890 |
| dust_total        |    0.386 |     1.106 |   0.305 |         7 |   -0.374 |     0.836 |

There is not a significant correlation between dust amount and stochastic assembly processes for any faction. (Note: Difference in dust amount for fractions under 6.35 µm show significant (p\<0.05) correlations when non-consecutive samples are included).

### Ice Volume differences and Assembly process

``` r
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
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-25-1.png" style="width:98.0%;height:98.0%" />

As ice amounts increase in difference between samples, dispersal processes do not appear to increase

#### Does the relative abundance of Cryobacterium appear to influence assembly processes?

``` r
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
```

``` r
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
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-27-1.png" style="width:98.0%;height:98.0%" />

-   The relative abundance of Cryobacterium does appear to influence assembly processes (pearson correlation ~ 0.4 in each case, significant at the 0.05 level).
-   Samples with more similar Cryobacterium communities have stronger homogenizing selective processes than those with very different cyobacterium relative abundances. Similarly, for stochastic processes, very different cryobacterium relative abundances seem related to dispersal limitation, while similar abundances have more homogenizing selection.

## Do the assembly processes change when rare taxa are excluded?

I ran the assembly analysis on the table where rare lineages had been filtered out.

| Table name | Filtering description |
|------------------------------------|------------------------------------|
| All | No filtering applied |
| max_mt_1.0 | only taxa with maximum rel. abund. \> 1% in at least 1 sample |
| mean_mt_1.0 | only taxa with mean rel. abund. \> 1% in the 33 samples |
| max_mt_0.1 | only taxa with maximum rel. abund. \> 0.1% in at least 1 sample |
| mean_mt_0.1 | only taxa with mean rel. abund. \> 0.1% in the 33 samples |

Description of the different taxa tables used and their filtering

``` r

calculate_percent_assembly_process <- function(filter_name = "max_mt_1.0", betanull_long_form = betanull_consecutive.lf.max_mt_1.0) {
  dt <- betanull_long_form %>%
  select(Site1, Site2, BetaNTI, RCBC, Assembly_Process) %>%
  group_by(Assembly_Process) %>%
  tally() %>%
  mutate(Total = sum(n),
         Percent = round(100*n/Total, digits = 2)) %>%
  arrange(desc(Percent)) %>%
  mutate(filter_level = filter_name)
  
  return(dt)
}

percentage_table <- purrr::map(otu_filter_levels, ~calculate_percent_assembly_process(filter_name = .x, 
                                                                  betanull_long_form = betanull_list_consecutive[[.x]])) %>%
  purrr::reduce(bind_rows)

percentage_table %>%
  select(Assembly_Process, Percent, filter_level) %>%
  pivot_wider(names_from = "filter_level", values_from = "Percent") %>%
  kbl() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| Assembly_Process | all | mean_mt_0.1 | max_mt_0.1 | mean_mt_1.0 | max_mt_1.0 |
|:---|---:|---:|---:|---:|---:|
| Homogenous selection | 59.38 | NA | 41.94 | NA | 9.68 |
| Drift | 28.12 | 58.06 | 48.39 | 68.75 | 64.52 |
| Homogenizing dispersal | 9.38 | 16.13 | 6.45 | 9.38 | 16.13 |
| Dispersal limitation and drift | 3.12 | 19.35 | 3.23 | 15.62 | 9.68 |
| Heterogenous selection | NA | 6.45 | NA | 6.25 | NA |

``` r
percentage_table %>%
  select(Assembly_Process, n, filter_level) %>%
  pivot_wider(names_from = "filter_level", values_from = "n") %>%
  kbl() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| Assembly_Process | all | mean_mt_0.1 | max_mt_0.1 | mean_mt_1.0 | max_mt_1.0 |
|:---|---:|---:|---:|---:|---:|
| Homogenous selection | 19 | NA | 13 | NA | 3 |
| Drift | 9 | 18 | 15 | 22 | 20 |
| Homogenizing dispersal | 3 | 5 | 2 | 3 | 5 |
| Dispersal limitation and drift | 1 | 6 | 1 | 5 | 3 |
| Heterogenous selection | NA | 2 | NA | 2 | NA |

Selection becomes less prominent (or non-existent) in datasets with rare taxa filtered out, instead drift takes over as the most prominent process.

Additionally there are interesting differences between the “max” and “mean” filtering strategy. In the “max” filtering strategy, taxa have to meet a less stringent ubiquity bar than the “mean” filtering strategy. Taxa don’t have to be present in multiple samples at high abundance, just 1 of the 33. So taxa that are typically rare, but can take over occasionally may make the cut. We could imagine that something like this might happen in the case of either high migration (dust deposition) or high growth (a temporary bloom). We see more selection occurring in the max-filtered datasets than in the mean-filtered ones. Furthermore, the less-strict the filtering level (0.1% is less strict than 1.0%), the more homogenizing selection dominates. This is likely because we are observing the large selective pressures that follow a bloom, or large migration event in one or a few samples. In contrast, the “mean” filtering strategy, experiences no homogenizing selection, and primarily is dominated by drift. I suspect this might be because the mean better captures the “core” community of ice-dwelling organisms that do not experience as strong of a selective pressure year over year.

When rare taxa (the likely ones putatively being blown in are filtered out, the process that dominates is drift). ***This suggests that Karna’s theory has merit: That there is a strong biotic filter imposed on the rare taxa that is not imposed on the more abundant taxa, possibly due to the deposition process?.***

#### How do these processes change our understandings of temperature transitions?

``` r
plot_tempchange <- function(betanull_long_form = betanull_consecutive.lf.max_mt_1.0, plot_title = "") {
  epochtemptype.plot.prop.df <- betanull_long_form %>%
  select(Site1, Site2, Assembly_Process, EpochType, TempConditionType, AvgDepth, AvgDepth.Site1) %>%
  mutate(DepthRange = paste0(round(AvgDepth, 1), "-", round(AvgDepth.Site1, 1)),
         DepthRange = fct_reorder(DepthRange, AvgDepth)) %>%
  group_by(DepthRange, EpochType, TempConditionType, Assembly_Process) %>%
  tally() %>%  
  ungroup() %>%
  mutate(DescreteDepth = fct_rev(DepthRange))

#temperatureinfo <- epochtemptype.plot.prop.df %>% 
#  select(DescreteDepth, EpochType, TempConditionType, n) %>%
#  distinct()

# Plot figure
temptype.plot.prop.consec <- ggplot(epochtemptype.plot.prop.df, aes(x = TempConditionType)) +
  guides(fill = guide_legend(nrow = 2)) + 
  scale_fill_manual(name = "Assembly Process", values = fill_assembly, 
                    breaks = assembly_levels, 
                    labels = assembly_labels, drop = FALSE) +
  geom_bar( aes(fill = Assembly_Process), position = "fill") +
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

# legend_plot <- cowplot::get_legend(temptype.plot.prop.consec)

temptype.plot.prop.consec <- temptype.plot.prop.consec + 
  theme(legend.position = "none") + ggtitle(plot_title)

return(list(figure = temptype.plot.prop.consec, 
            data = epochtemptype.plot.prop.df))
}

temp_change_plot_list <- purrr::map(otu_filter_levels, ~plot_tempchange(betanull_long_form = betanull_list_consecutive[[.x]], plot_title = .x) %>% purrr::pluck("figure")) 

names(temp_change_plot_list) <- otu_filter_levels
```

``` r
ggpubr::ggarrange(plotlist = temp_change_plot_list, common.legend = T) # need to fix legend....
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-30-1.png" style="width:98.0%;height:98.0%" />

-   The decrease in homogeneous selection that we observed from the previous table is not evenly split across warm-cold transitions. And the patterns seem to be different depending if the mean or max filtering option is applied.

-   I am not totally sure how to interpret the patterns we see. First, when filtering is at its most strict against rare taxa (mean, \>1%), drift dominates all transitions, with some homogenizing dispersal (high dispersal) in transitions that include warm years (warm-warm or warm-cold) and heterogeneous selection in cold years (recall that we did not previously observe heterogeneous selection –\> often it’s interpreted to mean competition among sister taxa). However, I would not read too much into these patterns as both the heterogeneous selection and homogenizing dispersal are only 2 and 3 observations, respectively.

-   The patterns for the mean-filtering do not hold when filtering is relaxed to 0.1%. The most we can say for them is that they are dominated by stochastic processes (drift, dispersal/dispersal limitation)

-   I think part of the challenge here is that it is impossible to know what constitutes the “ephemeral” dust-deposited community and what part of the community is the “resident” community.

-   Therefore, I am not particularly comfortable drawing much from these plots alone.

### Does dust show a relationship to deterministic processes when max filtering is applied?

If the strong filtering we see in the max is due to years in which dust is high, we might expect that when max filtering is applied, a relationship can be observed between dust amount and deterministic filtering.

``` r
dust_data_betanull <- dust_data %>% nest() %>%
  # prepare data for plotting and testing
  mutate(dust_data_betanull = purrr::map(data, ~compute_dust_difference(betanull = betanull_list_consecutive$max_mt_0.1, # use consecutive comparisons only
                                                                     sample_dust_data = .x,
                                                                     dust_column = "dust_amount"))) %>%
  # compute correlation of betaNTI
  mutate(betanti_corr = purrr::map(dust_data_betanull, ~filter(.x, AssemblyMetric == "BetaNTI") %>%
                                     select(AssemblyValue, DustDiff) %>%
                                     cor.test(~ AssemblyValue + DustDiff, data = ., method = "pearson"))) %>%
  # compute correlation of RCBC
  mutate(RCBC_corr = purrr::map(dust_data_betanull, ~filter(.x, AssemblyMetric == "RCBC") %>%
                                     select(AssemblyValue, DustDiff) %>%
                                     cor.test(~ AssemblyValue + DustDiff, data = ., method = "pearson"))) %>%
  # Plot
  mutate(dust_size_plot = purrr::map2(dust_data_betanull, dust_size_class, ~ggplot(.x, aes(x = DustDiff, y = AssemblyValue)) +
                                       geom_point(aes(color = Assembly_Process), size = rel(2)) +
                                       scale_color_manual(name = "Assembly Process", values = fill_assembly,
                                                          breaks = assembly_levels,
                                                          labels = assembly_labels) +
                                       facet_wrap(~AssemblyMetric, scales = "free_y") +
                                       theme_bw() + ggtitle(.y) + theme(legend.position = "none", aspect.ratio = 1)))
```

``` r
# now pull out the correlations

#response_plots <- purrr::map(flatten(dust_data_betanull$dust_size_plot), ~cowplot::plot_grid(plotlist = .x))

ggpubr::ggarrange(plotlist = dust_data_betanull$dust_size_plot, ncol = 1)
```

<img src="04_interpret_beta_null_results.markdown_strict_files/figure-markdown_strict/unnamed-chunk-32-1.png" style="width:98.0%;height:98.0%" />

``` r

# now pull out the correlations
dust_levels <- dust_data_betanull$dust_size_class
dust_betanti_corr <- dust_data_betanull %>% 
  mutate(betanti_tidy = purrr::map(betanti_corr, ~broom::tidy(.x))) %>%
  unnest(betanti_tidy) %>%
  select(dust_size_class, estimate:conf.high) %>%
  mutate(dust_size_class = factor(dust_size_class, levels = dust_levels))
  
dust_betanti_corr %>% 
  kbl(digits = 3) %>%
  column_spec(4, bold = ifelse((dust_betanti_corr$p.value > 0.05 | is.na(dust_betanti_corr$p.value)), FALSE, TRUE)) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| dust_size_class   | estimate | statistic | p.value | parameter | conf.low | conf.high |
|:------------------|---------:|----------:|--------:|----------:|---------:|----------:|
| dust_0.63_0.8µm   |    0.392 |     1.954 |   0.064 |        21 |   -0.024 |     0.692 |
| dust_0.8_1µm      |    0.392 |     1.952 |   0.064 |        21 |   -0.024 |     0.692 |
| dust_1_1.26µm     |    0.388 |     1.929 |   0.067 |        21 |   -0.029 |     0.690 |
| dust_1.26_1.59µm  |    0.385 |     1.914 |   0.069 |        21 |   -0.032 |     0.688 |
| dust_1.59_2.02µm  |    0.368 |     1.816 |   0.084 |        21 |   -0.052 |     0.678 |
| dust_2.02_2.52µm  |    0.344 |     1.678 |   0.108 |        21 |   -0.080 |     0.662 |
| dust_2.52_3.14µm  |    0.303 |     1.457 |   0.160 |        21 |   -0.125 |     0.636 |
| dust_3.14_4µm     |    0.239 |     1.128 |   0.272 |        21 |   -0.192 |     0.593 |
| dust_4_5.04µm     |    0.135 |     0.626 |   0.538 |        21 |   -0.293 |     0.519 |
| dust_5.04_6.35µm  |    0.032 |     0.146 |   0.886 |        21 |   -0.386 |     0.438 |
| dust_6.35_8µm     |   -0.084 |    -0.386 |   0.703 |        21 |   -0.480 |     0.340 |
| dust_8_10.08µm    |   -0.114 |    -0.527 |   0.604 |        21 |   -0.503 |     0.313 |
| dust_10.08_12.7µm |   -0.057 |    -0.260 |   0.797 |        21 |   -0.458 |     0.364 |
| dust_12.7_16.0    |   -0.024 |    -0.110 |   0.914 |        21 |   -0.432 |     0.392 |
| dust_total        |    0.386 |     1.915 |   0.069 |        21 |   -0.032 |     0.688 |

I think the answer is a resounding “no”. When dust differences are greater, the max-filtering method is not more likely to show strong abiotic filtering, even for small sizes of dust. This could mean that the rare taxa that are being filtered are not coming from dust, or that dust is a poor proxy for aeolian deposition of microbes. Other possible sources of rare taxa might include other forms of dispersal for example hydrologically-mediated dispersal (from rivulets on the surface of the glacier), blooms of rare taxa that occur due to reasons we cannot detect at this temporal resolution, relic DNA.

## Concluding Thoughts

-   We see strong evidence of abiotic selective processes driving the structure of these communities. Ecological drift is a second dominant assembly process. There is little evidence of biotic structuring within the community.
-   This doesn’t mean that there is no biotic structuring, just that it doesn’t show up as a sister-taxa competition. Other forms of biotic structuring, such as food-web type structuring cannot be ruled out via this method. What this does imply is that the conditions at hand, strongly favor particular sets of clades.
-   Ecological drift increases with time, although when the analysis is restricted to consecutive comparisons only, this increase is less linear. This is hard to attribute. We cannot rule out the influence of relic DNA, or active entrained cells in this pattern. However one piece of evidence against active entrained cells is a lack of “dispersal limitation with drift”, which I would expect to see more of if cells were actively growing in these ice cores. What is clear is that there are stronger abiotic selective pressures nearer the top of the core.
-   Climate epoch has a stronger influence on assembly process than temperature variation within an epoch. This may be due to different environmental/ecological forces being at play during each epoch. For example, climatic shifts, deposition shifts, or different sources of microbial communities.
-   In particular, the Holocene is a very selection-dominated time, despite the capture of multiple and Warm-Warm, Warm-Cold transitions. In the LGS, by contrast, Warm-Warm and Warm-Cold transitions are dominated by more stochastic processes, drift, and high dispersal. A similar pattern is seen in the Pre-LGS period, but the Warm-cold/Warm-Warm transitions are not consistent between them in amount of selection at play.
-   Microbial communities are significantly different in warm vs. cold periods. Interestingly cold seems to impart a greater selective pressure than warmth (though this is really only apparent prior to the Holocene). This may align with the role of photorophs in structuring the community, since colder years likely include less available light for phototrophs which could have downstream selective pressures on the community.
-   Rare taxa seem to contribute strongly to the homogenizing (abiotic) filtering signal that we see. When removed, drift and other stochastic processes dominate. This suggests that the rare taxa are those that experience the strongest abiotic filtering of the community. Although why this is a question that we may only be able to hypothesize about. Potential options include:
    -   1\) strong filtering of the dust depositional community
    -   2\) intermittent blooms from rare taxa that experience strong filtering
