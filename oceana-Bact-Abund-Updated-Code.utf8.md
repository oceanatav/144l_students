---
title: "Bact Abund Update"
author: "Oceana Tavasieff"
date: "11/6/2020"
output: html_document
---

---
title: "Bacterial Abundance Oceana"
author: "Oceana Tavasieff"
date: "10/21/2020"
output: github_document
---



# Intro

In which data on individual bottle bacterial abundance data from ACIDD experiements is processed, quality controlled, and analyzed. 


```r
library(tidyverse)
library(readxl)
library(lubridate)
```

# Import data


```r
excel_sheets("~/Documents/github_144l/144l_students/Input_Data/week3/144L_2018_BactAbund.xlsx")
```

```
## [1] "Metadata" "Data"
```

```r
metadata <- read_excel("~/Documents/github_144l/144l_students/Input_Data/week3/144L_2018_BactAbund.xlsx", sheet = "Metadata")

#glimpse(metadata)

unique(metadata$Experiment) #not unique
```

```
## [1] "144L_2018"
```

```r
unique(metadata$Location) #campus point only
```

```
## [1] "Campus Point"
```

```r
unique(metadata$Bottle) #a-h
```

```
## [1] "A" "B" "C" "D" "E" "F" "G" "H"
```

```r
unique(metadata$Treatment) # 4 treatments
```

```
## [1] "Control"                   "Ash Leachate"             
## [3] "Mud Leachate"              "Glucose_Nitrate_Phosphate"
```

```r
data <- read_excel("~/Documents/github_144l/144l_students/Input_Data/week3/144L_2018_BactAbund.xlsx", sheet = "Data")

glimpse(data)
```

```
## Rows: 72
## Columns: 3
## $ Bottle    <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B...
## $ Timepoint <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, ...
## $ Cells_ml  <dbl> 332531.5, 523943.1, 859019.9, 906998.9, 933025.2, 861129....
```

```r
# Shared columns between these two sheets (Data and Metadata).
# Combine into one datasheet

joined <- left_join (metadata, data) #left_join joins the right dataset to the left dataset
```

```
## Joining, by = c("Bottle", "Timepoint")
```

```r
names(joined) # Check our join function
```

```
##  [1] "Experiment"           "Location"             "Temperature"         
##  [4] "Depth"                "Bottle"               "Timepoint"           
##  [7] "Treatment"            "Target_DOC_Amendment" "Inoculum_L"          
## [10] "Media_L"              "Datetime"             "TOC_Sample"          
## [13] "Parallel_Sample"      "Cell_Sample"          "DNA_Sample"          
## [16] "DNA_SampleID"         "Cells_ml"
```

```r
glimpse(joined)
```

```
## Rows: 80
## Columns: 17
## $ Experiment           <chr> "144L_2018", "144L_2018", "144L_2018", "144L_2...
## $ Location             <chr> "Campus Point", "Campus Point", "Campus Point"...
## $ Temperature          <dbl> 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20...
## $ Depth                <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1...
## $ Bottle               <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "...
## $ Timepoint            <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5...
## $ Treatment            <chr> "Control", "Control", "Control", "Control", "C...
## $ Target_DOC_Amendment <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0...
## $ Inoculum_L           <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1...
## $ Media_L              <dbl> 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3...
## $ Datetime             <chr> "2018-10-15T16:30", "2018-10-16T08:00", "2018-...
## $ TOC_Sample           <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,...
## $ Parallel_Sample      <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,...
## $ Cell_Sample          <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE...
## $ DNA_Sample           <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,...
## $ DNA_SampleID         <chr> "144_A0_S6", NA, NA, NA, "144_A4_S7", NA, NA, ...
## $ Cells_ml             <dbl> 332531.5, 523943.1, 859019.9, 906998.9, 933025...
```

```r
#summary(joined)
```

# Prepare Data

Convert date and time columns from character to date values. Add columns with time elapsed for each experiment (subsetting data by experiment).
Cells/mL to Cells/L (DOC units are in umol C/L). Remove extraneous columns. Drop NA's



```r
cells <- joined %>%  
  mutate(Datetime= ymd_hm(Datetime),
         cells = Cells_ml * 1000,
         sd_cells = sd(cells)) %>% 
  group_by(Treatment, Bottle) %>%  
  mutate(interv = interval(start = first(Datetime), end = Datetime),
         hours = interv/3600, 
         days = hours/24) %>% ungroup() %>%
  select(Location:days, hours, days, cells, sd_cells) %>% drop_na(cells)


glimpse(cells)
```

```
## Rows: 72
## Columns: 21
## $ Location             <chr> "Campus Point", "Campus Point", "Campus Point"...
## $ Temperature          <dbl> 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20...
## $ Depth                <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1...
## $ Bottle               <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "...
## $ Timepoint            <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6...
## $ Treatment            <chr> "Control", "Control", "Control", "Control", "C...
## $ Target_DOC_Amendment <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0...
## $ Inoculum_L           <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1...
## $ Media_L              <dbl> 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3...
## $ Datetime             <dttm> 2018-10-15 16:30:00, 2018-10-16 08:00:00, 201...
## $ TOC_Sample           <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,...
## $ Parallel_Sample      <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,...
## $ Cell_Sample          <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE...
## $ DNA_Sample           <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,...
## $ DNA_SampleID         <chr> "144_A0_S6", NA, NA, NA, "144_A4_S7", NA, NA, ...
## $ Cells_ml             <dbl> 332531.5, 523943.1, 859019.9, 906998.9, 933025...
## $ cells                <dbl> 332531522, 523943125, 859019934, 906998856, 93...
## $ sd_cells             <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA...
## $ interv               <Interval> 2018-10-15 16:30:00 UTC--2018-10-15 16:30...
## $ hours                <Interval> 2018-10-15 16:30:00 UTC--2018-10-15 16:30...
## $ days                 <Interval> 2018-10-15 16:30:00 UTC--2018-10-15 16:30...
```

```r
#interval fails to return an interval... or rather a period
```

# Plot Growth Curves 


```r
custom.colors <- c("Control" = "#9eb0a3", "Ash Leachate" = "#ebb688", "Mud Leachate" = "#55a3a2", "Glucose_Nitrate_Phosphate" = "#bec46a", "Vial" = "#e0e067", "Bottle" = "#8797fa" )
levels <- c("Control", "Ash Leachate", "Mud Leachate", "Glucose_Nitrate_Phosphate", "Bottle", "Vial")



cells %>% 
  ggplot(aes(x = days, y = cells, group = interaction(Treatment, Bottle))) +
  geom_errorbar(aes(ymin = cells-sd_cells, ymax= cells + sd_cells, color= factor(Treatment, levels = levels)), width=0.1) + 
  geom_line(aes(color= factor (Treatment, levels= levels)), size = 1) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size= 3, color= "black", shape = 21) +
  labs(x = "Days", y= expression(paste("Cells, L"^-1)), fill = "") +
  guides(fill = F) +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  #facet_grid(rows = "", scales = "free") +
  theme_bw()
```

<img src="oceana-Bact-Abund-Updated-Code_files/figure-html/unnamed-chunk-4-1.png" width="672" />



```r
# Lets add a marker for where DNA sample = T (points where we know what the community composition is)
# Set * as the symbol (in geom_text; we also had to create a new column as the * if dna sample was taken; if not, then it printes NA

cells %>% 
  mutate(dna =  ifelse(DNA_Sample == T, "*", NA)) %>% 
  ggplot(aes(x = days, y = cells, group = interaction(Treatment, Bottle))) +
  geom_errorbar(aes(ymin = cells-sd_cells, ymax= cells + sd_cells, color= factor(Treatment, levels = levels)), width=0.1) + 
  geom_line(aes(color= factor (Treatment, levels= levels)), size = 1) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size= 3, color= "black", shape = 21) +
  geom_text(aes(label = dna), size = 10, color = "#e41a1c") + 
  labs(x = "Days", y= expression(paste("Cells, L"^-1)), fill = "") +
  guides(fill = F) +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  facet_grid("Bottle", scales = "free") +
  theme_bw()
```

```
## Warning: Removed 48 rows containing missing values (geom_text).
```

<img src="oceana-Bact-Abund-Updated-Code_files/figure-html/unnamed-chunk-5-1.png" width="576" />

We can calculate:
  -Total change in cells from initial conditions to end of expt
  - Specific growth rate as slope of ln(abundance) vs time during exptl growth phase
  - Doubling time as ln(2) divided by spec growth rate
  -mean of each of these parameters for each treatment

First we need to det. where exptl growth phase occurs for each expt (if it does). Let's plot ln(abund) vs time. 

## Identify Exponential Growth Rates


Bottle
A ~0-1d (T0-2)
B ~ 0-2 (T0-4)
C ~ 0-1.5 (T0-3)
D ~ 0-2 (T0-4)
E ~ 0-1.5 (T0-3)
F ~ 0-1
G ~ 0-1
H ~ 0-1

Recall:
log(x) give ln(x)
log10(x)  gives log base 10
log2(x) gives log base 2


```r
ln_cells <- cells %>% group_by(Treatment, Bottle) %>% 
  mutate(ln_cells = log(cells),
         diff_ln_cells = ln_cells - lag(ln_cells, default = first(ln_cells))) %>% ungroup()
```


```r
# View exponential growth


ln_cells %>% 
  mutate(dna =  ifelse(DNA_Sample == T, "*", NA)) %>% 
  ggplot(aes(x = days, y = diff_ln_cells, group = interaction(Treatment, Bottle))) + 
  geom_line(aes(color= factor(Treatment, levels = levels)), size = 1) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size= 3, color= "black", shape = 21) +
  geom_text(aes(label = dna), size = 10, color = "#e41a1c") + 
  labs(x = "Days", y= expression(paste("Change in ln(Cells)")), fill = "") +
  guides(fill = F) +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  facet_grid(rows = "Treatment", scales = "free") +
  theme_bw()
```

```
## Warning: Removed 48 rows containing missing values (geom_text).
```

<img src="oceana-Bact-Abund-Updated-Code_files/figure-html/unnamed-chunk-7-1.png" width="672" />

# Calculate Growth Rates, Doubling Times, and Change in Cell Abundance

```r
# ASH171 = Expts in SD, ASH172 = Expts in SB
#growth <- ln_cells %>% 
 # mutate(exp_start = ifelse(Experiment == "ASH171" & Bottle == "A", 4, NA))
      # Adds a column that adds the Timepoint of Exptl growth. for Bottle A, Experiment ash171 (in San Diego), we visually determined that Exptl Growth occurred between 3-5d, then we set midpt/ Timepoint to be 4d. 
# ended exp_start lines with a 4 because many bottles end at Timepoint 4
growth <- ln_cells %>% 
  mutate(exp_start = ifelse(Bottle %in% c("C", "D", "E", "F"), 0, NA),
         exp_start = ifelse(Bottle %in% c("A", "B", "G", "H"), 1, exp_start),
         exp_end = ifelse(Bottle %in% c("A", "F", "G", "H"), 2, NA),
         exp_end = ifelse(Bottle %in% c("B", "C", "E"), 4, exp_end),
         exp_end = ifelse(Bottle == "D", 3, exp_end))


check <- growth %>% select (Bottle, exp_start, exp_end) %>% distinct() #  double check our work above
```


```r
# MUTATE to input ln(cell abund) for the start of every exponential growth curve. recall this start corresponds with a particular Timepoint value
growth <- ln_cells %>% 
  mutate(exp_start = ifelse(Bottle %in% c("C", "D", "E", "F"), 0, NA),
         exp_start = ifelse(Bottle %in% c("A", "B", "G", "H"), 1, exp_start),
         exp_end = ifelse(Bottle %in% c("A", "F", "G", "H"), 2, NA),
         exp_end = ifelse(Bottle %in% c("B", "C", "E"), 4, exp_end),
         exp_end = ifelse(Bottle == "D", 3, exp_end)) %>% 
  mutate( ln_cells_exp_start = ifelse(Timepoint == exp_start, ln_cells, NA),
          ln_cells_exp_end = ifelse(Timepoint == exp_end, ln_cells, NA),
          cells_exp_start = ifelse (Timepoint == exp_start, cells, NA),
          cells_exp_end = ifelse(Timepoint == exp_end, cells, NA),
          days_exp_start = ifelse (Timepoint == exp_start, days, NA),
          days_exp_end = ifelse(Timepoint == exp_end, days, NA)) %>% 
  fill(ln_cells_exp_start:days_exp_end, .direction= "updown") %>% 
  mutate(mew = (ln_cells_exp_end - ln_cells_exp_start) / ( days_exp_end - days_exp_start),
         doubling = log(2)/mew,
         delta_cells = cells_exp_end - first(cells)) %>% 
  ungroup()

# recall our variables are grouped and so this fill fxn fills only along experiment in this use
```

# Convert Bacterial Abundance and Change in Bact. Abund. to Carbon Units

Apply a Carbon Conversion Factor (CCF) to bact. abundances (cells L^-1^) to generate bacterial carbon (umol C L^-1^).

We will apply average carbon conversion of bacterioplankton cells from Coastal JP (~30fg C cell^-1^), reported by Fukuda et al., 1998. 


```r
bactcarbon <- growth %>% 
  mutate(bc = cells * (2.5 * 10^-9),
         delta_bc = delta_cells * (2.5 * 10^-9))

glimpse(bactcarbon)
```

```
## Rows: 72
## Columns: 36
## $ Location             <chr> "Campus Point", "Campus Point", "Campus Point"...
## $ Temperature          <dbl> 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20...
## $ Depth                <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1...
## $ Bottle               <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "...
## $ Timepoint            <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6...
## $ Treatment            <chr> "Control", "Control", "Control", "Control", "C...
## $ Target_DOC_Amendment <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0...
## $ Inoculum_L           <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1...
## $ Media_L              <dbl> 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3...
## $ Datetime             <dttm> 2018-10-15 16:30:00, 2018-10-16 08:00:00, 201...
## $ TOC_Sample           <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,...
## $ Parallel_Sample      <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,...
## $ Cell_Sample          <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE...
## $ DNA_Sample           <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,...
## $ DNA_SampleID         <chr> "144_A0_S6", NA, NA, NA, "144_A4_S7", NA, NA, ...
## $ Cells_ml             <dbl> 332531.5, 523943.1, 859019.9, 906998.9, 933025...
## $ cells                <dbl> 332531522, 523943125, 859019934, 906998856, 93...
## $ sd_cells             <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA...
## $ interv               <Interval> 2018-10-15 16:30:00 UTC--2018-10-15 16:30...
## $ hours                <Interval> 2018-10-15 16:30:00 UTC--2018-10-15 16:30...
## $ days                 <Interval> 2018-10-15 16:30:00 UTC--2018-10-15 16:30...
## $ ln_cells             <dbl> 19.62225, 20.07689, 20.57130, 20.62565, 20.653...
## $ diff_ln_cells        <dbl> 0.000000000, 0.454648479, 0.494408990, 0.05434...
## $ exp_start            <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1...
## $ exp_end              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4...
## $ ln_cells_exp_start   <dbl> 20.07689, 20.07689, 20.10741, 20.10741, 20.107...
## $ ln_cells_exp_end     <dbl> 20.57130, 20.57130, 20.57130, 20.82448, 20.824...
## $ cells_exp_start      <dbl> 523943125, 523943125, 540180187, 540180187, 54...
## $ cells_exp_end        <dbl> 859019934, 859019934, 859019934, 1106511579, 1...
## $ days_exp_start       <dbl> 0.6458333, 0.6458333, 0.6458333, 0.6458333, 0....
## $ days_exp_end         <dbl> 1.145833, 1.145833, 1.145833, 2.145833, 2.1458...
## $ mew                  <dbl> 0.9888180, 0.9888180, 0.9277787, 0.4780432, 0....
## $ doubling             <dbl> 0.7009856, 0.7009856, 0.7471040, 1.4499675, 1....
## $ delta_cells          <dbl> 526488412, 526488412, 526488412, 773980057, 77...
## $ bc                   <dbl> 0.8313288, 1.3098578, 2.1475498, 2.2674971, 2....
## $ delta_bc             <dbl> 1.316221, 1.316221, 1.316221, 1.934950, 1.9349...
```



# Calculate Treatment Averages


```r
#when calculating averages, group data 1st!

averages <- bactcarbon %>% 
  group_by(Treatment, Bottle, Timepoint) %>% 
  mutate(ave_bc = mean(bc),
         sd_bc = sd(bc)) %>% 
  ungroup() %>% 
  group_by(Treatment, Bottle) %>% 
  mutate(ave_mew = mean(mew),
         sd_mew = sd(mew),
         ave_doubling =  mean(doubling),
         sd_doubling = sd(doubling),
         ave_delta_cells = mean(delta_cells),
         sd_delta_cells = sd(delta_cells),
         ave_delta_bc = mean(delta_bc),
         sd_delta_bc = sd(delta_bc),
         ave_lag = mean(days_exp_start),
         sd_lag = sd(days_exp_start)) %>% 
  ungroup()

subset <- averages %>% select(Treatment,  Bottle, Timepoint,bc, ave_bc, sd_bc)
```

## Plot Treatment Averages


```r
averages %>% 
  ggplot(aes(x = days, y = ave_bc), group = interaction(Experiment, Treatment)) +
  geom_errorbar(aes(ymin = ave_bc - sd_bc, ymax = ave_bc + sd_bc, color= factor(Treatment, levels = levels)), width = 0.1) +
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), color = "black", shape = 21, size = 2) +
  facet_grid(rows = "Bottle", scales = "free") +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) + labs(x = "Days", y = expression("Bacterial Carbon, umol C L"^-1), fill = "", color = "") + guides(color = F) +
  theme_bw()
```

<img src="oceana-Bact-Abund-Updated-Code_files/figure-html/unnamed-chunk-12-1.png" width="576" />

# Barplots


```r
bar.data <- averages %>% 
  select(Location, Treatment, ave_mew:sd_lag)
# many duplicated rows, repeat with distinct() fxn to cut down dataset

bar.data <- averages %>% 
  select(Location, Treatment, ave_mew:sd_lag) %>% distinct()
```



```r
mew <- bar.data %>%
  ggplot(aes(x = factor(Treatment, levels = levels), y =  ave_mew),
         group = interaction(Location, Treatment)) + 
  geom_col(color = "black", fill = "white") +
  geom_errorbar(aes(ymin = ave_mew - sd_mew, ymax = ave_mew + sd_mew), width = 0.1) +
  facet_grid(~factor(Location, levels = levels), scales = "free") +
  labs(x = "", y = expression("μ, d"^-1)) +
  theme_bw()

# FACET_GRID to grid by location
# use tilde before factor in place of rows = or col =  argument
```


```r
doubling <- bar.data %>%
  ggplot(aes(x = factor(Treatment, levels = levels), y =  ave_doubling),
         group = interaction(Location, Treatment)) + 
  geom_col(color = "black", fill = "white") +
  geom_errorbar(aes(ymin = ave_doubling - sd_doubling, ymax = ave_doubling + sd_doubling), width = 0.1) +
  facet_grid(~factor(Location, levels = levels), scales = "free") +
  labs(x = "", y = expression("Doubling Time, d"^-1)) +
  theme_bw()
```



```r
delta_bc <- bar.data %>%
  ggplot(aes(x = factor(Treatment, levels = levels), y =  ave_delta_bc),
         group = interaction(Location, Treatment)) + 
  geom_col(color = "black", fill = "white") +
  geom_errorbar(aes(ymin = ave_delta_bc - sd_delta_bc, ymax = ave_delta_bc + sd_delta_bc), width = 0.1) +
  facet_grid(~factor(Location, levels = levels), scales = "free") +
  labs(x = "", y = expression("Δ Bacterial Carbon, μmol C L"^-1)) +
  theme_bw()
```


```r
lag <- bar.data %>%
  ggplot(aes(x = factor(Treatment, levels = levels), y =  ave_lag),
         group = interaction(Location, Treatment)) + 
  geom_col(color = "black", fill = "white") +
  geom_errorbar(aes(ymin = ave_lag - sd_lag, ymax = ave_lag + sd_lag), width = 0.1) +
  facet_grid(~factor(Location, levels = levels), scales = "free") +
  labs(x = "", y = "Lag Phase, days") +
  theme_bw()
```

## Attach all barplots together


```r
#install.packages("patchwork")
library(patchwork)
```

```
## Warning: package 'patchwork' was built under R version 4.0.3
```


```r
lag + delta_bc + mew + doubling + plot_annotation(tag_levels = "a")
```

<img src="oceana-Bact-Abund-Updated-Code_files/figure-html/unnamed-chunk-15-1.png" width="672" />

# Save Data


```r
saveRDS(averages, "~/Documents/github_144l/144l_students/Output_Data/Week 3/ACIDD_Exp_Processed_BactAbund_rds")

#save_csv
```

