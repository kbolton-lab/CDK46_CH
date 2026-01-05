# CDK4/6 inhibition mitigates chemotherapy-induced expansion of TP53-mutant clonal hematopoiesis
Therapy-related myeloid neoplasm (tMN) represents a fatal consequence of exposure to cytotoxic therapy administered in the treatment of cancer. Individuals with pre-existing TP53 clonal hematopoiesis (CH) are at high risk of tMN but currently avoidance of therapy is the only strategy to reduce tMN risk. Here, we show in four randomized clinical trials of the CDK4/6 inhibitor, trilaciclib, given in conjunction with a variety of chemotherapeutic regimens and across diverse cancer patient populations that trilaciclib mitigates chemotherapy-related expansion of CH clones with mutations in DNA damage response (DDR) genes (including TP53, PPM1D, and CHEK2 mutations). This finding was also observed in a syngeneic murine model of TP53 mutant CH demonstrating that trilaciclib blocks platinum-induced TP53 competitive repopulation through promoting hematopoietic stem and progenitor quiescence and decreasing TP53 mutant stemness advantage. This represents the first demonstration of a pharmacologic strategy to block chemotherapy-induced expansion of pre-leukemic TP53-mutant clones.

The full manuscript can be found here: <link>

# System Requirements
## OS Requirements
This code is supported for all operating systems that can run the R coding language. This code has been tested on:
- macOS: Sequoia v15.5
- Linux: Ubuntu v20.04.6

## R Version
```
R v4.5.1 (2025-06-13)
```

## R Packages:
```
seqinr v4.2-36
biomaRt v2.64.0
GenVisR v1.39.0
conflicted v1.2.0
RColorBrewer v1.1-3
maftools v2.24.0
readxl v1.4.5
sjPlot v2.9.0
cowplot v1.2.0
ggsignif v0.6.4
ggsci v3.2.0
ggpubr v0.6.1
ComplexHeatmap v2.24.1
table1 v1.4.3
data.table v1.17.8
lubridate v1.9.4
forcats v1.0.0
stringr v1.5.1
dplyr v1.1.4
purrr v1.1.0
readr v2.1.5
tidyr v1.3.1
tibble v3.3.0
ggplot2 v3.5.2
tidyverse v2.0.0
```

# Installation
## Installing R
R can be installed quickly (few minutes) following the instructions from: https://cran.rstudio.com/

## Installing Packages
The R packages above can be installed through a combination of `install.packages("package_name")` or using the Bioconductor package `BiocManager::install("package_name")`

# Using this R Markdown File
```
git clone https://github.com/kbolton-lab/CDK46_CH/
cd CDK46_CH
R -e "rmarkdown::render('CDK46_CH.Rmd')"
```
Depending on various processing time and connection to the biomaRt server, can take around 10 to 20 minutes.

# Data
## Clonal Hematopoiesis (CH) Variant Calls:
- untreated_df.minimum.csv
- sclc_df.minimum.csv
- crc_df.minimum.csv
- tnbc_df.minimum.csv

## Serial Timepoint Values:
- untreated_timepoints.under17000.minimum.csv
- sclc_timepoints.under17000.minimum.csv
- crc_timepoints.under17000.minimum.csv
- tncb_timepoints.under17000.minimum.csv

## Publication Workflow Images
- G1_SCLC.png
- G1_CRC.png
- G1_TNBC.png

## Demographic Data
- g1_demo.csv
- g1_labresults.csv

## Data used for Supplemental Tables
- cdk46_ch_pilot_harvard.filtered.csv

# Output
After installation and running this R Markdown file, users can expect to have all of the figures as well as Supplemental Tables found in the manuscript generated and placed into the `figures` directory

