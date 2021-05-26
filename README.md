# GlobalChangeSpace
Code and data supplement for the manuscript "Representation of global change across biodiversity databases" by Daskalova et al.

This repository is maintained by Gergana Daskalova. For any questions, please contact gndaskalova@gmail.com.

### Abstract 
Global change has altered biodiversity and impacted ecosystem functions and services around the planet. Understanding the effects of anthropogenic drivers like human use and climate change and predicting trajectories of future biodiversity change have become key challenges for science and policy. However, our knowledge of biodiversity change is limited by the available data and their biases. In both the terrestrial and marine realms, we test the representation of three worldwide biodiversity databases (Living Planet, BioTIME and PREDICTS) across geographic and temporal variation in global change and across the tree of life. We find that variation in global change drivers is better captured over space than over time and in the marine realm versus on land. We provide recommendations to improve the use of existing data, better target future ecological monitoring and capture different combinations of global change.

### Repository structure and code organisation
All R scripts can be found in the code folder. They are sorted based on order of operations, but each script is standalone and can be run without running the previous ones. The data are sorted in input and output folders. A pre-registration for this study is available in the doc folder, as well a draft manuscript. Note that this manuscript represents a side project that developed out of the pre-registration (i.e., the pre-registration covers content for this project as well as another project). The figures folder contains way too many figures and will be organised later.

### Scripts

```
- 01-extract-drivers-over-space.R  # To extract global change driver intensity for five drivers across the sites represented in the Living Planet, BioTIME and PREDICTS databases, as well as randomly around the world 
- 02-extract-terr-temp-over-time.R # To extract annual temperature over the duration of each time series in the Living Planet and BioTIME databases (terrestrial realm) 
- 03-extract-marine-temp-over-time.R # To extract annual temperature over the duration of each time series in the Living Planet and BioTIME databases (marine realm) 
- 04-extract-marine-ecoregions-lpd.R # To extract the marine ecoregions represented in the Living Planet Database 
- 05-run-models.R # To run statistical models 
- 06-calculate-predictions.R # To calculate model predictions 
- 07-plot-figures.R # To create figures
```

### Data

#### Living Planet Database

The population time series came from the Living Planet Database, publicly available from http://www.livingplanetindex.org

#### BioTIME

The biodiversity time series came from the BioTIME Database. Approximately 92% of the biodiversity studies analysed here are available as part of the published BioTIME Database. The data are openly available, and can be accessed on Zenodo (https://doi.org/10.5281/zenodo.1211105) or through the BioTIME website (http://biotime.st-andrews.ac.uk/). The remaining 8% of the studies were used with permission, with details on how to download those data available in the supplementary information.

The public studies that were included in the version of BioTIME we analyzed can be downloaded from http://biotime.st-andrews.ac.uk/BioTIME_download.php

For more information about the BioTIME database, please see:

Dornelas, M., L.H. Antao, F. Moyes, A.E. Bates, A.E. Magurran, and BioTIME consortium (200+ authors). 2018. BioTIME: a database of biodiversity time-series for the Anthropocene. Global Ecology and Biogeography. 10.1111/geb.12729

#### Global change driver data

We used the 16 marine and terrestrial global change driver layers compiled by Bowler et al. 2020. We selected these layers because they had been harmonized across both realms and hence were most suitable for our global analysis. As in Bowler et al., these layers were grouped into five focal drivers: human use (land-use for the terrestrial realm, and exploitation for the marine realm), climate change, human population density, pollution and invasion potential. The driver data are spatially-explicit, but current data limitations prevent us from extracting temporally-explicit values of the magnitudes of most of these drivers over time. For details on the individual layers forming the global change data, including their resolutions and temporal coverage, see Table S1 in Bowler et al. 2020.

Bowler, D. E., Bjorkman, A. D., Dornelas, M., Myers‐Smith, I. H., Navarro, L. M., Niamir, A., ... & Bates, A. E. (2020). Mapping human pressures on biodiversity across the planet uncovers anthropogenic threat complexes. People and Nature, 2(2), 380-394.

#### Land Use Harmonisation Dataset

Available from http://luh.umd.edu

#### Sea surface temperatures from NOAA

Available from https://psl.noaa.gov/repository/entry/show/PSD+Climate+Data+Repository/Public/PSD+Datasets/NOAA+OI+SST/Weekly+and+Monthly/sst.mnmean.nc?entryid=cac1c2a6-a864-4409-bb77-1fdead8eeb6e&output=default.html

#### Air temperatures from CRU

Available from https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.01/cruts.1709081022.v4.01/tmp/

### R and package versions

The session information is attached below. Some of the scripts take a while to run because of the large amounts of data.

```
R version 3.6.2 (2019-12-12)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] broom_0.7.2              cluster_2.1.0            ggbiplot_0.55            scales_1.1.1            
 [5] plyr_1.8.6               readxl_1.3.1             speciesgeocodeR_2.0-10   cowplot_1.0.0           
 [9] mapdata_2.3.0            maps_3.3.0               hrbrthemes_0.6.0         ggridges_0.5.2          
[13] RColorBrewer_1.1-2       ggfortify_0.4.11         stargazer_5.2.2          viridis_0.5.1           
[17] viridisLite_0.3.0        treemapify_2.5.3         ncdf4_1.17               ggthemes_4.2.0          
[21] ggalt_0.6.2              dggridR_2.0.6            rgdal_1.5-18             gridExtra_2.3           
[25] raster_3.3-13            sp_1.4-4                 factoextra_1.0.7.999     CoordinateCleaner_2.0-18
[29] ggeffects_0.14.1         sjPlot_2.8.2             parameters_0.5.0         sjstats_0.17.9          
[33] tidybayes_2.0.1          bayesplot_1.7.1          modelr_0.1.6             brms_2.12.0             
[37] Rcpp_1.0.5               rstan_2.19.3             StanHeaders_2.21.0-1     forcats_0.4.0           
[41] stringr_1.4.0            dplyr_1.0.2              purrr_0.3.4              readr_1.3.1             
[45] tidyr_1.1.2              tibble_3.0.4             ggplot2_3.3.2            tidyverse_1.3.0         

loaded via a namespace (and not attached):
  [1] tidyselect_1.1.0     rgbif_3.3.0          lme4_1.1-21          htmlwidgets_1.5.2    munsell_0.5.0       
  [6] codetools_0.2-16     effectsize_0.2.0     units_0.6-7          DT_0.17              miniUI_0.1.1.1      
 [11] withr_2.3.0          Brobdingnag_1.2-6    colorspace_2.0-0     knitr_1.28           uuid_0.1-4          
 [16] rstudioapi_0.13      stats4_3.6.2         Rttf2pt1_1.3.8       emmeans_1.4.4        conditionz_0.1.0    
 [21] oai_0.3.0            bridgesampling_1.0-0 coda_0.19-3          vctrs_0.3.5          generics_0.1.0      
 [26] xfun_0.12            R6_2.5.0             markdown_1.1         assertthat_0.2.1     promises_1.1.1      
 [31] rgeos_0.5-5          gtable_0.3.0         ash_1.0-15           processx_3.4.4       rlang_0.4.8         
 [36] systemfonts_0.1.1    splines_3.6.2        extrafontdb_1.0      lazyeval_0.2.2       inline_0.3.15       
 [41] yaml_2.2.1           reshape2_1.4.3       abind_1.4-5          threejs_0.3.3        crosstalk_1.1.0.1   
 [46] backports_1.2.0      httpuv_1.5.2         rsconnect_0.8.16     extrafont_0.17       tools_3.6.2         
 [51] geoaxe_0.1.0         ellipsis_0.3.1       base64enc_0.1-3      rnaturalearth_0.1.0  classInt_0.4-3      
 [56] ps_1.4.0             prettyunits_1.1.1    zoo_1.8-7            haven_2.2.0          ggrepel_0.8.1       
 [61] fs_1.3.1             magrittr_2.0.1       data.table_1.13.2    colourpicker_1.0     reprex_0.3.0        
 [66] mvtnorm_1.1-0        whisker_0.4          sjmisc_2.8.3         matrixStats_0.55.0   evaluate_0.14       
 [71] hms_0.5.3            shinyjs_1.1          mime_0.9             arrayhelpers_1.1-0   xtable_1.8-4        
 [76] shinystan_2.5.0      rstantools_2.0.0     compiler_3.6.2       KernSmooth_2.23-16   crayon_1.3.4        
 [81] minqa_1.2.4          htmltools_0.5.0      mgcv_1.8-31          later_1.1.0.1        lubridate_1.7.4     
 [86] DBI_1.1.0            sjlabelled_1.1.3     dbplyr_1.4.2         proj4_1.0-10         MASS_7.3-51.4       
 [91] sf_0.9-6             boot_1.3-23          Matrix_1.2-18        permute_0.9-5        cli_2.2.0           
 [96] parallel_3.6.2       insight_0.8.1        igraph_1.2.4.2       pkgconfig_2.0.3      geosphere_1.5-10    
[101] xml2_1.3.2           svUnit_0.7-12        dygraphs_1.1.1.6     picante_1.8.1        estimability_1.3    
[106] rvest_0.3.5          callr_3.5.1          digest_0.6.27        vegan_2.5-6          rmarkdown_2.1       
[111] cellranger_1.1.0     gdtools_0.2.1        shiny_1.4.0          gtools_3.8.1         nloptr_1.2.1        
[116] lifecycle_0.2.0      nlme_3.1-142         jsonlite_1.7.1       fansi_0.4.1          pillar_1.4.7        
[121] lattice_0.20-38      loo_2.2.0            fastmap_1.0.1        httr_1.4.2           pkgbuild_1.1.0      
[126] glue_1.4.1.9000      xts_0.12-0           bayestestR_0.5.2     shinythemes_1.1.2    class_7.3-15        
[131] stringi_1.5.3        performance_0.4.4    ggfittext_0.9.0      ape_5.3              e1071_1.7-4         
```
