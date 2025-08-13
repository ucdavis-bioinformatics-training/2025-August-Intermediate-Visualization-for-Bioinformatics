# Advanced and custom plots

## TSS enrichments plots

## Heat maps

## Dendrograms

## Genomic coverage plots

# Set up

## Packages

In addition to ggplot2, we will use the ggplot extension libraries
ggtree and ggcoverage as well as ComplexHeatmap, which uses a completely
different “grammar.”

``` r
library(ggplot2)
library(ggtree)
library(ggcoverage)
library(ComplexHeatmap)
```

## Data

We will be making use of a few different data sets in this section: one
for each of the plots.

``` r
# ATAC-seq for TSS enrichment
# scRNA for heat map
# dendrogram data
# genomic coverage data
```

# TSS enrichment

# Heat map

# Dendrogram

# Genomic coverage

# Final thoughts

There are many types of visualizations we haven’t had the chance to
cover, and more visualizations are developed all the time. We hope that
what you’ve learned within this course will allow you to customize your
plots and figures to communicate your findings.

Here are a few things to think about as you choose how to present your
work using graphics:

- **Color** Make sure to choose colors that provide sufficient contrast
  between relevant data. Try to keep color palettes small; using a
  limited number of colors (typically no more than 12 within one figure)
  increases the readability of figures.

- **Size** Accurately estimating differences in area is difficult, so
  consider alternatives to pie charts and other area-based
  visualizations. Additionally, differences in size (e.g. of points)
  imply quantitative differences within the data; point size should not
  be mapped to a qualitative measure. Make sure that any plot elements
  are sized so that they are easily distinguishable. Ideally, they will
  be large enough to see clearly, but small enough to avoid extensive
  overlaps.

- **Redundancy** Each graphical attribute of a plot communicates
  information. When choosing which to use in any given figure, consider
  whether redundancy (e.g. using both point color and point shape to
  denote treatment group) increases clarity or introduces unnecessary
  complexity. Sometimes, using fewer aesthetics speeds reader
  comprehension.

- **Consistency** Maintain the order of categorical values on axes
  across figures. Any time a color palette is re-used within a
  publication, poster, or presentation, it should map to the same values
  to avoid confusion.

We hope this workshop has been helpful for you! Please email us at
<bioinformatics-training@ucdavis.edu> or reach out on GitHub with any
questions or comments about the course or course materials.

``` r
sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Monterey 12.4
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] ComplexHeatmap_2.24.1 ggcoverage_1.4.1      ggtree_3.16.3        
    ## [4] ggplot2_3.5.2        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1            dplyr_1.1.4                
    ##  [3] farver_2.1.2                Biostrings_2.76.0          
    ##  [5] bitops_1.0-9                fastmap_1.2.0              
    ##  [7] lazyeval_0.2.2              RCurl_1.98-1.17            
    ##  [9] GenomicAlignments_1.44.0    XML_3.99-0.18              
    ## [11] digest_0.6.37               lifecycle_1.0.4            
    ## [13] cluster_2.1.8.1             tidytree_0.4.6             
    ## [15] magrittr_2.0.3              compiler_4.5.1             
    ## [17] rlang_1.1.6                 tools_4.5.1                
    ## [19] yaml_2.3.10                 rtracklayer_1.68.0         
    ## [21] knitr_1.50                  S4Arrays_1.8.1             
    ## [23] curl_6.4.0                  DelayedArray_0.34.1        
    ## [25] RColorBrewer_1.1-3          aplot_0.2.8                
    ## [27] abind_1.4-8                 BiocParallel_1.42.1        
    ## [29] withr_3.0.2                 purrr_1.1.0                
    ## [31] BiocGenerics_0.54.0         ggh4x_0.3.1                
    ## [33] stats4_4.5.1                colorspace_2.1-1           
    ## [35] iterators_1.0.14            scales_1.4.0               
    ## [37] SummarizedExperiment_1.38.1 cli_3.6.5                  
    ## [39] rmarkdown_2.29              crayon_1.5.3               
    ## [41] treeio_1.32.0               generics_0.1.4             
    ## [43] rstudioapi_0.17.1           httr_1.4.7                 
    ## [45] rjson_0.2.23                ape_5.8-1                  
    ## [47] parallel_4.5.1              ggplotify_0.1.2            
    ## [49] XVector_0.48.0              restfulr_0.0.16            
    ## [51] matrixStats_1.5.0           yulab.utils_0.2.0          
    ## [53] vctrs_0.6.5                 Matrix_1.7-3               
    ## [55] jsonlite_2.0.0              GetoptLong_1.0.5           
    ## [57] gridGraphics_0.5-1          IRanges_2.42.0             
    ## [59] patchwork_1.3.1             S4Vectors_0.46.0           
    ## [61] ggrepel_0.9.6               clue_0.3-66                
    ## [63] foreach_1.5.2               tidyr_1.3.1                
    ## [65] glue_1.8.0                  codetools_0.2-20           
    ## [67] shape_1.4.6.1               gtable_0.3.6               
    ## [69] GenomeInfoDb_1.44.1         GenomicRanges_1.60.0       
    ## [71] BiocIO_1.18.0               UCSC.utils_1.4.0           
    ## [73] tibble_3.3.0                pillar_1.11.0              
    ## [75] htmltools_0.5.8.1           circlize_0.4.16            
    ## [77] GenomeInfoDbData_1.2.14     R6_2.6.1                   
    ## [79] ggpattern_1.1.4             doParallel_1.0.17          
    ## [81] evaluate_1.0.4              lattice_0.22-7             
    ## [83] Biobase_2.68.0              png_0.1-8                  
    ## [85] Rsamtools_2.24.0            ggfun_0.2.0                
    ## [87] Rcpp_1.1.0                  gridExtra_2.3              
    ## [89] SparseArray_1.8.1           nlme_3.1-168               
    ## [91] xfun_0.52                   GlobalOptions_0.1.2        
    ## [93] fs_1.6.6                    MatrixGenerics_1.20.0      
    ## [95] pkgconfig_2.0.3
