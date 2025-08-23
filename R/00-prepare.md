For this course, you will need an up-to-date version of R (4.4.0 or
higher) and RStudio. In the course of the workshop, we will be using a
number of packages. Let’s get them all installed before we start.

# Set up and installations

## BiocManager

``` r
if (!any(rownames(installed.packages()) == "BiocManager")){
  install.packages("BiocManager")
}
```

## Rmarkdown

``` r
if (!any(rownames(installed.packages()) == "Rmarkdown")){
  BiocManager::install("Rmarkdown")
}
```

# Cleaning and managing data

## lattice

``` r
if (!any(rownames(installed.packages()) == "dplyr")){
  install.packages("lattice")
}
```

## reshape2

``` r
if (!any(rownames(installed.packages()) == "dplyr")){
  install.packages("reshape2")
}
```

## dplyr

``` r
if (!any(rownames(installed.packages()) == "dplyr")){
  BiocManager::install("dplyr")
}
```

## tidyr

``` r
if (!any(rownames(installed.packages()) == "tidyr")){
  BiocManager::install("tidyr")
}
```

## magrittr

``` r
if (!any(rownames(installed.packages()) == "magrittr")){
  BiocManager::install("magrittr")
}
```

## kableExtra

``` r
if (!any(rownames(installed.packages()) == "kableExtra")){
  BiocManager::install("kableExtra")
}
```

# Basic visualizations

## ggplot2

``` r
if (!any(rownames(installed.packages()) == "ggplot2")){
  BiocManager::install("ggplot2")
}
```

## viridis

``` r
if (!any(rownames(installed.packages()) == "viridis")){
  BiocManager::install("viridis")
}
```

## ggExtra

``` r
if (!any(rownames(installed.packages()) == "ggExtra")){
  BiocManager::install("ggExtra")
}
```

## ggsignif

``` r
if (!any(rownames(installed.packages()) == "ggsignif")){
  BiocManager::install("ggsignif")
}
```

## ggrepel

``` r
if (!any(rownames(installed.packages()) == "ggrepel")){
  BiocManager::install("ggrepel")
}
```

## ggalluvial

``` r
if (!any(rownames(installed.packages()) == "ggalluvial")){
  BiocManager::install("ggalluvial")
}
```

# Advanced visualizations

## circlize

``` r
if (!any(rownames(installed.packages()) == "circlize")){
  BiocManager::install("circlize")
}
```

## Signac

``` r
if (!any(rownames(installed.packages()) == "Signac")){
  BiocManager::install("Signac")
}
```

## Seurat

``` r
if (!any(rownames(installed.packages()) == "Seurat")){
  BiocManager::install("Seurat")
}
```

## ggtree

``` r
if (!any(rownames(installed.packages()) == "ggtree")){
  BiocManager::install("ggtree")
}
```

## ggcoverage

``` r
if (!any(rownames(installed.packages()) == "ggcoverage")){
  install.package("remotes")
  remotes::install_github("showteeth/ggcoverage")
}
```

## Complex Heatmap

``` r
if (!any(rownames(installed.packages()) == "ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap")
}
```

# Download markdown documentation

``` r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2025-August-Intermediate-Visualization-for-Bioinformatics/R/01-boxplot.Rmd")
```
