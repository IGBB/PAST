---
layout: page
permalink: /running_rshiny/
title: Running PAST in Shiny
---

## Installing

PAST Shiny is available through the Shiny branch of our Github repository. The Shiny application can be downloaded <a href="https://github.com/IGBB/PAST/releases/download/v1.0.0/app.R">here</a> or from the [Releases page](https://github.com/IGBB/PAST/releases/). Once the Shiny application has been downloaded, some dependencies need to be installed. In an R Console, run the following line of code to install PAST and its dependencies as well as the dependencies of the PAST Shiny.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PAST")
install.packages(c("DT", "zip", "shiny", "shinydashboard", "gridExtra"))
```

## Running PAST Shiny

PAST Shiny can be run easily with two methods. First, PAST Shiny can be run through the R Console with the following line of R code. `path/to/folder/with/shiny/` should be replaced with the real path to the downloaded application.

```r
shiny::runApp('path/to/folder/with/shiny/app.R')
```

PAST Shiny can also be run through RStudio. If `app.R` is opened in RStudio, RStudio detects that PAST Shiny is an R Shiny application and provides a "Run App" button at the top of the editor. PAST Shiny works best with RStudio 1.2 or higher, as an error in RStudio's browser prevents the results from being downloaded as a ZIP archive.


