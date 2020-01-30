# MS-PEP
A shiny app to analyze mass spectrometry data for peptide cleavage patterns

## Requirements
R 3.5 or higher

Optionally: RStudio or other interface

## Required data
MS data in table form:

  Collumn 1: M/Z
  
  Collumn 2: Intensity

## How to use:

### First use:
Uncomment line 28-35 in app.R in case some packages still need to be installed.
(Optionally but not necesarry: comment line 39-45)

### All uses:
Download both app.R and mzR-functions.R files in single directory.
Run app.R using RStudio or by running following command in R:
shiny::runApp('Path to directory where saved')

Upload MS data and insert sequence. Click on submit.

## Packages used:
    MALDIquant
    MALDIquantForeign
    xlsx
    stringr
    shiny
    shinyWidgets
    plotly
