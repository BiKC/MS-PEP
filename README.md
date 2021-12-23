# MS-PEP
A shiny app to analyze mass spectrometry data for peptide cleavage patterns

## Requirements
R 3.5 or higher

Optionally: RStudio or other interface

## Required data
1. MS data txt file in tab-separated table form:

  Column 1: M/Z
  
  Column 2: Intensity
  
2. Peptide sequence that was analyzed

## How to use:

### First use:
Uncomment line 28-35 in app.R in case some packages still need to be installed.
(Optionally but not necesarry: comment line 39-45)

### All uses:
Download both app.R and mzR-functions.R files in single directory.
Run app.R using RStudio or by running following command in R:
shiny::runApp('Path to directory where app.R is saved')

Upload MS data and insert sequence. Click on submit.

## Packages used:
    MALDIquant
    MALDIquantForeign
    xlsx
    stringr
    shiny
    shinyWidgets
    plotly

## Online tool:
The tool has been published online at:
https://bioit.shinyapps.io/PARs/

I added a line!