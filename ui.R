###############################################################################
## Filename: ui.R
## Created: July 5, 2022
## Author: Stephanie Vartany
## Purpose: Data viewer for "MONTE enables serial immunopeptidome, ubiquitylome,
## proteome, phosphoproteome, acetylome analyses of sample-limited tissues"
## This file controls the user interface. It mostly calls functions from 
## shinyModules.R
## Note: code is heavily adapted from Karsten Krug's previous shiny apps
###############################################################################

## Define UI
shinyUI(fluidPage(
  
  ## Application title
  titlePanel(HTML(TITLESTRING), windowTitle = WINDOWTITLE), 
  
  tabsetPanel(
    ## MONTE dataset tab
    tabPanel("MONTE Dataset",
      viewerTabUI("MonteTab", "Monte Tab", monte.params),
      HLATableUI("MonteTabHLA", params = monte.params)),
    
    ## Discovery dataset tab
    tabPanel("Discovery Dataset",
             viewerTabUI("DiscoTab", "Discovery Tab", disco.params))
  )
))
