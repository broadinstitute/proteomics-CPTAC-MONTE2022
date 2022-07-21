###############################################################################
## Filename: server.R
## Created: July 5, 2022
## Author: Stephanie Vartany
## Purpose: Data viewer for "MONTE enables serial immunopeptidome, ubiquitylome,
## proteome, phosphoproteome, acetylome analyses of sample-limited tissues"
## This file controls the server. It mostly calls functions from shinyModules.R
## Note: code is heavily adapted from Karsten Krug's previous shiny apps
###############################################################################

shinyServer( function(input, output, session) {
  viewerTabServer("MonteTab", monte_table, monte.params)
  viewerTabServer("DiscoTab", disco_table, disco.params)
  HLATableServer("MonteTabHLA", hla.table, monte.params)
})