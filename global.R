###############################################################################
## Filename: global.R
## Created: July 5, 2022
## Author: Stephanie Vartany
## Purpose: Data viewer for "Workflow enabling deepscale immunopeptidome, 
## proteome, ubiquitylome, phosphoproteome, and acetylome analyses of 
## sample-limited tissues"
## This file imports the data, completes necessary preprocessing, and defines
## global parameters. It also defines functions from src/shinyModules.R and
## src/helperFunctions.R
## Note: code is adapted from Karsten Krug's previous shiny apps
###############################################################################

###############################################################################
# # Use this when publishing to RStudio Connect
# options(repos = c(options('repos')$repos,
#                   bioc = 'https://bioconductor.org/packages/3.15/bioc/'))
###############################################################################


## import libraries
library(shiny)
library(dplyr)
library(circlize)
library(ggplot2)
library(WriteXLS)
library(ComplexHeatmap)
library(glue)
library(reactable)
library(tippy)


source("src/shinyModules.R")
source("src/helperFunctions.R")

## import the data
load("data/MONTE_data_2022-11-28.RData")

## global parameters
TITLESTRING <<- '<font size="5" face="times"><i><b>"Workflow enabling deepscale immunopeptidome, proteome, ubiquitylome, phosphoproteome, and acetylome analyses of sample-limited tissues"</b></i> (<a href="https://www.biorxiv.org/content/10.1101/2021.06.22.449417v1" target="_blank_">citation here</a>)</font><br>'

WINDOWTITLE <<- 'MONTE-CPTAC-LUAD-2022' # TODO: edit this
FILENAMESTRING <<- 'MONTE-CPTAC-LUAD-2022' # TODO: edit this

GENEMAX <<- 20

## other parameters
monte_participants <<- unique(monte_table$participant)

#########################################
## dataset preprocessing

# convert HLA counts to binary
MONTE.HLA.counts[, -(1:2)] <- as.integer(MONTE.HLA.counts[,-(1:2)] > 0)

##########################################
## annotaion tracks shown in heatmap 
anno.all <- rev(c('Multi.omic.subtype'='NMF.consensus', 
                  'RNA.subtype.TCGA'='mRNA.Expression.Subtype.TCGA',
                  'Tumor.Stage'='Stage',
                  'Smoking.Score.WGS'='Smoking.Score.WGS',
                  'TSNet.Purity'='TSNet.Purity',
                  'ESTIMATE.ImmuneScore'='ESTIMATE.ImmuneScore',
                  'CIMP.status'='CIMP.status',
                  'TP53.mutation'='TP53.mutation.status',
                  "KRAS.mutation"="KRAS.mutation.status" ,
                  "EGFR.mutation"="EGFR.mutation.status",
                  "STK11.mutation"="STK11.mutation.status",
                  "ALK.fusion"="ALK.fusion"
))

###########################################
## dataset-specific parameters
disco.params <- list( #oms = c("Transcriptome", "Ubiquitylome", "Proteome", "Phosphoproteome", "Acetylome"),
                     genes.all = sort(unique(disco.row.anno$geneSymbol)),
                     genes.start = c('EGFR', 'RB1', 'KRAS', 'STK11'),
                     annotations = c("In.MONTE", names(rev(anno.all))),
                     annotations.start = "In.MONTE")
monte.params <- list( #oms = c("Immunopeptidome", "Ubiquitylome", "Proteome", "Phosphoproteome", "Acetylome"),
                     # genes.all = sort(unique(monte_table$geneSymbol)),
                     genes.all = sort(unique(c(disco.row.anno$geneSymbol, monte_table$geneSymbol))),
                     genes.HLA = sort(unique(c(hla.table$geneSymbol, monte_table$geneSymbol))),
                     genes.start = c('EGFR', 'RB1', 'KRAS', 'STK11'),
                     annotations = c("Experiment", names(rev(anno.all))),
                     annotations.start = "Experiment")

##############################
## dataset descriptions
disco.params$data_description <- glue( '<p>Copy number aberrations are relative to matching normal blood sample and are on log2(CNA)-1 scale. For other data types the heatmap depicts abundances observed in tumor relative to normal adjacent tissue (NAT).</p>\n
                   <table style="width:80%"><tr><th>Type</th><th># features</th><th># samples</th></tr>\n
                   <tr><td>CNA</td><td>{sum(grepl("_CNA", disco.row.anno$DataType))}</td><td>{ncol(disco_table)}</td></tr>\n
                   <tr><td>mRNA</td><td>{sum(grepl("_RNAseq", disco.row.anno$DataType))}</td><td>{ncol(disco_table)}</td></tr>\n
                   <tr><td>Protein</td><td>{sum(grepl("_Protein", disco.row.anno$DataType))}</td><td>{ncol(disco_table)}</td></tr>\n
                   <tr><td>Phosphosites </td><td>{sum(grepl("_pSTY", disco.row.anno$DataType))}</td><td>{ncol(disco_table)}</td></tr>\n
                   <tr><td>Acetylsites</td><td>{sum(grepl("_acK", disco.row.anno$DataType))}</td><td>{ncol(disco_table)}</td></tr>\n
                   </table>\n
                   <br><p>For more details please see our reference publication <a href="https://www.cell.com/cell/fulltext/S0092-8674(20)30744-3" target="_blank_">Gillette <i>et al.</i> Cell, 2020</a></p>')
monte.params$data_description <- glue('<p>The heatmap depicts multi-ome data from the 
                   MONTE workflow ("hlaft") and the serial multi-omic enrichment workflow ("noip"). 
                   CNA and RNAseq are copied from the LUAD 2020 discovery dataset.
                   HLA tracks show whether at least one class 1 or class 2 peptide was detected in the MONTE workflow.
                   For other data types the heatmap depicts abundances observed in tumor relative to normal adjacent tissue (NAT).<p>\n
                   <table style="width:100%"><tr><th>Type</th><th># MONTE  </th><th># serial emrichment  </th><th># samples</th></tr>\n
                   <tr><td>Protein </td><td>{sum(monte_table$expt == "hlaft" & monte_table$ome == "prot")}</td>
                                       <td>{sum(monte_table$expt == "noip" & monte_table$ome == "prot")}
                                       </td><td>{length(monte_participants)}</td></tr>\n
                   <tr><td>Phosphosites </td><td>{sum(monte_table$expt == "hlaft" & monte_table$ome == "phos")}</td>
                                        <td>{sum(monte_table$expt == "noip" & monte_table$ome == "phos")}</td>
                                        <td>{length(monte_participants)}</td></tr>\n
                   <tr><td>Acetylsites </td><td>{sum(monte_table$expt == "hlaft" & monte_table$ome == "acet")}</td>
                                        <td>{sum(monte_table$expt == "noip" & monte_table$ome == "acet")}</td>
                                        <td>{length(monte_participants)}</td></tr>\n
                   <tr><td>Ubiquitylsites </td><td>{sum(monte_table$expt == "hlaft" & monte_table$ome == "ubiq")}</td>
                                        <td>{sum(monte_table$expt == "noip" & monte_table$ome == "ubiq")}</td>
                                        <td>{length(monte_participants)}</td></tr>\n
                   </table>')


##############################
## color mappings for 'anno.all'
column.anno.col <<- list(
  Multi.omic.subtype=c(C1=rgb(142, 210, 198, maxColorValue = 255), 
                       C2=rgb(250, 247, 182, maxColorValue = 255), 
                       C3=rgb(190, 187, 219, maxColorValue = 255), 
                       C4=rgb(244, 127, 114, maxColorValue = 255)),
  RNA.subtype.TCGA=c( 'Proximal-proliferative'=rgb(130, 178, 212, maxColorValue = 255),
                      'Proximal-inflammatory'=rgb(187, 129, 184, maxColorValue = 255),
                      'Terminal Respiratory Unit'=rgb(252, 180, 98, maxColorValue = 255)),
  CIMP.status=c('CIMP-1'=rgb(179, 205, 227, maxColorValue = 255),
                'CIMP-2'=rgb(139, 150, 199, maxColorValue = 255),
                'CIMP+'=rgb(136, 69, 153, maxColorValue = 255)),
  Tumor.Stage=c('1'=blues9[2], '1A'=blues9[3], '1B'=blues9[4], 
                '2A'=blues9[5], '2B'=blues9[6], '3'=blues9[7], '3A'=blues9[8],
                'NA' = 'grey90'),
  TP53.mutation=c('0'='grey90', '1'='black'),
  KRAS.mutation=c('0'='grey90', '1'='black'),
  EGFR.mutation=c('0'='grey90', '1'='black'),
  STK11.mutation=c('0'='grey90', '1'='black'),
  ALK.fusion=c('0'='grey90', '1'='darkblue'),
  
  Smoking.Score.WGS=circlize::colorRamp2( c(0, 0.5, 1), c('grey90', 'grey50','black')),
  
  TSNet.Purity=circlize::colorRamp2( c(0, 0.5, 1), c('grey90', 'grey50','black')),
  ESTIMATE.ImmuneScore=circlize::colorRamp2( c(676, 4700, 6800), 
                                             c(rgb(255,245,240, maxColorValue = 255),  
                                               rgb(251,106,74, maxColorValue = 255),  
                                               rgb(165,21,22, maxColorValue = 255)))
  
)
