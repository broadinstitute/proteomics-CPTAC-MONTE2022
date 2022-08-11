library(dplyr)
library(cmapR)

# loaf discovery dataset from original data viewer
load("original-data/data-luad-v3.2-tumor-over-nat.RData")
disco.column.anno <- column.anno
disco.row.anno <- row.anno
disco_table <- tab.expr.all

## extracct MONTE data from mega_table
mega_table <- readRDS("original-data/MONTE_mega-table_2022-07-05.rds")

# select only the columns you want
mega_table <- select(mega_table, expt, ome, id, variableSites, sequence,
                     geneSymbol, ratio, participant, ratio, in_cls1, in_cls2)

# exclude genes with mutations (genes with no geneSymbol)
mega_table <- mega_table[!is.na(mega_table$geneSymbol), ]

# get only monte samples
monte_table <- subset(mega_table, expt %in% c("hlaft", "noip"))


# ## get HLA information
# hla.table <- readRDS("data/hla_table_combined_filtered_2022-07-13.rds")
# 
# # small edit for HLA table
# hla.table$directory <- sub('-', '.', hla.table$directory)
# 
# # get MONTE HLA counts
# MONTE.HLA.counts <- data.frame(matrix(nrow=2*length(unique(monte_table$geneSymbol)),
#                                            ncol = 2 + length(unique(monte_table$participant))))
# names(MONTE.HLA.counts) <- c("geneSymbol", "HLA.class", sort(unique(monte_table$participant)))
# MONTE.HLA.counts$geneSymbol <- rep(sort(unique(monte_table$geneSymbol)), each=2)
# MONTE.HLA.counts$HLA.class <- rep(c(1,2), nrow(MONTE.HLA.counts)/2)
# 
# for (gene in unique(MONTE.HLA.counts$geneSymbol)) {
#   t_cls1 <- filter(hla.table, geneSymbol == gene, hla_class == "cls1")
#   t_cls2 <- filter(hla.table, geneSymbol == gene, hla_class == "cls2")
#   
#   gene_row_idx <- sort(which(MONTE.HLA.counts$geneSymbol == gene)) [1]
#   
#   # deal with class 1
#   if (nrow(t_cls1) == 0) {
#     MONTE.HLA.counts[gene_row_idx, -(1:2)] <- 0
#   } else {
#     for (participant in unique(monte_table$participant)) {
#       MONTE.HLA.counts[gene_row_idx, participant] <- 
#         sum(t_cls1$directory == participant)
#     }
#   }
# 
#   # deal with class 2
#   if (nrow(t_cls2) == 0) {
#     MONTE.HLA.counts[gene_row_idx+1, -(1:2)] <- 0
#   } else {
#     for (participant in unique(monte_table$participant)) {
#       MONTE.HLA.counts[gene_row_idx+1, participant] <- 
#         sum(t_cls2$directory == participant)
#     }
#   }
# }

# read in annotation file
anno_table <- read.csv("original-data/luad-v3.2-sample-annotation-all.csv")

anno.keep <- c('NMF.consensus',
               'mRNA.Expression.Subtype.TCGA',
               'Stage',
               'Smoking.Score.WGS',
               'TSNet.Purity',
               'ESTIMATE.ImmuneScore',
               'CIMP.status',
               'TP53.mutation.status',
               "KRAS.mutation.status" ,
               "EGFR.mutation.status",
               "STK11.mutation.status",
               "ALK.fusion")

anno_table <- select(anno_table, all_of(c("Sample.ID", anno.keep)))


###############################################################
## rearrange new HLA tables

# read in tables
hla_cls1_table <- read.csv('original-data/HLA_I_Combined_ForViewer.txt', sep='\t')
hla_cls2_table <- read.csv('original-data/HLA_II_Combined_ForViewer.txt', sep='\t')

# select the columns you want
hla_cls1_table <- select(hla_cls1_table, c('directory','geneSymbol','sequence','PTM_nuORF'))
hla_cls2_table <- select(hla_cls2_table, c('directory','geneSymbol','sequence','PTM_nuORF'))

# add in the "class" column
hla_cls1_table$hla_class <- rep('cls1', dim(hla_cls1_table)[1])
hla_cls2_table$hla_class <- rep('cls2', dim(hla_cls2_table)[1])

# small chance in 'directory' syntax
hla_cls1_table$directory <- sub('-', '.', hla_cls1_table$directory)
hla_cls2_table$directory <- sub('-', '.', hla_cls2_table$directory)

# change columns around so hla_class is first
hla_cls1_table <- hla_cls1_table[, c(5, 1:4)]
hla_cls2_table <- hla_cls2_table[, c(5, 1:4)]

hla.table <- rbind(hla_cls1_table, hla_cls2_table)

# get MONTE HLA counts
MONTE.HLA.counts <- data.frame(matrix(nrow=2*length(unique(monte_table$geneSymbol)),
                                      ncol = 2 + length(unique(monte_table$participant))))
names(MONTE.HLA.counts) <- c("geneSymbol", "HLA.class", sort(unique(monte_table$participant)))
MONTE.HLA.counts$geneSymbol <- rep(sort(unique(monte_table$geneSymbol)), each=2)
MONTE.HLA.counts$HLA.class <- rep(c(1,2), nrow(MONTE.HLA.counts)/2)

for (gene in unique(MONTE.HLA.counts$geneSymbol)) {
  t_cls1 <- filter(hla.table, geneSymbol == gene, hla_class == "cls1")
  t_cls2 <- filter(hla.table, geneSymbol == gene, hla_class == "cls2")
  
  gene_row_idx <- sort(which(MONTE.HLA.counts$geneSymbol == gene)) [1]
  
  # deal with class 1
  if (nrow(t_cls1) == 0) {
    MONTE.HLA.counts[gene_row_idx, -(1:2)] <- 0
  } else {
    for (participant in unique(monte_table$participant)) {
      MONTE.HLA.counts[gene_row_idx, participant] <- 
        sum(t_cls1$directory == participant)
    }
  }
  
  # deal with class 2
  if (nrow(t_cls2) == 0) {
    MONTE.HLA.counts[gene_row_idx+1, -(1:2)] <- 0
  } else {
    for (participant in unique(monte_table$participant)) {
      MONTE.HLA.counts[gene_row_idx+1, participant] <- 
        sum(t_cls2$directory == participant)
    }
  }
}


########################

save(anno_table, disco_table, disco.column.anno, disco.row.anno, hla.table, monte_table, MONTE.HLA.counts,
     file = "data/MONTE_data_2022-08-11.RData")




