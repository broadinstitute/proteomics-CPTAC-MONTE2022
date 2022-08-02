library(dplyr)
library(cmapR)

# loaf discovery dataset from original data viewer
load("data/data-luad-v3.2-tumor-over-nat.RData")
disco.column.anno <- column.anno
disco.row.anno <- row.anno
disco_table <- tab.expr.all

## extracct MONTE data from mega_table
mega_table <- readRDS("data/MONTE_mega-table_2022-07-05.rds")

# select only the columns you want
mega_table <- select(mega_table, expt, ome, id, variableSites, sequence,
                     geneSymbol, ratio, participant, ratio, in_cls1, in_cls2)

# exclude genes with mutations (genes with no geneSymbol)
mega_table <- mega_table[!is.na(mega_table$geneSymbol), ]

# get only monte samples
monte_table <- subset(mega_table, expt %in% c("hlaft", "noip"))


## get HLA information
hla.table <- readRDS("data/hla_table_combined_filtered_2022-07-13.rds")

# small edit for HLA table
hla.table$directory <- sub('-', '.', hla.table$directory)

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

# read in annotation file
anno_table <- read.csv("data/luad-v3.2-sample-annotation-all.csv")

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

save(anno_table, disco_table, disco.column.anno, disco.row.anno, hla.table, monte_table, MONTE.HLA.counts,
     file = "data/MONTE_data_2022-08-02.RData")



# # load in mega table
# mega_table <- readRDS("MONTE_mega-table_2022-07-05.rds")
# 
# # select only the columns you want
# # can add more columns later
# mega_table <- select(mega_table, expt, ome, id, variableSites, sequence,
#                      geneSymbol, ratio, participant, ratio, in_cls1, in_cls2)
# 
# # exclude genes with mutations (genes with no geneSymbol)
# mega_table <- mega_table[!is.na(mega_table$geneSymbol), ]
# 
# # get annotations for HLA class presence
# anno_HLA <- data.frame(matrix(nrow=2, ncol=length(unique(mega_table$participant))))
# for (p in unique(mega_table$participant)) {
#   t <- filter(mega_table, participant == p)
#   stopifnot(all(t$in_cls1 == t$in_cls1[1]))
#   stopifnot(all(t$in_cls2 == t$in_cls2[1]))
#   anno_HLA[1, p] <- t$in_cls1
#   anno_HLA[2, p] <- t$in_cls2
# }




# # split into discovery and MONTE datasets
# mega_table$expt <- as.factor(mega_table$expt)
# disco_table <- subset(mega_table, expt=="disco")
# monte_table <- subset(mega_table, expt %in% c("hlaft", "noip"))
# 
# # get annotations for experiment type for monte
# anno_expt <- data.frame(matrix(nrow=1, ncol=length(unique(monte_table$participant))))
# colnames(anno_expt) <- unique(monte_table$participant)
# for (p in unique(monte_table$participant)) {
#   t <- filter(monte_table, participant == p)
#   print(paste(p, ": ", unique(t$expt)))
# }

# ## read in transcriptome data
# trans_gct_file_path <- "luad-v3.2-rnaseq-prot-uq-rpkm-log2-NArm.gct"
# trans_gct <- parse_gctx(trans_gct_file_path)
# 
# participants <- c("C3N.02145", "C3N.00199", "C3N.00169", "C3N.00547",
#                   "C3N.00579", "C3N.01016", "C3L.01632", "C3N.01024",
#                   "C3N.01416", "C3L.02549")
# 
# anno_table <- read.csv("luad-v3.2-sample-annotation-all.csv")
# 
## annotaion tracks shown in heatmap
# anno.keep <- c('NMF.consensus',
#                'mRNA.Expression.Subtype.TCGA',
#                'Stage',
#                'Smoking.Score.WGS',
#                'TSNet.Purity',
#                'ESTIMATE.ImmuneScore',
#                'CIMP.status',
#                'TP53.mutation.status',
#                "KRAS.mutation.status" ,
#                "EGFR.mutation.status",
#                "STK11.mutation.status",
#                "ALK.fusion")
# 
# anno_table <- select(anno_table, all_of(c("Sample.ID", anno.keep)))
# 
# save(monte_table, disco_table, anno_table, file="MONTE_data_2022-07-06.RData")

###################
# # original discovery dataset
# hla.table <- readRDS("data/hla_table_combined_filtered_2022-07-13.rds")

# load("data/MONTE_data_2022-07-06.RData")
# load("data/data-luad-v3.2-tumor-over-nat.RData")
# hla.table <- readRDS("data/hla_table_combined_filtered_2022-07-13.rds")
# 
# save(anno_table, disco_table, disco.column.anno, disco.row.anno, hla.table, monte_table,
#      file = "data/MONTE_data_2022-07-18.RData")
