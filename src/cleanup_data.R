library(dplyr)
library(cmapR)

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
