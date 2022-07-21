params <- list(
  genes.char = "MACF1, BRDT, SLK, PIK3CA, PIK3CD, SOX2, TP63, VIM",
  zscore = 'row',
  PTMsites = 'most variable',
  min.val = -3,
  max.val = 3,
  sort.after = 'Multi.omic.subtype',
  id = "MonteTab")

table <- monte_table

HM.out <- myComplexHeatmap(table, params)


########

HM.Table <- HM.out$Table
HM <- HM.out$HM

geneSymbols <- HM.Table[!is.na(HM.Table$geneSymbol), ]$geneSymbol


# lets just make the first gene for now
gene <- geneSymbols[1]
HM.subset <- HM[geneSymbols == gene]

# get HLA class data for this gene
hla.classes.table <- filter(table, geneSymbol == gene)



###########

hla.table <- readRDS("data/hla_table_combined_filtered_2022-07-13.rds")
hla.table$participant <- sub('-','.', hla.table$directory)

n.genes <- length(unique(hla.table$geneSymbol))
n.participants <- length(unique(hla.table$participant))

HLA.classes <- data.frame(matrix(nrow = n.genes*2, ncol = 2 + n.participants))
names(HLA.classes) <- c('geneSymbol', 'HLA.Class', sort(unique(hla.table$participant)))
HLA.classes$geneSymbol <- rep(sort(unique(hla.table$geneSymbol)), each=2)
HLA.classes$HLA.Class <- rep(c(1, 2), times = nrow(HLA.classes)/2)

for (gene in sort(unique(hla.table$geneSymbol))) {
  for (part in unique(hla.table$participant)) {
    t <- filter(hla.table, geneSymbol == gene & participant == part)
    row_idx <- which(HLA.classes$geneSymbol == gene & HLA.classes$HLA.Class == 1)
    HLA.classes[row_idx, part] <- sum(t$hla_class == 'cls1')
    HLA.classes[row_idx+1, part] <- sum(t$hla_class == 'cls2')
  }
  print(paste(gene, ' is done!'))
}