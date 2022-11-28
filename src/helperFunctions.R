###############################################################################
## Filename: helperFunctions.R
## Created: July 6, 2022
## Author: Stephanie Vartany
## Purpose: Data viewer for "MONTE enables serial immunopeptidome, ubiquitylome,
## proteome, phosphoproteome, acetylome analyses of sample-limited tissues"
## This file contains helper functions to visualize heatmaps and HLA tables
## Note: code is adapted from Karsten Krug's previous shiny apps
###############################################################################

## function to extract gene names from a string input
extractGenes <- function(genes.char, table){
  
  if(is.null(genes.char))
    return(NULL)
  
  if( length(genes.char) == 0 ){
    return(NULL)
  }
  ## extract genes
  genes.vec= unlist(strsplit(genes.char, ','))
  if(length(genes.vec)==1)
    genes.vec=unlist(strsplit(genes.char, ' '))
  if(length(genes.vec)==1)
    genes.vec=unlist(strsplit(genes.char, ';'))
  
  ## unique gene names
  genes.vec <- unique(genes.vec)
  
  ## limit to 'GENEMAX' genes
  if(length(genes.vec) > GENEMAX){
    warning(paste('more than', GENEMAX,'gene ids submitted! Showing results for the first 20 genes in the list.\n'))
    genes.vec <- genes.vec[1:GENEMAX]
  }
  
  ## remove spaces
  genes.vec <- sapply(genes.vec, function(x) gsub(' ', '', x))
  
  ## exclude genes that are not in the dataset
  genes.notInTable <- setdiff(genes.vec, table$geneSymbol)
  if (length(genes.notInTable) > 0) {
    warning("Some genes were not found in the dataset. These are being excluded.\n")
    genes.vec <- setdiff(genes.vec, genes.notInTable)
  }
  
  return(list(genes.vec=genes.vec, genes.notInTable=genes.notInTable))
}

## make row label
makeRowLabel <- function(gene.table) {
  VMsites <- gene.table$variableSites
  VMsites[is.na(VMsites)] <- ''
  gene.table$row_label <- paste(paste(toupper(substr(gene.table$ome,1,1)),
                                      substr(gene.table$ome, 
                                             start=2, 
                                             stop=nchar(gene.table$ome)),
                                      sep=''), 
                                sapply(gene.table$id, function(x) strsplit(x,'_')[[1]][1]),
                                VMsites,
                                sep = ' ')
  return(gene.table)
}

## function to create a heatmap matrix for a single gene
getHMTable <- function(gene, table, params) {
  
  zscore <- params$zscore
  PTMsites <- params$PTMsites
  
  # extract rows for that gene
  gene.table <- filter(table, geneSymbol == gene)
  
  # come up with correct row labels
  gene.table <- makeRowLabel(gene.table)
  
  # make dataframe with the entries for the heatmap
  row_label.unique <- sort(unique(gene.table$row_label))
  col_labels.unique <- sort(unique(table$col_label))
  HM.table <- data.frame(matrix(nrow = length(row_label.unique), 
                                ncol = length(col_labels.unique) + 3))
  names(HM.table) <- c("geneSymbol", "ome", "row_label", col_labels.unique)
  HM.table$row_label <- row_label.unique
  
  # find the correct spot in the data frame and copy the data over
  for (i in 1:dim(gene.table)[1]) {
    df.row <- gene.table[i,]
    col_idx <- which(names(HM.table) == df.row$col_label)
    row_idx <- which(HM.table$row_label == df.row$row_label)
    stopifnot(length(row_idx) == 1, length(col_idx) == 1)
    HM.table[row_idx, col_idx] <- df.row$ratio
    HM.table$geneSymbol[row_idx] <- df.row$geneSymbol
    HM.table$ome[row_idx] <- df.row$ome
  }
  
  # TODO: pick only one protein abundance to show 
  
  return (HM.table)
}

## extract most variable PTM
extractMostVariable <- function(HM.table) {
  
  omes <- c("phos", "acet", "ubiq")
  for (ome in omes) {
    for (gene in unique(HM.table$geneSymbol)) {
      ## extract MS phospho
      vm.idx <- which(HM.table$ome == ome & HM.table$geneSymbol == gene)
      if( length(vm.idx) > 1 ){
        
        vm.sd <- apply(HM.table[vm.idx, -(1:3)], 1, sd, na.rm=T)
        max.idx <- vm.idx[which.max(vm.sd)]
        rm.idx <- setdiff(vm.idx, max.idx)
        
        HM.table <- HM.table[-(rm.idx), ]
      }
    }
  }
  
  return (HM.table)
}

getHMTableDisco <- function(genes.Table, row.anno, params) {
  # make row labels
  row_labels <- paste(row.anno[ , 'geneSymbol'],
                      gsub( '5_acK', 'acK',
                            gsub( '4_pSTY', 'pSTY', 
                                  gsub('3_Protein', 'Protein', 
                                       gsub('2_RNAseq', 'RNA-Seq', 
                                            gsub('1_CNA', 'CNA', 
                                                 row.anno[ , 'DataType'])
                                       ) 
                                  )
                            )
                      )
  )
  
  ## MS
  ms.psty.idx <- grep('pSTY', row_labels)
  if(length(ms.psty.idx)>0){
    row_labels[ ms.psty.idx ] <- paste( sub('pSTY', '',
                                            row_labels[ ms.psty.idx ]),
                                        paste('p', 
                                              sub('.*_([S|T|Y][0-9]*)[s|t|y].*', 
                                                  '\\1', 
                                                  row.anno[ ms.psty.idx, 'ID']), 
                                              sep=''),
                                        sep='' )
  }
  ms.ack.idx <- grep('acK', row_labels)
  if(length(ms.ack.idx)>0){
    row_labels[ ms.ack.idx ] <- paste( sub('acK', '',
                                           row_labels[ ms.ack.idx ]),
                                       paste('ac', 
                                             sub('.*_([K][0-9]*)[k].*',
                                                 '\\1', 
                                                 row.anno[ ms.ack.idx, 'ID']), 
                                             sep=''), 
                                       sep='' )
  }
  # get omes
  ome <- row.anno$DataType
  ome[ome == '1_CNA'] <- 'CNA'
  ome[ome == '2_RNAseq'] <- 'RNAseq'
  ome[ome == '3_Protein'] <- 'prot'
  ome[ome == '4_pSTY'] <- 'phos'
  ome[ome == '5_acK'] <- 'acet'
  
  # combine into one genes.Table
  genes.Table <- cbind(data.frame(geneSymbol = row.anno$geneSymbol,
                                  ome = ome,
                                  row_label = row_labels),
                       genes.Table)
  rownames(genes.Table) <- make.unique(genes.Table$row_label, sep= ' ')
  
  return (genes.Table)
}

# function to add in HLA data to genes table for heatmap
addHLAToTable <- function(genes.Table, genes.vec, params) {
  hla.table <- filter(MONTE.HLA.counts, geneSymbol %in% genes.vec)
  participant.names <- names(hla.table[, 3:12])
  names(hla.table) <- c('geneSymbol','HLA.Class', 
                        paste(participant.names, ': hlaft', sep = ''))
  hla.table[paste(participant.names, ': noip', sep = '')] <- NA
  hla.table$ome <- paste('HLA.cls', hla.table$HLA.Class, sep='')
  hla.table$row_label <- paste('HLA class', hla.table$HLA.Class, sep = ' ')
  hla.table <- hla.table[sort(names(hla.table))]
  hla.table <- hla.table[, c(21, 23, 24, 1:20)]
  genes.Table <- rbind(genes.Table, hla.table)
  
  return(genes.Table)
}

## function to make a complex heatmap
myComplexHeatmap <- function(table, params) {
  # extract parameters
  genes.char <- params$genes.char
  zscore <- params$zscore
  PTMsites <- params$PTMsites
  min.val <- params$min.val
  max.val <- params$max.val
  sort.after <- params$sort.after
  dataset_id <- params$id
  
  ## preprocessing if its the MONTE dataset
  if (dataset_id == "MonteTab") {
    # general parameters
    num.participants <- length(unique(table$participant))
    participants.unique <- unique(table$participant)
    
    # change participant names to include experiment
    table <- mutate(table, col_label = paste(participant, expt, sep=": "))
    params$expt <- unique(table$expt)
    
    # extract genes
    genes.all <- extractGenes(genes.char, select(table, geneSymbol))
    genes.vec <- genes.all$genes.vec
    genes.notInMONTE <- genes.all$genes.notInTable
    
    # make table of NAs for genes.notInMONTE (to show that they're not in MONTE data)
    if (length(genes.notInMONTE) > 0) {
      genes.notInMONTE.Table <- data.frame(matrix(nrow=length(genes.notInMONTE),
                                                  ncol = 3 + num.participants*2))
      names(genes.notInMONTE.Table) <- c("geneSymbol", "ome", "row_label", sort(unique(table$col_label)))
      genes.notInMONTE.Table$geneSymbol <- genes.notInMONTE
      genes.notInMONTE.Table$row_label <- paste(genes.notInMONTE.Table$geneSymbol,
                                                "no data avaliable",
                                                sep = ' ')
    }
    
    # make genes.Table for the genes actually in the dataset
    if (length(genes.vec) > 0){
      # extract rows for selected genes
      table <- filter(table, geneSymbol %in% genes.vec)
      
      # make gene.Table and gene.Matrix for input to heatmap
      genes.Table <- lapply(genes.vec, function(x) getHMTable(x, table, params))
      if (is.list(genes.Table)) {
        genes.Table <- do.call(rbind, genes.Table)
      }
      genes.Table <- genes.Table[rowSums(is.na(genes.Table)) != ncol(genes.Table), ]
      
      # add in HLA data
      genes.Table <- addHLAToTable(genes.Table, genes.vec, params)
      
      # add RNA-seq and CNA by extracting from disco_table
      monte.participants.in.disco <- disco_table[, sort(participants.unique)]
      disco.genes.Table <- monte.participants.in.disco[which(disco.row.anno$geneSymbol %in% genes.vec), ]
      row.anno <- disco.row.anno[which(disco.row.anno$geneSymbol %in% genes.vec), ]
      disco.genes.Table <- getHMTableDisco(disco.genes.Table, row.anno, params) 
      disco.genes.Table <- filter(disco.genes.Table, ome %in% c('CNA', 'RNAseq'))
      disco.genes.Table <- cbind(disco.genes.Table, disco.genes.Table[, -(1:3)])
      names(disco.genes.Table) <- c(names(disco.genes.Table)[1:3],
                                    paste(names(disco.genes.Table)[4:13], ': hlaft', sep=''),
                                    paste(names(disco.genes.Table)[4:13], ': noip', sep=''))
      column.order <- names(genes.Table)
      disco.genes.Table <- disco.genes.Table[, column.order]
      disco.genes.Table$row_label <- disco.genes.Table$ome
      genes.Table <- rbind(genes.Table, disco.genes.Table)
      
      # order by gene symbol
      genes.Table <- genes.Table[order(genes.Table$geneSymbol), ]
      
      # use only most variable VM sites (if requested)
      if (params$PTMsites == "most variable") {
        genes.Table <- extractMostVariable(genes.Table)
      }
      
    }
      
    # combine genes.Table with genes.notInMONTE.Table if they both exist
    if (length(genes.vec) > 0 & length(genes.notInMONTE) > 0) {
      genes.Table <- rbind(genes.Table, genes.notInMONTE.Table)
    } else if (length(genes.notInMONTE) > 0) {
      genes.Table <- genes.notInMONTE.Table
    }
    
  }
  
  ## preprocessing if discovery dataset
  else if (dataset_id == "DiscoTab") {
    participants.unique <- unique(disco.column.anno$Sample.ID)
    num.participants <- length(participants.unique)
    params$expt <- "disco"
    
    # extract genes
    genes.all <- extractGenes(genes.char, select(disco.row.anno, geneSymbol))
    genes.vec <- genes.all$genes.vec
    
    # extract rows for that gene
    genes.Table <- table[which(disco.row.anno$geneSymbol %in% genes.vec), ]
    row.anno <- disco.row.anno[which(disco.row.anno$geneSymbol %in% genes.vec), ]
    
    genes.Table <- getHMTableDisco(genes.Table, row.anno, params) 
    genes.Table <- genes.Table[,c(1:3, order(names(genes.Table)[-(1:3)])+3)]
    
    # use only most variable VM sites (if requested)
    if (params$PTMsites == "most variable") {
      genes.Table <- extractMostVariable(genes.Table)
    }
  }
  
  else {
    stop("Invalid dataset parameter")
  }
  
  # make levels for omes (for heatmap ordering)
  genes.Table$ome <- factor(genes.Table$ome,
                 levels = c('CNA', 'RNAseq','prot','phos','ubiq','acet','HLA.cls1','HLA.cls2'))
  
  # make matrix to feed into heatmap
  genes.Matrix <- as.matrix(genes.Table[,-(1:3)])
  rownames(genes.Matrix) <- genes.Table$row_label
  
  # z-score rows if requested
  if (params$zscore == "row") {
    exclude.idx <- which(genes.Table$ome %in% c('CNA', 'HLA.cls1', 'HLA.cls2')) #exclude CNA and HLA from z-scoring
    genes.Matrix[-(exclude.idx),] <- t(apply(genes.Table[-(exclude.idx),-(1:3)], 1, 
                                              function(x)(x-mean(x, na.rm=T))/sd(x, na.rm=T)))
  }
  
  # define color map
  col_fun_ratios <- colorRamp2(c(min.val, 0, max.val), 
                               c("blue", "white", "red"))
  
  # make heatmap for annotations
  anno.fig <- getHeatmapAnnotations(participants.unique,
                                   anno_table = anno_table,
                                   anno.all = anno.all,
                                   params = params)
  if (params$id == "MonteTab") {
    column.anno.col$Experiment <- c("darkolivegreen2", "grey90")
    names(column.anno.col$Experiment) <- params$expt
  } else {
    column.anno.col$In.MONTE <- c('Yes' = "darkolivegreen2", 'No' = "grey90")
  }
  HM.anno <- HeatmapAnnotation(df = anno.fig,
                               col = column.anno.col,
                               show_legend = T,
                               show_annotation_name = T,
                               annotation_name_side = "left",
                               na_col = "grey90",
                               annotation_legend_param = list(
                                 direction = 'horizontal',
                                 max_width = 350),
                               height = unit(0.5, 'cm') * nrow(anno.fig))
  column.to.sort <- anno.fig[, which(colnames(anno.fig) == sort.after)]
  
  # make final heatmap for genes
  HM <- Heatmap(genes.Matrix,
                col = col_fun_ratios,
                column_title = "Participant",
                row_title_rot = 0,
                row_order = order(genes.Table$ome, genes.Table$row_label),
                cluster_rows = F,
                cluster_columns = F,
                row_split = genes.Table$geneSymbol,
                row_names_max_width = unit(8, 'cm'),
                column_split = column.to.sort, 
                top_annotation = HM.anno,
                show_row_names = T,
                show_column_names = T,
                name = 'relative abundance',
                heatmap_legend_param = list(direction='horizontal',
                                            max_width = 350,
                                            legend_width = unit(4, 'cm')),
                height = unit(0.5, 'cm') * nrow(genes.Matrix),
                column_names_side = "top",
                column_names_rot = ifelse(params$id == "MonteTab", 45, 90),
                column_names_gp = gpar(fontsize = ifelse(params$id == "MonteTab", 12, 8)),
                column_gap = if (is.double(column.to.sort)) {unit(0, "mm")} else {unit(1, "mm")})
  
  # combine annotation table with genes.Table for final output
  rownames(anno.fig) <- colnames(genes.Matrix)
  anno.fig.new <- t(anno.fig)
  names(anno.fig.new) <- colnames(anno.fig.new)
  
  anno.fig.new <- data.frame(matrix(ncol=ncol(genes.Table),
                                    nrow=ncol(anno.fig)))
  names(anno.fig.new) <- names(genes.Table)
  rownames(anno.fig.new) <- names(anno.fig)
  anno.fig.new[, -(1:3)] <- t(anno.fig)
  final.Table <- rbind(anno.fig.new, genes.Table)
  
  # make extra hla legend
  if (params$id == 'MonteTab') {
    hla.legend <- list(Legend(title = "HLA",
                              labels=c("detected", "not detected", "NA"),
                              legend_gp = gpar(fill = c(col_fun_ratios(c(1,0)), "grey")),
                              at = c('detected','not detected','NA')))
  } else {
    hla.legend = list()
  }
  
  return(list(HM=HM, 
              Table=final.Table, 
              hla.legend = hla.legend, 
              complexLegend = makeComplexLegend(column.anno.col)))
}

## function to dynamically determine the height (in px) of the heatmap
## depending on the number of genes
dynamicHeightHM <- function(n.entries){
  height <- 0.3*(n.entries+4) + 3 ## height in inch
  height <- height * 48             ## inch  to pixel
  
  return(height)
}

## function to generate heatmap annotations
getHeatmapAnnotations <- function(participants, anno_table, anno.all, params) {
  
  # filter for just the samples you want
  anno <- filter(anno_table, Sample.ID %in% participants)
  
  # edit column names
  colnames(anno) <- c("Sample.ID", names(rev(anno.all)))
  
  # deal with more than one experiment
  if (params$id == "MonteTab") {
    anno_new <- data.frame(matrix(nrow = 0, ncol = length(names(anno))))
    for (expt in params$expt) {
      df <- mutate(anno, Sample.ID = paste(Sample.ID, expt, sep = ": "), Experiment = expt)
      anno_new <- rbind(anno_new, df)
    }
    col_order <- c("Sample.ID", "Experiment", names(rev(anno.all)))
    anno <- anno_new[ , col_order]
  } else {
    In.MONTE <- anno$Sample.ID %in% monte_participants
    anno_new <- rep('No', length(In.MONTE))
    anno_new[In.MONTE] <- 'Yes'
    anno$In.MONTE <- anno_new
    col_order <- c("Sample.ID", "In.MONTE", names(rev(anno.all)))
    anno <- anno[ , col_order]
  }
  
  # make sure annotations dataframe is ordered the same as the heatmap
  anno <- anno[order(anno$Sample.ID), ] 
  rownames(anno) <- anno$Sample.ID
  anno.fig <- anno[, -(1)]
  
  return(anno.fig)
}

# ## function to get valid row labels for a single gene
# getGeneRowLabels <- function(gene, table) {
#   gene.table <- filter(table, geneSymbol == gene)
#   gene.table <- makeRowLabel(gene.table)
#   row_labels <- unique(gene.table$row_label)
#   return(row_labels)
# }

##### function to generate scatter plot
#  myScatterPlot <- function(table, params) {
# 
#   # extract parameters
#   gene.x <- params$gene.x
#   gene.y <- params$gene.y
#   x.label <- params$x
#   y.label <- params$y
#   
#   # extract rows for selected genes and row labels
#   table <- filter(table, geneSymbol %in% c(gene.x, gene.y))
#   table <- makeRowLabel(table)
#   table.x <- filter(table, geneSymbol == gene.x, row_label == x.label)
#   table.y <- filter(table, geneSymbol == gene.y, row_label == y.label)
#   
#   # combine both tables together (renaming columns)
#   scatter.table <- full_join(table.x, table.y, by=c("participant", "expt"))
#   
#   if (length(unique(scatter.table$expt)) > 1) {
#     point <- geom_point(data = scatter.table,
#                         aes(x=ratio.x, y=ratio.y, col=expt), na.rm=T)
#   } else {
#     point <- geom_point(data = scatter.table,
#                         aes(x=ratio.x, y=ratio.y), na.rm=T)
#   }
#   # make plot
#   scatter.plot <- ggplot(data = scatter.table) +
#     point +
#     labs(x = paste(scatter.table$geneSymbol.x, scatter.table$row_label.x, sep=': '),
#          y = paste(scatter.table$geneSymbol.y, scatter.table$row_label.y, sep=': '),
#          title = "Plot two -omes against each other")
#   
#   return(scatter.plot)
# 
# }

## function to make HLA tables
makeHLATables <- function (gene, hla.table) {
  participants.unique <- unique(hla.table$directory)

  # filter for sequences matching gene
  table <- filter(hla.table, geneSymbol == gene)
  
  # create appropriate row labels
  table$row_label <- ifelse(table$PTM_nuORF == "", table$sequence, 
                            paste(table$sequence, ' (', table$PTM_nuORF, ')', sep=''))

  # class 1
  class = 'cls1'
  class1.table <- filter(table, hla_class == class)
  row_label.unique <- unique(class1.table$row_label)
  final.class1.table <- data.frame(matrix(nrow = length(row_label.unique),
                                         ncol = length(participants.unique)))
  rownames(final.class1.table) <- sort(row_label.unique)
  names(final.class1.table) <- sort(participants.unique)
  for (p in participants.unique) {
    p.row_labels <- class1.table$row_label[which(class1.table$directory == p)]
    stopifnot(length(unique(p.row_labels)) == length(p.row_labels))
    final.class1.table[p.row_labels, p] <- "\U2713"
  }
  final.class1.table[is.na(final.class1.table)] <- ' '


  # class 2
  class = 'cls2'
  class2.table <- filter(table, hla_class == class)
  row_labels.unique <- unique(class2.table$row_label)
  final.class2.table <- data.frame(matrix(nrow = length(row_labels.unique),
                                         ncol = length(participants.unique)))
  rownames(final.class2.table) <- sort(row_labels.unique)
  names(final.class2.table) <- sort(participants.unique)
  for (p in participants.unique) {
    p.row_labels <- class2.table$row_label[which(class2.table$directory == p)]
    stopifnot(length(unique(p.row_labels)) == length(p.row_labels))
    final.class2.table[p.row_labels, p] <- "\U2713"
  }
  final.class2.table[is.na(final.class2.table)] <- ' '

  # add column for sequence
  final.class1.table$Sequence <- rownames(final.class1.table)
  final.class2.table$Sequence <- rownames(final.class2.table)
  
  # reorder columns
  final.class1.table <- final.class1.table[ , c(11, 1:10)]
  final.class2.table <- final.class2.table[ , c(11, 1:10)]


  return (list(cls1 = final.class1.table, cls2 = final.class2.table))
}

makeComplexLegend <- function(column.anno.col) {
  legend_list <- list()
  for (i in c(length(column.anno.col), 1:(length(column.anno.col)-1))) {
    title = names(column.anno.col)[i]
    color = column.anno.col[[title]]
    
    if (typeof(color) == 'character') {
      lgd <- Legend(labels = setdiff(names(color), 'NA'), title = title, legend_gp=gpar(fill=color))
      
    } else {
      lgd <- Legend(col_fun = color, title = title, direction = 'horizontal')
    }
    legend_list <- c(legend_list, lgd)
    
    
  }
  
  pd = packLegend(list = legend_list, 
                  column_gap = unit(8, "mm"),
                  direction = 'horizontal', 
                  max_width = unit(20, 'cm'))
  return(pd)
}

###############################################################################
## FOR TESTING PURPOSES ONLY

# parameters
# params <- list(
#   genes.char = c('EGFR', 'RB1', 'KRAS', 'STK11', "TP53"),
#   zscore = 'row',
#   PTMsites = 'most variable',
#   min.val = -3,
#   max.val = 3,
#   sort.after = 'Multi.omic.subtype',
#   id = "MonteTab")
# 
# table <- monte_table
# 
# HM.out <- myComplexHeatmap(table, params)
# 
# ## test download PDF
# file <- paste(FILENAMESTRING,  '-',  gsub(' |\\:','-', Sys.time()), '.pdf', sep='')
# pdf(
#   file = file,
#   width = 1400/72,
#   height = dynamicHeightHM(nrow(HM.out$Table))/72 # TODO: add parenthesis
# )
# HM.out$HM # TODO: add parentheses
# dev.off()

## test scatter plot
# params <- list(gene.x = "EGFR",
#                gene.y = "BAK1",
#                x = "prot ENSP00000275493.2",
#                y = "prot ENSP00000363591.3")
