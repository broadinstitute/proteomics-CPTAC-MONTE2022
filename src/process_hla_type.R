# hla.type.table <- read.csv("original-data/Table_S4_11-28-2022.csv")
# 
# hla.type.table$Participant <- gsub('-','.',hla.type.table$Participant)
# hla.type.table <- hla.type.table[, c(2, 14:19)]
# 
# load("data/MONTE_data_2022-08-11.RData")
# 
# save(list = ls(), file = 'data/MONTE_data_2022-11-28.RData')

hla.athena <- read.csv('original-data/HLA_1_data_example_11-28-2022.csv')
hla.athena$directory <- gsub('-', '.', hla.athena$directory)

# COMBINATIONS OF geneSymbol + sequence DO NOT ALL HAVE THE SAME TYPE
# out <- df %>%
#   select(geneSymbol, sequence, allele, Msi_HLAthena, Binder) %>%
#   group_by(geneSymbol, sequence) %>% 
#   summarise(n_distinct(allele, Msi_HLAthena, Binder), .groups = 'keep')

hla.athena <- select(hla.athena, directory, geneSymbol, sequence, allele, Msi_HLAthena, Binder)

load("data/MONTE_data_2022-11-28.RData")

# merge hla tables
hla.athena$hla_class <- rep('cls1', nrow(hla.athena))

o <- merge(hla.table, hla.athena, by = c('directory', 'geneSymbol', 'sequence'), all.x = T)



save(list = ls(), file = 'data/MONTE_data_2022-11-28.RData')