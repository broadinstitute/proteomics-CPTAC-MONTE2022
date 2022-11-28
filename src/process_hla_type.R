hla.type.table <- read.csv("original-data/Table_S4_11-28-2022.csv")

hla.type.table$Participant <- gsub('-','.',hla.type.table$Participant)
hla.type.table <- hla.type.table[, c(2, 14:19)]

load("data/MONTE_data_2022-08-11.RData")

save(list = ls(), file = 'data/MONTE_data_2022-11-28.RData')