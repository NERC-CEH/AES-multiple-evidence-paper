#Correlations

buttabund <- read.csv("LandSpAES butterfly abundance.csv")

UKBMS_all_data <- read.csv("UKBMS butterfly data.csv")

WCBS_all_data <- read.csv("WCBS butterfly data.csv")

LS_corrs <- cor(buttabund[,c(4,5,10,13,16)], method = "spearman")

UKBMS_corrs <- cor(UKBMS_all_data[,c(18,19,9,12,15)], method = "spearman")

WCBS_corrs <- cor(WCBS_all_data[,c(18,19,9,12,15)], method = "spearman")

##extract AES and PCA gradient correlations

corr_table <- rbind(LS_corrs[3:5,1:2], UKBMS_corrs[3:5,1:2])

corr_table <- data.frame(rbind(corr_table, WCBS_corrs[3:5,1:2]))

corr_table$Dataset <- c(rep("LandSpAES",3),rep("UKBMS",3), rep("WCBS", 3))

corr_table$PCA_variable <- rep(c("Climate PC1", "Landscape PC1", "Habitat PC1"),3)

corr_table <- corr_table[,c(3,4,1,2)]
row.names(corr_table) <- 1:9
corr_table[,3:4] <- round(corr_table[3:4],3)

write.csv(corr_table, "Table of correlations for SI.csv")
