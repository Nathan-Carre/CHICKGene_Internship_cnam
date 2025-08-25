df <- read.csv("resBatchEffectSex_Tecnica_CHIKR_oui_CC_DMS_ChikRxCC_non-DMS_ChikR.csv", header = TRUE)
colnames(df)[1] <- "gene"

df_sorted <- df[order(df$pvalue, decreasing = FALSE), ]

write.csv(df_sorted, "SORTED_p-value_resBatchEffectSex_Tecnica_CHIKR_oui_CC_DMS_ChikRxCC_non-DMS_ChikR.csv", row.names = FALSE)

head(df_sorted)

