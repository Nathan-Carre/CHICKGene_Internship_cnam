#WITH BATCH EFFECT

countdata <- read.table("CountsChikgeneFinal133.txt",header=TRUE,row.names=1,check.names = FALSE)
countdata <- as.matrix(countdata)

coldata <- read.table("pheno_133_DMS_CHIKR.csv", header=TRUE, sep=",", check.names=FALSE, stringsAsFactors=FALSE)
coldata$individuo <- factor(coldata$individuo)
# coldata$condition <- factor(coldata$condition, levels = c("0","1"), labels = c("NC", "C"))
coldata$condition_dms <- factor(coldata$condition_dms, levels = c("0","1"), labels = c("CC_non-DMS", "CC_DMS"))
# coldata$CC_DMS <- factor(coldata$CC_DMS, levels = c("0","1", "2"), labels = c("CC_non-DMS", "CC_DMS", "NCC"))
coldata$cohorte <- factor(coldata$cohorte)
coldata$age <- factor(coldata$age)
coldata$year <- factor(coldata$year)
# coldata$year2 <- factor(coldata$year2)
coldata$sex <- factor(coldata$sex)
coldata$tecnica <- factor(coldata$tecnica)
coldata$maternite <- factor(coldata$maternite)
coldata$diabete <- factor(coldata$diabete)
coldata$diabete_pre <- factor(coldata$diabete_pre)
coldata$diabete_inc <- factor(coldata$diabete_inc)
coldata$cardiopathie <- factor(coldata$cardiopathie)
coldata$cardio_pre <- factor(coldata$cardio_pre)
coldata$cardio_inc <- factor(coldata$cardio_inc)
coldata$tabac <- factor(coldata$tabac)
coldata$imc <- factor(coldata$imc)
coldata$chikr_oui <- factor(coldata$chikr_oui, levels = c("0","1", "2"), labels = c("CC_non-chikr", "CC_chikr", "Others"))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ sex + tecnica + chikr_oui + condition_dms)

dds
dds <- dds[ rowSums(counts(dds)) > 100, ]
dds
dds <- DESeq(dds)
resultsNames(dds)

ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Bien faire attention à la condition sélectionné : condition_DMS_NCC_DMS_vs_NCC_non.DMS

ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

res <- results(dds, contrast = c("condition_dms", "CC_DMS", "CC_non-DMS"))
res
res_df <- as.data.frame(res)
save(res_df, file = "res_for_volcano_DMS.RData")


write.csv(as.data.frame(res),file="resBatchEffectSex_Tecnica_CHIKR_oui_CC_DMS_ChikRxCC_non-DMS_ChikR.csv")

norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)
write.csv(as.data.frame(norm.counts),file="norm_counts_batch_CC_DMS_non-chikrxCC_non-DMS_non-chikr.csv")
write.csv(as.data.frame(log.norm.counts),file="log_norm_counts_batch_CC_DMS_non-chikrxCC_non-DMS_non-chikr.csv")
counts <- counts(dds, normalized=FALSE)
write.csv(as.data.frame(norm.counts),file="counts_batch_CC_DMS_non-chikrxCC_non-DMS_non-chikr..csv")

# Pas de batch dans le design de Deseq2 et gérer les Batch que dans le Limma


load("res_for_volcano_DMS.RData")

res <- res_df  # si tu veux garder le nom "res"


if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}

BiocManager::install('EnhancedVolcano', force = TRUE)


library(EnhancedVolcano)


jpeg("volcano_Batch_Sex_Tecnica_CHIKR_oui_CC_DMSxCC_non-DMS_p-value_0.7_res300.jpg", width = 10, height = 10, units = "in", res = 300)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-5,
                FCcutoff = 0.7,
                pointSize = 1.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')
dev.off()

jpeg("volcano_Batch_Sex_Tecnica_CHIKR_oui_CC_DMSxCC_non-DMS_p_adj_1_res300.jpg", width=1000, height=1000)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-5,
                FCcutoff = 1.0,
                pointSize = 1.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')
dev.off()