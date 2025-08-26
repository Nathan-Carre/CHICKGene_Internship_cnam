#WITH BATCH EFFECT

countdata <- read.table("filtered_cut_219_CountsChikgene.txt",header=TRUE,row.names=1,check.names = FALSE)
countdata <- as.matrix(countdata)

coldata <- read.csv("pheno_CC_cardio_inc_tabac.csv", header=TRUE)
coldata$condition <- factor(coldata$condition, levels = c("0","1"), labels = c("NC", "C"))
coldata$cardio_inc_tabac <- factor(coldata$cardio_inc_tabac, levels = c("0","1", "2"), labels = c("CC_cardio_inc_no_tabac", "CC_cardio_inc_tabac", "Others"))
coldata$year <- factor(coldata$year)
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
coldata$chikr_oui <- factor(coldata$chikr_oui, levels = c("0","1"), labels = c("CC_non-chikr", "CC_chikr"))
coldata$age_45 <- factor(coldata$age_45)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ sex + tecnica + cardio_inc_tabac)

dds
dds <- dds[ rowSums(counts(dds)) > 100, ]
dds
dds <- DESeq(dds)
resultsNames(dds)

ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

res <- results(dds, contrast = c("cardio_inc_tabac", "CC_cardio_inc_tabac", "CC_cardio_inc_no_tabac"))
res
res_df <- as.data.frame(res)
save(res_df, file = "res_Cardio_inc_tabac.RData")


write.csv(as.data.frame(res),file="resBatchEffect_SexTecnicaCardioInc_CC_cardio_inc_tabac_vs_CC_cardio_inc_no_tabac.csv")

write.csv(as.data.frame(assay(ntd)),file="ntd_BatchEffect_expression_SexTecnicaCardioInc_CC_cardio_inc_tabac_vs_CC_cardio_inc_no_tabac.csv")
write.csv(as.data.frame(assay(vsd)),file="vsd_BatchEffect_expression_SexTecnicaCardioInc_CC_cardio_inc_tabac_vs_CC_cardio_inc_no_tabac.csv")

norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)
write.csv(as.data.frame(norm.counts),file="norm_counts_batch_CC_cardio_inc_tabac_vs_CC_cardio_inc_no_tabac.csv")
write.csv(as.data.frame(log.norm.counts),file="log_norm_counts_batch_CC_cardio_inc_tabac_vs_CC_cardio_inc_no_tabac.csv")
counts <- counts(dds, normalized=FALSE)
write.csv(as.data.frame(norm.counts),file="counts_batch_CC_cardio_inc_tabac_vs_CC_cardio_inc_no_tabac.csv")


load("res_Cardio_inc_tabac.RData")

res <- res_df  # si tu veux garder le nom "res"


if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}

BiocManager::install('EnhancedVolcano', force = TRUE)


library(EnhancedVolcano)


jpeg("volcano_BatchEffect_SexTecnicaCardioInc_CC_cardio_inc_tabac_vs_CC_cardio_inc_no_tabac_p-value_0.7_res300.jpg", width = 10, height = 10, units = "in", res = 300)
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

jpeg("volcano_BatchEffect_SexTecnicaCardioInc_CC_cardio_inc_tabac_vs_CC_cardio_inc_no_tabac_p_adj_1_res300.jpg", width=1000, height=1000)
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