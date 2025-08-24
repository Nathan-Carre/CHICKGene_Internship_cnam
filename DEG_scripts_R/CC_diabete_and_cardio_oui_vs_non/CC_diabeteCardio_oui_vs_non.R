#WITH BATCH EFFECT

countdata <- read.table("filtered_cut_219_CountsChikgene.txt",header=TRUE,row.names=1,check.names = FALSE)
countdata <- as.matrix(countdata)

coldata <- read.csv("pheno_CC_diabeteCardio_oui_vs_non.csv", header=TRUE)
coldata$condition <- factor(coldata$condition, levels = c("0","1"), labels = c("NC", "C"))
coldata$diabeteCardio_oui_vs_non <- factor(coldata$diabeteCardio_oui_vs_non, levels = c("0","1", "2"), labels = c("CC_diabeteCardio_non", "CC_diabeteCardio_oui", "Others"))
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
                              design = ~ sex + tecnica + diabeteCardio_oui_vs_non)

dds
dds <- dds[ rowSums(counts(dds)) > 100, ]
dds
dds <- DESeq(dds)
resultsNames(dds)

ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

res <- results(dds, contrast = c("diabeteCardio_oui_vs_non", "CC_diabeteCardio_oui", "CC_diabeteCardio_non"))
res
write.csv(as.data.frame(res),file="resBatchEffect_SexTecnicaDiabeteCardioCC_diabeteCardio_nonVSCC_diabeteCardio_non.csv")

write.csv(as.data.frame(assay(ntd)),file="ntd_BatchEffect_expression_SexTecnicaDiabeteCardioCC_diabeteCardio_nonVSCC_diabeteCardio_non.csv")
write.csv(as.data.frame(assay(vsd)),file="vsd_BatchEffect_expression_SexTecnicaDiabeteCardioCC_diabeteCardio_nonVSCC_diabeteCardio_non.csv")

norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)
write.csv(as.data.frame(norm.counts),file="norm_counts_batch_ConditionDiabeteCardioCC_diabeteCardio_nonVSCC_diabeteCardio_non.csv")
write.csv(as.data.frame(log.norm.counts),file="log_norm_counts_batch_ConditionDiabeteCardioCC_diabeteCardio_nonVSCC_diabeteCardio_non.csv")
counts <- counts(dds, normalized=FALSE)
write.csv(as.data.frame(norm.counts),file="counts_batch_ConditionDiabeteCardioCC_diabeteCardio_nonVSCC_diabeteCardio_non.csv")





#Volcano plot
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
devtools::install_github('kevinblighe/EnhancedVolcano')

library(EnhancedVolcano)

jpeg("volcano_Batch_CFSxC.jpg", width=1000, height=1000)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
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

jpeg("volcano_Batch_CFSxC_padj.jpg", width=1000, height=1000)
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