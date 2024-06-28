library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(harmony)
library(tidyr)
library(ggpubr)
library(tidyverse)
library(scico)

setwd("~/Desktop/Howe")

# howe <- readRDS(file = "objects/seurat_object_final.rds")

howe <- readRDS(file = "objects/seurat_object.rds")
DefaultAssay(howe) <- 'RNA'
howe <- NormalizeData(object = howe)
howe <- FindVariableFeatures(object=howe)
howe <- ScaleData(object = howe)

new.cluster.ids <- c(
  "DN_T_1","EC_1","VSMC_1","CD8T_1","Macrophage_1","Macrophage_2",
  "VSMC_FB","VSMC_2","CD8T_2","Dendritic","B_Cell","EC_2","EC_3",
  "NKT","Unknown_immune","DN_T_2","VSMC_3","Plasma_cell","Mast_cell",
  "CD8_TRM")

names(new.cluster.ids) <- levels(howe)
howe <- RenameIdents(howe, new.cluster.ids)
DimPlot(howe, label = TRUE, pt.size = 0.5)

howe@meta.data$cell_type = howe@active.ident


my_col <- setNames(c("#D24D29", "#CC2841", "#8A1837", "#FAF6C7", "#FDD167",
                     "#F0EA37", "#D2D1E6", "#8987BD", "#563B93", "#322767",
                     "#D5E9F5", "#BAD3EC", "#8FBEDB", "#589ECE", "#2D80BC",
                     "#055EA8", '#D9EBD2', "#C5DFBA", "#9BD0A3", "#6AC3B8"),
                   c("EC_1", "EC_2", "EC_3", "Macrophage_1", "Macrophage_2", 
                     "Dendritic", "VSMC_1", "VSMC_2", "VSMC_3", "VSMC_FB", 
                     "CD8T_1", "CD8T_2", "DN_T_1", "DN_T_2", "NKT", "CD8_TRM", 
                     "B_Cell", "Plasma_cell", "Mast_cell", "Unknown_immune"))
            
            
# DimPlot Figures
svg("images/UMAP/umap_labelled.svg", width=12, height=12)
DimPlot(howe, pt.size = 0.5)
dev.off()
svg("images/UMAP/umap_labelled_col.svg", width=12, height=12)
DimPlot(howe, cols=my_col, pt.size = 0.5)
dev.off()
svg("images/UMAP/umap_condition.svg", width=12, height=12)
DimPlot(howe, group.by='condition', pt.size = 0.5)
dev.off()
svg("images/UMAP/umap_sample.svg", width=12, height=12)
DimPlot(howe, group.by='sample', pt.size = 0.5)
dev.off()

# Cell Annotation Figures
svg("images/UMAP/ft_CD3E_TCell.svg", width=12, height=12)
FeaturePlot(howe, features = 'CD3E')
dev.off()
svg("images/UMAP/ft_AIF1_Monocyte_Dendritic.svg", width=12, height=12)
FeaturePlot(howe, features = 'AIF1')
dev.off()
svg("images/UMAP/ft_MYL9_VSMC.svg", width=12, height=12)
FeaturePlot(howe, features = 'MYL9')
dev.off()
svg("images/UMAP/ft_VWF_EC.svg", width=12, height=12)
FeaturePlot(howe, features = 'VWF')
dev.off()
svg("images/UMAP/ft_MS4A1_BCell.svg", width=12, height=12)
FeaturePlot(howe, features = 'MS4A1')
dev.off()
svg("images/UMAP/ft_MZB1_PlasmaCell.svg", width=12, height=12)
FeaturePlot(howe, features = 'MZB1')
dev.off()
svg("images/UMAP/ft_TPSAB1_Mast.svg", width=12, height=12)
FeaturePlot(howe, features = 'TPSAB1')
dev.off()
svg("images/UMAP/dotplot_TCell_annot.svg", width=12, height=6)
DotPlot(object = howe, features = c('CD2','CD3E','CD8A','CD4','NKG7','TRDC','ITGAE'), 
        idents=c('DN_T_1','DN_T_2','CD8T_1','CD8T_2','CD8_TRM','NKT','Unknown_immune'),
        dot.min=.22)
dev.off()
svg("images/UMAP/dotplot_VSMC_annot.svg", width=12, height=6)
DotPlot(object = howe, features = c('MYL9',"MYH11",'LUM'), 
        idents=c('VSMC_1','VSMC_2','VSMC_3','VSMC_FB'),
        dot.min=.3)
dev.off()
svg("images/UMAP/dotplot_MonoDC_annot.svg", width=12, height=6)
DotPlot(object = howe, features = c('CD14',"APOBEC3A",'ITGAX','SIRPA'), 
        idents=c('Macrophage_1','Macrophage_2','Dendritic'),
        dot.min=.09)
dev.off()


# Loading in gene-lists for module scoring
df_all <- read.csv('data/all_pm.csv')

df_pm <- subset(df_all, df_all$q.value<0.05 & df_all$log2FC > 0)
df_pm_stringent <- subset(df_all, df_all$q.value<0.05 & df_all$log2FC > 1)
df_mp <- subset(df_all, df_all$q.value<0.05 & df_all$log2FC < 0)
df_mp_stringent <- subset(df_all, df_all$q.value<0.05 & df_all$log2FC < -1)

# Extract gene names into a list
all_genes <- c(df_all$T..GeneID)
pm_genes <- c(df_pm$T..GeneID)
pm_s_genes <- c(df_pm_stringent$T..GeneID)
mp_genes <- c(df_mp$T..GeneID)
mp_s_genes <- c(df_mp_stringent$T..GeneID)

howe <- AddModuleScore(object = howe,
                       features = list(all_genes),
                       name = 'all_pm', search=TRUE,assay='RNA')

howe <- AddModuleScore(object = howe,
                       features = list(pm_genes),
                       name = 'DE_pm', search=TRUE,assay='RNA')

howe <- AddModuleScore(object = howe,
                       features = list(pm_s_genes),
                       name = 'DE_pm_stringent', search=TRUE,assay='RNA')

howe <- AddModuleScore(object = howe,
                       features = list(mp_genes),
                       name = 'DE_mp', search=TRUE,assay='RNA')

howe <- AddModuleScore(object = howe,
                       features = list(mp_s_genes),
                       name = 'DE_mp_stringent', search=TRUE,assay='RNA')

# ModuleScore Figures (BoxPlot)

svg("images/module_score_box/all_genes.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$all_pm1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "all_pm1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

svg("images/module_score_box/DE_pm.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$DE_pm1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "DE_pm1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/DE_pm_stringent.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$DE_pm_stringent1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "DE_pm_stringent1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/DE_pm_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$DE_pm1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "DE_pm1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()
svg("images/module_score_box/DE_pm_stringent_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$DE_pm_stringent1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "DE_pm_stringent1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()

svg("images/module_score_box/DE_mp.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$DE_mp1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "DE_mp1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/DE_mp_stringent.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$DE_mp_stringent1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "DE_mp_stringent1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/DE_mp_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$DE_mp1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "DE_mp1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()
svg("images/module_score_box/DE_mp_stringent_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$DE_mp_stringent1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "DE_mp_stringent1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()


# Asymptomatic Module Scoring
df_as <- read.csv('data/as_pm.csv')
df_as_pm <- subset(df_as, df_as$q.value<0.05 & df_as$log2FC > 0)
df_as_pm_stringent <- subset(df_as, df_as$q.value<0.05 & df_as$log2FC > 1)
df_as_mp <- subset(df_as, df_as$q.value<0.05 & df_as$log2FC < 0)
df_as_mp_stringent <- subset(df_as, df_as$q.value<0.05 & df_as$log2FC < -1)

a_pm_genes <- c(df_as_pm$T..GeneID)
a_pm_s_genes <- c(df_as_pm_stringent$T..GeneID)
a_mp_genes <- c(df_as_mp$T..GeneID)
a_mp_s_genes <- c(df_as_mp_stringent$T..GeneID)

howe <- AddModuleScore(object = howe,
                       features = list(a_pm_genes),
                       name = 'asymp_DE_pm', search=TRUE,assay='RNA')

howe <- AddModuleScore(object = howe,
                       features = list(a_pm_s_genes),
                       name = 'asymp_DE_pm_stringent', search=TRUE,assay='RNA')

howe <- AddModuleScore(object = howe,
                       features = list(a_mp_genes),
                       name = 'asymp_DE_mp', search=TRUE,assay='RNA')

howe <- AddModuleScore(object = howe,
                       features = list(a_mp_s_genes),
                       name = 'asymp_DE_mp_stringent', search=TRUE,assay='RNA')

svg("images/module_score_box/asymp_DE_pm.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$asymp_DE_pm1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "asymp_DE_pm1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/asymp_DE_pm_stringent.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$asymp_DE_pm_stringent1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "asymp_DE_pm_stringent1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/asymp_DE_pm_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$asymp_DE_pm1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "asymp_DE_pm1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()
svg("images/module_score_box/asymp_DE_pm_stringent_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$asymp_DE_pm_stringent1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "asymp_DE_pm_stringent1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()


svg("images/module_score_box/asymp_DE_mp.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$asymp_DE_mp1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "asymp_DE_mp1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/asymp_DE_mp_stringent.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$asymp_DE_mp_stringent1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "asymp_DE_mp_stringent1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/asymp_DE_mp_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$asymp_DE_mp1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "asymp_DE_mp1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()
svg("images/module_score_box/asymp_DE_mp_stringent_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$asymp_DE_mp_stringent1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "asymp_DE_mp_stringent1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()



# Symptomatic module scoring

df_s <- read.csv('data/s_pm.csv')
df_s_pm <- subset(df_s, df_s$q.value<0.05 & df_s$log2FC > 0)
df_s_pm_stringent <- subset(df_s, df_s$q.value<0.05 & df_s$log2FC > 1)
df_s_mp <- subset(df_s, df_s$q.value<0.05 & df_s$log2FC < 0)
df_s_mp_stringent <- subset(df_s, df_s$q.value<0.05 & df_s$log2FC < -1)

s_pm_genes <- c(df_s_pm$T..GeneID)
s_pm_s_genes <- c(df_s_pm_stringent$T..GeneID)
s_mp_genes <- c(df_s_mp$T..GeneID)
s_mp_s_genes <- c(df_s_mp_stringent$T..GeneID)

howe <- AddModuleScore(object = howe,
                       features = list(s_pm_genes),
                       name = 'symp_DE_pm', search=TRUE,assay='RNA')
howe <- AddModuleScore(object = howe,
                       features = list(s_pm_s_genes),
                       name = 'symp_DE_pm_stringent', search=TRUE,assay='RNA')
howe <- AddModuleScore(object = howe,
                       features = list(s_mp_genes),
                       name = 'symp_DE_mp', search=TRUE,assay='RNA')
howe <- AddModuleScore(object = howe,
                       features = list(s_mp_s_genes),
                       name = 'symp_DE_mp_stringent', search=TRUE,assay='RNA')

svg("images/module_score_box/symp_DE_pm.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$symp_DE_pm1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "symp_DE_pm1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/symp_DE_pm_stringent.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$symp_DE_pm_stringent1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "symp_DE_pm_stringent1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/symp_DE_pm_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$symp_DE_pm1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "symp_DE_pm1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()
svg("images/module_score_box/symp_DE_pm_stringent_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$symp_DE_pm_stringent1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "symp_DE_pm_stringent1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()


svg("images/module_score_box/symp_DE_mp.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$symp_DE_mp1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "symp_DE_mp1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/symp_DE_mp_stringent.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$symp_DE_mp_stringent1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "symp_DE_mp_stringent1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/symp_DE_mp_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$symp_DE_mp1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "symp_DE_mp1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()
svg("images/module_score_box/symp_DE_mp_stringent_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$symp_DE_mp_stringent1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "symp_DE_mp_stringent1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()



# miRNA analysis
df_m <- read.csv('data/miRNA.csv',header=TRUE)
df_m <- subset(df_m, df_m$FDR<0.01)
df_m_genes <- c(df_m$Gene.Symbol)

df_m <- read.csv('data/miRNA.csv',header=TRUE)
df_m_s <- subset(df_m, df_m$FDR<0.001)
df_m_s_genes <- c(df_m_s$Gene.Symbol)

howe <- AddModuleScore(object = howe,
                       features = list(df_m_genes),
                       name = 'miRNA_pm')
howe <- AddModuleScore(object = howe,
                       features = list(df_m_s_genes),
                       name = 'miRNA_pm_stringent')

svg("images/module_score_box/mirna/miRNA_pm.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$miRNA_pm1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "miRNA_pm1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/mirna/miRNA_pm_stringent.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$miRNA_pm_stringent1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "miRNA_pm_stringent1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/mirna/miRNA_pm_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$miRNA_pm1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "miRNA_pm1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()
svg("images/module_score_box/mirna/miRNA_pm_stringent_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$miRNA_pm_stringent1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "miRNA_pm_stringent1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()


# Asymptomatic miRNA
df_asymp_m <- read.csv('data/miRNA_asymp.csv',header=TRUE)
df_asymp_m <- subset(df_asymp_m, df_asymp_m$FDR<0.01)
df_asymp_m_genes <- c(df_asymp_m$Gene.Symbol)

df_asymp_m_s <- read.csv('data/miRNA_asymp.csv',header=TRUE)
df_asymp_m_s <- subset(df_asymp_m_s, df_asymp_m_s$FDR<0.001)
df_asymp_m_s_genes <- c(df_asymp_m_s$Gene.Symbol)

howe <- AddModuleScore(object = howe,
                       features = list(df_asymp_m_genes),
                       name = 'miRNA_asymp_pm')
howe <- AddModuleScore(object = howe,
                       features = list(df_asymp_m_s_genes),
                       name = 'miRNA_asymp_pm_stringent')

svg("images/module_score_box/mirna/asymp_miRNA_pm.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$miRNA_asymp_pm1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "miRNA_asymp_pm1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/mirna/asymp_miRNA_pm_stringent.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$miRNA_asymp_pm_stringent1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "miRNA_asymp_pm_stringent1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/mirna/asymp_miRNA_pm_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$miRNA_asymp_pm1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "miRNA_asymp_pm1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()
svg("images/module_score_box/mirna/asymp_miRNA_pm_stringent_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$miRNA_asymp_pm_stringent1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "miRNA_asymp_pm_stringent1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()


# Symptomatic miRNA
df_symp_m <- read.csv('data/miRNA_symp.csv',header=TRUE)
df_symp_m <- subset(df_symp_m, df_symp_m$FDR<0.01)
df_symp_m_genes <- c(df_symp_m$Gene.Symbol)

df_symp_m_s <- read.csv('data/miRNA_symp.csv',header=TRUE)
df_symp_m_s <- subset(df_symp_m_s, df_symp_m_s$FDR<0.001)
df_symp_m_s_genes <- c(df_symp_m_s$Gene.Symbol)

howe <- AddModuleScore(object = howe,
                       features = list(df_symp_m_genes),
                       name = 'miRNA_symp_pm')
howe <- AddModuleScore(object = howe,
                       features = list(df_symp_m_s_genes),
                       name = 'miRNA_symp_pm_stringent')

svg("images/module_score_box/mirna/symp_miRNA_pm.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$miRNA_symp_pm1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "miRNA_symp_pm1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/mirna/symp_miRNA_pm_stringent.svg",  width=12, height=9)
order = reorder(howe$cell_type, -howe@meta.data$miRNA_symp_pm_stringent1)
howe$celltype <- factor(howe$cell_type, levels = levels(order))
ggboxplot(howe@meta.data, x = "celltype", y = "miRNA_symp_pm_stringent1", fill='celltype',
          palette = my_col,legend = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
svg("images/module_score_box/mirna/symp_miRNA_pm_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$miRNA_symp_pm1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "miRNA_symp_pm1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()
svg("images/module_score_box/mirna/symp_miRNA_pm_stringent_condition.svg",  width=12, height=9)
order = reorder(howe$condition, -howe@meta.data$miRNA_symp_pm_stringent1)
howe$condition <- factor(howe$condition, levels = levels(order))
ggboxplot(howe@meta.data, x = "condition", y = "miRNA_symp_pm_stringent1", fill='condition',
          legend = "none") +
  stat_compare_means(comparisons = list(c('AC','PA')), method = "wilcox.test")
dev.off()






##### DELETE?

# Statistical tests
# https://gist.github.com/makloufg/8e0456fddc56ebef2d9b67576b131482

order = reorder(howe$cell_type, -howe@meta.data$all_pm1)
howe$celltype_M1 <- factor(howe$cell_type, levels = levels(order))

# Perform the test
compare_means(all_pm1 ~ celltype_M1, data = howe@meta.data,
              ref.group = ".all.", method = "wilcox")
compare_means(all_pm1 ~ celltype_M1, data = howe@meta.data,
              method = "kruskal.test")

p1 <- ggboxplot(howe@meta.data, x = "celltype_M1", y = "all_pm1", fill='celltype_M1',palette = my_col,
                legend = "none")
p2 <- ggviolin(howe@meta.data, x = "celltype_M1", y = "all_pm1", fill = "celltype_M1",palette = my_col,
               legend='none')

sp1 = p1 +
  stat_compare_means(label = "p.signif", method = "wilcox", ref.group = ".all.", hide.ns = FALSE, show.legend = FALSE) +
  stat_compare_means(method = "kruskal.test", label.x = 2, label.y = 0.6) + # Add global the p-value 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = mean(howe@meta.data$all_pm1), linetype = 2) # Add horizontal line at base mean

sp1
p1
p2

sp2 = p2 +
  stat_compare_means(label = "p.signif", method = "wilcox", ref.group = ".all.", hide.ns = FALSE, show.legend = FALSE) +
  stat_compare_means(method = "kruskal.test", label.x = 2, label.y = 0.6) + # Add global the p-value 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = mean(howe@meta.data$all_pm1), linetype = 2) # Add horizontal line at base mean
sp2


# https://stackoverflow.com/questions/68666942/how-to-change-p-value-label-in-ggpubr-ggplot2 
# https://www.biostars.org/p/458261/

compare_means(miRNA_ft1 ~ cell_type,
              data = howe@meta.data, 
              method = "kruskal.test", paired = FALSE)

# Adding statistical test/p-value

# svg("signaling_pathways.svg",  width=20, height=16)
VlnPlot(howe, group.by='cell_type',features="all_H_features1",sort=TRUE, y.max = 2) + 
  stat_compare_means(comparisons = list(c('DN_T_2','CD8_TRM')), method = "wilcox.test")





