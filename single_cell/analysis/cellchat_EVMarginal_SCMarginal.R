library(CellChat)
library(Seurat)
library(reshape2)
library(patchwork)
library(stringr)
library(RColorBrewer)
library(ComplexHeatmap)

# Load Seurat Object
howe <- readRDS(file = "objects/seurat_object.rds")
new.cluster.ids <- c(
  "DN_T_1","EC_1","VSMC_1","CD8T_1","Macrophage_1","Macrophage_2",
  "VSMC_FB","VSMC_2","CD8T_2","Dendritic","B_Cell","EC_2","EC_3",
  "NKT","Unknown_immune","DN_T_2","VSMC_3","Plasma_cell","Mast_cell",
  "CD8_TRM")
names(new.cluster.ids) <- levels(howe)
howe <- RenameIdents(howe, new.cluster.ids)
DimPlot(howe, label = TRUE, pt.size = 0.5)
howe@meta.data$cell_type = howe@active.ident

# Split Seurat Object

howe_ac <- subset(x = howe, subset = condition == "AC")
howe_pa <- subset(x = howe, subset = condition == "PA")

# Read Gene names
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


# CellChat Analysis
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

# Subsetting cellchat database
interaction_db <- CellChatDB$interaction

contains_any_gene <- function(string, genes) {
  any(sapply(genes, function(gene) str_detect(string, paste0("\\b", gene, "\\b"))))
}

# Apply this function to each row of the dataframe and create a logical vector
matches <- apply(interaction_db, 1, function(row) {
  contains_any_gene(row["ligand"], mp_genes)
})

# Subset the dataframe
subset_interactions <- interaction_db[matches, ]
CellChatDB.use$interaction <- subset_interactions


cc_full <- createCellChat(object = howe_pa, group.by = "ident", assay = "SCT")
cc_full@DB <- CellChatDB.use

cellchat <- subsetData(cc_full)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


my_col <- setNames(c("#8FBEDB", "#D24D29", "#D2D1E6", "#D5E9F5", "#FAF6C7", 
                     "#FDD167", "#322767", "#8987BD", "#BAD3EC", "#F0EA37",
                     "#D9EBD2", "#CC2841", "#8A1837", "#2D80BC", "#6AC3B8",
                     "#589ECE", "#563B93", "#C5DFBA", "#9BD0A3", "#055EA8"),
                   c("DN_T_1", "EC_1", "VSMC_1", "CD8T_1", "Macrophage_1",
                     "Macrophage_2", "VSMC_FB", "VSMC_2", "CD8T_2", "Dendritic",
                     "B_Cell", "EC_2", "EC_3", "NKT", "Unknown_immune",
                     "DN_T_2", "VSMC_3", "Plasma_cell", "Mast_cell", "CD8_TRM"))

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                         width=36, height=30, font.size = 14,
                                         color.use=my_col)
ht2


svg("images/cellchat/signaling_pathways_EV_marginal_sc_marginal.svg", width=20, height=16)
ht2
dev.off()
