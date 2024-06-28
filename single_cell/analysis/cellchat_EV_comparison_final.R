# Comparison analysis, but using a subsetted dataframe (hence alternate)

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


# CellChat Analysis
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

# Subsetting cellchat database
interaction_db <- CellChatDB$interaction

df_all <- read.csv('data/all_pm.csv')
all_genes <- c(df_all$T..GeneID)

contains_any_gene <- function(string, genes) {
  any(sapply(genes, function(gene) str_detect(string, paste0("\\b", gene, "\\b"))))
}
# Apply this function to each row of the dataframe and create a logical vector
matches <- apply(interaction_db, 1, function(row) {
  contains_any_gene(row["ligand"], all_genes)
})
# Subset the dataframe
subset_interactions <- interaction_db[matches, ]
CellChatDB.use$interaction <- subset_interactions


# Running cellchat
cc_ac <- createCellChat(object = howe_ac, group.by = "ident", assay = "SCT")
cc_ac@DB <- CellChatDB.use

cellchat_ac_a <- subsetData(cc_ac)
cellchat_ac_a <- identifyOverExpressedGenes(cellchat_ac_a)
cellchat_ac_a <- identifyOverExpressedInteractions(cellchat_ac_a)

cellchat_ac_a <- computeCommunProb(cellchat_ac_a, type = "triMean")
cellchat_ac_a <- filterCommunication(cellchat_ac_a, min.cells = 10)
cellchat_ac_a <- computeCommunProbPathway(cellchat_ac_a)
cellchat_ac_a <- aggregateNet(cellchat_ac_a)
cellchat_ac_a <- netAnalysis_computeCentrality(cellchat_ac_a, slot.name = "netP")

cc_pa <- createCellChat(object = howe_pa, group.by = "ident", assay = "SCT")
cc_pa@DB <- CellChatDB.use

cellchat_pa_a <- subsetData(cc_pa)
cellchat_pa_a <- identifyOverExpressedGenes(cellchat_pa_a)
cellchat_pa_a <- identifyOverExpressedInteractions(cellchat_pa_a)

cellchat_pa_a <- computeCommunProb(cellchat_pa_a, type = "triMean")
cellchat_pa_a <- filterCommunication(cellchat_pa_a, min.cells = 10)
cellchat_pa_a <- computeCommunProbPathway(cellchat_pa_a)
cellchat_pa_a <- aggregateNet(cellchat_pa_a)
cellchat_pa_a <- netAnalysis_computeCentrality(cellchat_pa_a, slot.name = "netP")

saveRDS(cellchat_ac_a, file = "cellchat_ac_alternate.rds")
saveRDS(cellchat_pa_a, file = "cellchat_pa_alternate.rds")


library(CellChat)
library(Seurat)
library(reshape2)
library(patchwork)
library(stringr)
library(RColorBrewer)
library(ComplexHeatmap)
library(dplyr)


# Analysis
cellchat_ac_a <- readRDS(file = "objects/cellchat_ac_alternate.rds")
cellchat_pa_a <- readRDS(file = "objects/cellchat_pa_alternate.rds")

object.list_a_a <- list(AC = cellchat_ac_a, PA = cellchat_pa_a)
cellchat_merged_a <- mergeCellChat(object.list_a_a, 
                                           add.names = names(object.list_a_a))

# Number of Interactions
netVisual_heatmap(cellchat_merged_a,color.heatmap = c("#2166ac","#b2182b"),
                  font.size = 14,font.size.title = 16)

compareInteractions(cellchat_merged_a, show.legend = F, group = c(1,2),)


count_interactions_ac <- as.data.frame(cellchat_ac_a@net[["count"]])
count_interactions_ac <- melt(count_interactions_ac)
count_interactions_ac <- count_interactions_ac %>%
  group_by(variable) %>%
  summarize(Sum = sum(value))
count_interactions_ac$source <- 'AC'

count_interactions_pa <- as.data.frame(cellchat_pa_a@net[["count"]])
count_interactions_pa <- melt(count_interactions_pa)
count_interactions_pa <- count_interactions_pa %>%
  group_by(variable) %>%
  summarize(Sum = sum(value))
count_interactions_pa$source <- 'PA'

count_interactions <- rbind(count_interactions_ac, count_interactions_pa)

svg("images/cellchat/intnumber_EV_all_comparison_.svg", units="in", width=12, height=9, res=300)
ggplot(count_interactions, aes(x=variable, y=Sum, fill=source))+
  geom_bar(stat="identity", color="black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

svg("images/cellchat/signaling_pathways_EV_all_comparison_.svg", units="in", width=20, height=16, res=300)
rankNet(cellchat_merged_a, mode = "comparison", stacked = T, do.stat = TRUE,
        font.size = 14,thresh=0.001)
dev.off()


# Compare outgoing/incoming signaling patterns between two datasets
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list_a_a[[i]]@netP$pathways, object.list_a_a[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list_a_a[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list_a_a)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list_a_a[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list_a_a)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# Aggregates incoming/outgoing signalling
ht1 = netAnalysis_signalingRole_heatmap(object.list_a_a[[i]], pattern = "all", signaling = pathway.union, title = names(object.list_a_a)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list_a_a[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list_a_a)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
my_col <- setNames(c("#8FBEDB", "#D24D29", "#D2D1E6", "#D5E9F5", "#FAF6C7", 
                     "#FDD167", "#322767", "#8987BD", "#BAD3EC", "#F0EA37",
                     "#D9EBD2", "#CC2841", "#8A1837", "#2D80BC", "#6AC3B8",
                     "#589ECE", "#563B93", "#C5DFBA", "#9BD0A3", "#055EA8"),
                   c("DN_T_1", "EC_1", "VSMC_1", "CD8T_1", "Macrophage_1",
                     "Macrophage_2", "VSMC_FB", "VSMC_2", "CD8T_2", "Dendritic",
                     "B_Cell", "EC_2", "EC_3", "NKT", "Unknown_immune",
                     "DN_T_2", "VSMC_3", "Plasma_cell", "Mast_cell", "CD8_TRM"))

svg("images/cellchat/signaling_pathways_EV_all_pm.svg", units="in", width=15, height=10, res=300)
ht1 = netAnalysis_signalingRole_heatmap(object.list_a_a[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list_a_a)[i], width = 15, height = 18, color.heatmap = "GnBu",color.use=my_col)
ht2 = netAnalysis_signalingRole_heatmap(object.list_a_a[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list_a_a)[i+1], width = 15, height = 18, color.heatmap = "GnBu",color.use=my_col)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


########################################################
library(circlize)
testchord <- function(object, slot.name = "net", color.use = NULL, color_labels = NULL,
                      signaling = NULL, pairLR.use = NULL, net = NULL,
                      sources.use = NULL, targets.use = NULL,
                      lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                      link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                      transparency = 0.4, link.border = NA,
                      title.name = NULL, legend.pos.x = 20, legend.pos.y = 20, show.legend = TRUE,
                      thresh = 0.05, custom_colors = "default",
                      ...){
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
      slot.name = "netP"
    }
  }
  
  if (!is.null(pairLR.use) & !is.null(signaling)) {
    stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  }
  
  if (is.null(net)) {
    prob <- slot(object, "net")$prob
    pval <- slot(object, "net")$pval
    prob[pval > thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source","target","interaction_name")
    
    pairLR = dplyr::select(object@LR$LRsig, c("interaction_name_2", "pathway_name",  "ligand",  "receptor" ,"annotation","evidence"))
    idx <- match(net$interaction_name, rownames(pairLR))
    temp <- pairLR[idx,]
    net <- cbind(net, temp)
  }
  
  if (!is.null(signaling)) {
    pairLR.use <- data.frame()
    for (i in 1:length(signaling)) {
      pairLR.use.i <- searchPair(signaling = signaling[i], pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
      pairLR.use <- rbind(pairLR.use, pairLR.use.i)
    }
  }
  
  if (!is.null(pairLR.use)){
    if ("interaction_name" %in% colnames(pairLR.use)) {
      net <- subset(net,interaction_name %in% pairLR.use$interaction_name)
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      net <- subset(net, pathway_name %in% as.character(pairLR.use$pathway_name))
    }
  }
  
  if (slot.name == "netP") {
    net <- dplyr::select(net, c("source","target","pathway_name","prob"))
    net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
    net <- net %>% dplyr::group_by(source_target, pathway_name) %>% dplyr::summarize(prob = sum(prob))
    a <- stringr::str_split(net$source_target, "sourceTotarget", simplify = T)
    net$source <- as.character(a[, 1])
    net$target <- as.character(a[, 2])
    net$ligand <- net$pathway_name
    net$receptor <- " "
  }
  
  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- levels(object@idents)[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  } else {
    sources.use <- levels(object@idents)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- levels(object@idents)[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  } else {
    targets.use <- levels(object@idents)
  }
  # remove the interactions with zero values
  df <- subset(net, prob > 0)
  
  if (nrow(df) == 0) {
    stop("No signaling links are inferred! ")
  }
  
  if (length(unique(net$ligand)) == 1) {
    message("You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway")
  }
  
  df$id <- 1:nrow(df)
  # deal with duplicated sector names
  ligand.uni <- unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i <- df[df$ligand == ligand.uni[i], ]
    source.uni <- unique(df.i$source)
    for (j in 1:length(source.uni)) {
      df.i.j <- df.i[df.i$source == source.uni[j], ]
      df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
    }
  }
  receptor.uni <- unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i <- df[df$receptor == receptor.uni[i], ]
    target.uni <- unique(df.i$target)
    for (j in 1:length(target.uni)) {
      df.i.j <- df.i[df.i$target == target.uni[j], ]
      df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
    }
  }
  
  cell.order.sources <- levels(object@idents)[levels(object@idents) %in% sources.use]
  cell.order.targets <- levels(object@idents)[levels(object@idents) %in% targets.use]
  
  df$source <- factor(df$source, levels = cell.order.sources)
  df$target <- factor(df$target, levels = cell.order.targets)
  # df.ordered.source <- df[with(df, order(source, target, -prob)), ]
  # df.ordered.target <- df[with(df, order(target, source, -prob)), ]
  df.ordered.source <- df[with(df, order(source, -prob)), ]
  df.ordered.target <- df[with(df, order(target, -prob)), ]
  
  df.ordered.source$source <- 'EV'
  df.ordered.source <- df.ordered.source %>%
    group_by(source, target, interaction_name) %>%
    summarise(
      prob = sum(prob),  # Sum of 'prob' for each group
      interaction_name_2 = first(interaction_name_2),
      pathway_name = first(pathway_name),
      ligand = first(ligand),
      receptor = first(receptor),
      annotation = first(annotation),
      evidence = first(evidence),
      id = first(id),
      .groups = 'drop'  
    )
  
  order.source <- unique(df.ordered.source[ ,c('ligand','source')])
  order.target <- unique(df.ordered.target[ ,c('receptor','target')])
  
  # define sector order
  order.sector <- c(order.source$ligand, order.target$receptor)
  
  # Remove duplicates for vector order
  order_duplicates <- function(vec) {
    occurrences <- table(vec) # Table to count occurrences of each element
    result_vec <- vec # Initialize the result vector
    # Loop through unique elements that have duplicates
    for (element in names(occurrences[occurrences > 1])) {
      # Get indices of the current element
      indices <- which(vec == element)
      # Append spaces to make each occurrence unique
      result_vec[indices] <- paste0(element, sapply(1:length(indices), function(i) paste(rep(" ", i - 1), collapse = "")))
    }
    return(result_vec)
  }
  order.sector <- order_duplicates(order.sector)
  
  # Pre-defined colour palettes
  color_palettes <- list(
    EC = c('EV' = "#B3B3B3", 'EC_1' = "#D24D29", 'EC_2' = "#CC2841", 'EC_3' = "#8A1837"),
    MC = c('EV' = "#B3B3B3", 'Macrophage_1' = "#FAF6C7", 'Macrophage_2' = "#FDD167"),
    VSMC = c('EV' = "#B3B3B3", 'VSMC_1' = "#D2D1E6", 'VSMC_2' = "#8987BD", 'VSMC_3' = "#563B93"),
    default = scPalette(length(color_labels))
  )
  
  # define cell type color
  if (!is.null(custom_colors) && custom_colors %in% names(color_palettes)) {
    color.use <- color_palettes[[custom_colors]]
  } else {
    color.use <- color_palettes[["default"]]
  }
  names(color.use) <- color_labels
  color.use <- color.use[color_labels]
  
  # define edge color
  edge.color <- color.use[as.character(df.ordered.source$source)]
  names(edge.color) <- as.character(df.ordered.source$source)
  
  # define grid colors
  grid.col.ligand <- color.use[as.character(order.source$source)]
  names(grid.col.ligand) <- as.character(order.source$source)
  grid.col.receptor <- color.use[as.character(order.target$target)]
  names(grid.col.receptor) <- as.character(order.target$target)
  grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
  names(grid.col) <- order.sector
  
  df.plot <- df.ordered.source[ ,c('ligand','receptor','prob')]
  
  # Remove duplicates in plotting dataframe
  df.plot <- df.plot %>%
    mutate(receptor = if_else(ligand == receptor, paste0(receptor, " "), receptor))
  
  if (directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }
  circos.clear()
  chordDiagram(df.plot,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight+arrows"),
               link.arr.type = link.arr.type,
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = scale,
               link.target.prop = link.target.prop,
               reduce = reduce,
               ...)
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
  }, bg.border = NA)
  
  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(color.use), type = "grid", legend_gp = grid::gpar(fill = color.use), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }
  
  circos.clear()
  if(!is.null(title.name)){
    text(-0, 1.02, title.name, cex=1)
  }
  gg <- recordPlot()
  return(gg)
}

svg("images/cellchat/chord_EV_all_upregulated_plaque_EC.svg", units="in", width=10, height=8, res=300)
testchord(object.list_a_a[[1]], targets.use = c('EC_1','EC_2','EC_3'),
                     color_labels = c('EV','EC_1','EC_2','EC_3'),
                     slot.name = 'net', net = net.up, lab.cex = 1,legend.pos.y = 30, 
                     small.gap = 3.5,
                     reduce =0.015, custom_colors = 'EC',
                     title.name = paste0("Up-regulated signaling in ", 
                                         names(object.list_a_a)[1]))
dev.off()

svg("images/cellchat/chord_EV_all_upregulated_marginal_EC.svg", units="in", width=10, height=8, res=300)
testchord(object.list_a_a[[2]], targets.use = c('EC_1','EC_2','EC_3'),
                     color_labels = c('EV','EC_1','EC_2','EC_3'),
                     slot.name = 'net', net = net.down, lab.cex = 1,legend.pos.y = 30, 
                     small.gap = 3.5,
                     reduce = 0.015, custom_colors = 'EC',
                     title.name = paste0("Up-regulated signaling in ", 
                                         names(object.list_a_a)[2]))
dev.off()

svg("images/cellchat/chord_EV_all_upregulated_plaque_VSMC.svg", units="in", width=10, height=8, res=300)
testchord(object.list_a_a[[1]], targets.use = c('VSMC_1','VSMC_2','VSMC_3'),
          color_labels = c('EV','VSMC_1','VSMC_2','VSMC_3'),
          slot.name = 'net', net = net.up, lab.cex = 1,legend.pos.y = 30, 
          small.gap = 3.5,
          reduce =0.015, custom_colors = 'VSMC',
          title.name = paste0("Up-regulated signaling in ", 
                              names(object.list_a_a)[1]))
dev.off()

svg("images/cellchat/chord_EV_all_upregulated_marginal_VSMC.svg", units="in", width=10, height=8, res=300)
testchord(object.list_a_a[[2]], targets.use = c('VSMC_1','VSMC_2','VSMC_3'),
          color_labels = c('EV','VSMC_1','VSMC_2','VSMC_3'),
          slot.name = 'net', net = net.down, lab.cex = 1,legend.pos.y = 30, 
          small.gap = 3.5,
          reduce = 0.015, custom_colors = 'VSMC',
          title.name = paste0("Up-regulated signaling in ", 
                              names(object.list_a_a)[2]))
dev.off()

svg("images/cellchat/chord_EV_all_upregulated_plaque_Macrophage.svg", units="in", width=10, height=8, res=300)
testchord(object.list_a_a[[1]], targets.use = c('Macrophage_1','Macrophage_2'),
          color_labels = c('EV','Macrophage_1','Macrophage_2'),
          slot.name = 'net', net = net.up, lab.cex = 1,legend.pos.y = 30, 
          small.gap = 3.5,
          reduce =0.015, custom_colors ='MC',
          title.name = paste0("Up-regulated signaling in ", 
                              names(object.list_a_a)[1]))
dev.off()

svg("images/cellchat/chord_EV_all_upregulated_marginal_Macrophage.svg", units="in", width=10, height=8, res=300)
testchord(object.list_a_a[[2]], targets.use = c('Macrophage_1','Macrophage_2'),
          color_labels = c('EV','Macrophage_1','Macrophage_2'),
          slot.name = 'net', net = net.down, lab.cex = 1,legend.pos.y = 30, 
          small.gap = 3.5,
          reduce = 0.015, custom_colors = 'MC',
          title.name = paste0("Up-regulated signaling in ", 
                              names(object.list_a_a)[2]))
dev.off()
