#' Plot Figure 4
#' @export
Figure4 <- function(){

  #### figure 3 study oct 2018
  #### Jing
  #### 2019.07.10
  #### 1st version Oct 15 2018
  #### updated in Oct 16 2018
  #### updated in Oct 30 2018
  #### updated in Nov 05 2018
  #### updated in Nov 08 2018
  #### updated in Nov 12 2018
  #### updated in July 10 2018

  #### objective: make nice final plots

  #### update Oct 30:
  ## use ComplexHeatmap to plot 2 heatmaps
  ## change sample order by put RB674 and RB704 next to c1

  #### update Nov 05:
  ## keep only 4 ganglion genes POU4F2, EBF3, EBF1, GAP43
  ## make heatmap width smaller (retinoblastoma part)
  ## add fold change in addition to significance

  #### update Nov 08:
  ## trial: change organoids color to grey and grey scale

  #### update Nov 12:
  ## A split cone ganglion genes by a blank line
  ## B change correlation plot color to green scale
  ## C modify grey scale to reach a better contrast

  #### update Nov 12:
  ## A order: OTX2, CRX, THRB, RXRG

  ####

  library(gtools)
  library(dplyr)
  library(ggpubr)
  library(gplots)
  library(ape)
  library(corrplot)
  library(GLAD)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)


  #### directories
  dir_data <- "/Volumes/LaCie/Work/Retinoblastoma/00_CancerCell/Figure_NanoString/NanoString_organoids_tumors/V20181016/data"
  dir_res <- "/Volumes/LaCie/Work/Retinoblastoma/00_CancerCell/Figure_NanoString/NanoString_organoids_tumors/V20181016/results/plots_20191112"



  #### prepare data files  ####

  file.copy("/Volumes/LaCie/Work/Retinoblastoma/6_Nanostring/NanoString2017_2batches/sample_description_2batches.csv", dir_data)
  file.copy("/Volumes/LaCie/Work/Database/Retinoblastoma/DataCancerCell/group103.RData", dir_data)
  file.copy("/Volumes/LaCie/Work/Retinoblastoma/6_Nanostring/NanoString2017_2batches/data/nanostring_2017_2batch_remove7sample_norm1_log2.txt", dir_data)
  file.copy("/Volumes/LaCie/Work/Retinoblastoma/1_Transcriptome/CRX_RNAseq/retinal.markers_ordered.csv", dir_data)


  #### load data ####
  setwd(dir_data)
  #### sample description
  sample.description = read.table("sample_description_2batches.csv", sep = ";", header = T, stringsAsFactors = F)
  rownames(sample.description) = sample.description$Description
  #### load group annotation
  load("group103.RData")
  rownames(group103) = group103$Sample

  ## Remove 7 samples - NORM1
  nano = read.table("nanostring_2017_2batch_remove7sample_norm1_log2.txt", sep = "\t", header = T, stringsAsFactors = F)
  nano.des = nano[,1:2]
  nano.exp = nano[1:96,3:121] %>% as.matrix
  rownames(nano.exp) = nano[1:96,1]
  sample.description.i = sample.description[colnames(nano.exp), ]
  rownames(nano.exp) = gsub(" ex ", "ex", rownames(nano.exp))
  rownames(nano.exp) = gsub(" ex", "ex", rownames(nano.exp))
  rownames(nano.exp) = gsub("-", "_", rownames(nano.exp))

  #### retinal marker genes ordered by chronological order
  retinal.markers = read.table("retinal.markers_ordered.csv", sep = ";", header = T, stringsAsFactors = F)
  nano.exp
  genes.retinal.markers = retinal.markers$Gene


  genes_cone <- c("OTX2", "CRX", "THRB", "RXRG", "PDE6H", "GNAT2", "ARR3", "GUCA1C")
  genes_ganglion <- c("POU4F2", "EBF3", "EBF1", "GAP43")
  # genes_ganglion <- c("POU4F2", "EBF3", "EBF1", "GAP43", "ATOH7", "SOX11", "NEFL")

  ####
  #### sample selection for nanostring

  ## organoids
  samples.organoids = sample.description.i[which(sample.description.i$Tissue == "Organoid"), "Description"] %>% mixedsort

  ## celllines
  samples.celllines = sample.description.i[which(sample.description.i$Tissue == "Cellline"), "Description"]
  samples.y79 = samples.celllines[which(startsWith(samples.celllines, "Y79"))]
  samples.weri = samples.celllines[which(startsWith(samples.celllines, "WERI"))] %>% sort
  samples.weri = samples.weri[c(7:9, 1:3, 13:15, 4:6, 10:12)]

  ## tumors
  samples.tumors = sample.description.i[which(sample.description.i$Tissue == "Tumor"), "Description"] %>% mixedsort %>% setdiff("RB102")

  samples.tumors.c1 = intersect(samples.tumors, group103$Sample[which(group103$FinalGroup == 1)]) %>% mixedsort
  samples.tumors.c2 = intersect(samples.tumors, group103$Sample[which(group103$FinalGroup == 2)]) %>% mixedsort
  samples.tumors.c1.c2 = c(samples.tumors.c1, samples.tumors.c2)

  samples.tumors.c1.m = intersect(samples.tumors, group103$Sample[which(group103$FinalGroupMYCN == 1)]) %>% mixedsort
  samples.tumors.c2.m = intersect(samples.tumors, group103$Sample[which(group103$FinalGroupMYCN == 2)]) %>% mixedsort
  samples.tumors.c3.m = intersect(samples.tumors, group103$Sample[which(group103$FinalGroupMYCN == 3)]) %>% mixedsort
  samples.tumors.c1.c2.c3 = c(samples.tumors.c1.m, samples.tumors.c2.m, samples.tumors.c3.m)

  ## retina
  samples.retina = sample.description.i[which(sample.description.i$Tissue == "Retina"), "Description"] %>% mixedsort
  samples.retina = c("RET215_110", "RET2", "RET1", "RET15_92")


  #### prepare mean expression data for different days of organoids and C1/C2 C1/C2woMYCN/MYCN tumors
  nano.exp.mean.organoids = lapply(list("D0", "D35", "D49", "D56", "D84", "D112", "D173"),
                                   function(z) rowMeans(nano.exp[ , samples.organoids[grep(z, samples.organoids)]])) %>% cbind.data.frame() %>% as.matrix()
  colnames(nano.exp.mean.organoids) = c("D0", "D35", "D49", "D56", "D84", "D112", "D173")
  nano.exp.mean.organoids


  nano.exp.mean.tumors = lapply(list(samples.tumors.c1, samples.tumors.c2),
                                function(z) rowMeans(nano.exp[, z] ) ) %>% cbind.data.frame() %>% as.matrix()
  colnames(nano.exp.mean.tumors) = c("C1", "C2")
  nano.exp.mean.tumors


  nano.exp.mean.tumors.mycn = lapply(list(samples.tumors.c1.m, samples.tumors.c2.m, samples.tumors.c3.m),
                                     function(z) rowMeans(nano.exp[, z] ) ) %>% cbind.data.frame() %>% as.matrix()
  colnames(nano.exp.mean.tumors.mycn) = c("C1m", "C2woMYCN", "CMYCN")
  nano.exp.mean.tumors.mycn


  nano.exp.all <- cbind(nano.exp, nano.exp.mean.organoids, nano.exp.mean.tumors, nano.exp.mean.tumors.mycn)


  # organoids + tumors C1 C2
  nano.exp.i = nano.exp[c(genes_cone, genes_ganglion), c(samples.organoids, samples.tumors.c1.c2)]
  nano.exp.i = nano.exp[genes_cone, c(samples.organoids, samples.tumors.c1.c2)]
  nano.exp.i = nano.exp[genes_ganglion, c(samples.organoids, samples.tumors.c1.c2)]


  # organoids + retina + tumors C1 C2
  nano.exp.i = nano.exp[c(genes_cone, genes_ganglion), c(samples.organoids, samples.retina, samples.tumors.c1.c2)]
  nano.exp.i = nano.exp[genes_cone, c(samples.organoids, samples.retina, samples.tumors.c1.c2)]
  nano.exp.i = nano.exp[genes_ganglion, c(samples.organoids, samples.retina, samples.tumors.c1.c2)]


  # organoids + retina + tumors C1 C2 C3
  nano.exp.i = nano.exp[c(genes_cone, genes_ganglion), c(samples.organoids, samples.retina, samples.tumors.c1.c2.c3)]
  nano.exp.i = nano.exp[genes_cone, c(samples.organoids, samples.retina, samples.tumors.c1.c2.c3)]
  nano.exp.i = nano.exp[genes_ganglion, c(samples.organoids, samples.retina, samples.tumors.c1.c2.c3)]


  #### setwd to dir_res
  # dir.create(dir_res)
  setwd(dir_res)


  #### correlation plot
  # colnames(nano.exp.i)
  # corrplot(cor(nano.exp.i[, 1:25], nano.exp.i[, 26:ncol(nano.exp.i)]))
  #
  # corrplot(cor(nano.exp.i.organoids.mean, nano.exp.i[, 22:ncol(nano.exp.i)]))
  #
  # corrplot(cor(nano.exp.i.organoids.mean, cbind(nano.exp.i[, 22:25], nano.exp.i.tumors.mean)))
  #
  # corrplot(cor(nano.exp.mean.tumors, nano.exp.mean.organoids), tl.col = "black", cl.lim = c(0, 1))
  # dev.print(png, "corplot.png", res=300, width= 1600, height=2400)

  # heatcolor.i <- colorRamp2( c(0, 0.2, 0.4, 0.6, 0.8, 1),
  #                            c("#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")) # blue scale
  # heatcolor.i <- colorRamp2( c(0, 0.2, 0.4, 0.6, 0.8, 1),
  #                            c("#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")) # red scale
  brewer.pal(9, "Greens")
  heatcolor.i <- colorRamp2( c(0, 0.2, 0.4, 0.6, 0.8, 1),
                             c("#FFFFFF", "#E5F5E0", "#C7E9C0", "#74C476", "#238B45", "#00441B")) # green scale

  matrix.correlation1 <- cor(nano.exp.mean.organoids[genes_cone,-1], nano.exp.mean.tumors[genes_cone,], method = "pearson")
  matrix.correlation2 <- cor(nano.exp.mean.organoids[genes_cone,-1], nano.exp.mean.tumors[genes_cone,], method = "spearman")
  summary(matrix.correlation1)
  summary(matrix.correlation2)

  capture.output(matrix.correlation1, file="cor_pearson.txt")
  capture.output(matrix.correlation2, file="cor_spearman.txt")


  Heatmap(matrix.correlation1, col = heatcolor.i, cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "black"),
          row_names_side = "left", column_names_side = "top", name = "r")
  dev.print(png, "correlation_pearson2.png", res=300, width= 550, height=750)

  Heatmap(matrix.correlation2, col = heatcolor.i, cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "black"),
          row_names_side = "left", column_names_side = "top", name = "r")
  dev.print(png, "correlation_spearman2.png", res=300, width= 750, height=750)

  corrplot(matrix.correlation, tl.col = "black", tl.srt = 0, cl.lim = c(0, 1))
  dev.print(png, "corplot22.png", res=300, width= 1000, height=1800)

  #
  ?RColorBrewer
  colorRamp(brewer.pal(9, "Greens")(100) )
  ?colorRampPalette
  # setwd("~/Downloads/")
  #### heatmap
  nano.exp.i <- cbind(nano.exp.mean.organoids[c(genes_cone, genes_ganglion),], nano.exp[c(genes_cone, genes_ganglion), c(samples.tumors.c1.c2)])

  nano.exp.i <- nano.exp.all[c(genes_cone, genes_ganglion),
                             c("D0",
                               "D35", "D49", "D56", "D84", "D112", "D173",  ## organoids
                               "RB63", "RB52", "RB603", "RB218", "RB47", "RB49", "RB28", "RB634",
                               "RB25", "RB23", "RB4", "RB9", "RB217", "RB35", "RB213", "RB11", "RB2", "RB30",
                               "RB10", "RB34", "RB219", "RB203", "RB216",  ## C1 end
                               "RB647", "RB704", "RB3", "RB51", "RB42", "RB593", "RB48", "RB21", "RB220", "RB43", "RB54", "RB39", "RB598",
                               "RB38", "RB27", "RB31", "RB716", "RB55", "RB58", "RB56", "RB223", "RB46", "RB45", "RB222", "RBsjd8", "RB200", "RB57",
                               "RB111", "RB215", "RB211", "RB225", "RB590", "RBsjd7", "RB37", "RB205", "RB224", "RB62", "RB59", "RB212", "RB15", "RB7", "RB14", "RB209", "RBsjd2" ## C2 end
                             )]

  #### by ComplexHeatmap

  matrix.organoids <- nano.exp.i[ , c("D35", "D49", "D56", "D84", "D112", "D173") ]
  matrix.C1 <- nano.exp.i[ , c("RB63", "RB52", "RB603", "RB218", "RB47", "RB49", "RB28", "RB634",
                               "RB25", "RB23", "RB4", "RB9", "RB217", "RB35", "RB213", "RB11", "RB2", "RB30",
                               "RB10", "RB34", "RB219", "RB203", "RB216") ]
  matrix.C2 <- nano.exp.i[ , c("RB647", "RB704", "RB3", "RB51", "RB42", "RB593", "RB48", "RB21", "RB220", "RB43", "RB54", "RB39", "RB598",
                               "RB38", "RB27", "RB31", "RB716", "RB55", "RB58", "RB56", "RB223", "RB46", "RB45", "RB222", "RBsjd8", "RB200", "RB57",
                               "RB111", "RB215", "RB211", "RB225", "RB590", "RBsjd7", "RB37", "RB205", "RB224", "RB62", "RB59", "RB212", "RB15", "RB7", "RB14", "RB209", "RBsjd2") ]

  # added nov 12
  # corrplot(cor(nano.exp.i[, samples.organoids], nano.exp.i[, samples.tumors.c1.c2]))
  matrix.correlation1 <- cor(matrix.organoids[genes_cone,], matrix.C1[genes_cone,])
  matrix.correlation2 <- cor(matrix.organoids[genes_cone,], matrix.C2[genes_cone,])

  Heatmap(matrix.correlation1, col = heatcolor.i, cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "black"),
          row_names_side = "left", column_names_side = "top", name = "r")
  dev.print(png, "correlation_pearson1.png", res=300, width= 1300, height=750)

  Heatmap(matrix.correlation2, col = heatcolor.i, cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "black"),
          row_names_side = "left", column_names_side = "top", name = "r")
  dev.print(png, "correlation_pearson2.png", res=300, width= 1500, height=750)

  dev.print(png, "heatmap_correlationfig3.png", res=300, height=1200, width=3900)


  ncol(matrix.organoids)
  ncol(matrix.C1)
  ncol(matrix.C2)

  summary(nano.exp.i)
  max(nano.exp.i)
  min(nano.exp.i)
  mean(nano.exp.i)
  median(nano.exp.i)


  heatcolor.i <- colorRamp2( c( 0, 9, 15 ), c("blue", "white", "red"))
  # brewer.pal(9, "BuGn")[4:9]


  # col.organoids <- gray.colors(n = 6, start = 0.90, end = 0.20, gamma = 2.2, alpha = NULL)
  col.organoids <- rep("white", 6)
  colnames(matrix.organoids)
  df1 <- data.frame( Organoids = colnames(matrix.organoids) )
  ha1 = HeatmapAnnotation(df = df1, annotation_legend_param = list(title="Retinal organoids", nrow = 3),
                          col = list(Organoids = c("D35" = col.organoids[1], "D49" = col.organoids[2], "D56" = col.organoids[3],
                                                   "D84" = col.organoids[4], "D112" = col.organoids[5], "D173" = col.organoids[6]) ) )

  df2 <- data.frame( RetinoblastomaC1 = rep("Subtype 1", ncol(matrix.C1)) )
  ha2 = HeatmapAnnotation(df = df2, annotation_legend_param = list(title="Retinoblastoma"),
                          col = list(RetinoblastomaC1 = c("Subtype 1" = "goldenrod") ) )

  df3 <- data.frame( RetinoblastomaC2 = rep("Subtype 2", ncol(matrix.C2)) )
  ha3 = HeatmapAnnotation(df = df3, annotation_legend_param = list(title=""),
                          col = list(RetinoblastomaC2 = c("Subtype 2" = "cornflowerblue") ) )

  # c(rep("yellow1", 4), rep("yellow3", 4), rep("aquamarine1", 3), rep("aquamarine3", 4))
  rownames(nano.exp.i)
  df_row <- data.frame( marker = c(rep("Early cone", 4), rep("Late cone", 4), rep("Ganglion", 4) ) )
  ha_row <- HeatmapAnnotation(df = df_row, annotation_legend_param = list(title="Retinal markers", nrow = 3),
                              col = list(marker = c("Early cone" = "yellow1", "Late cone" = "yellow3",
                                                    "Ganglion" = "skyblue1") ),
                              which = "row")

  ht1 = Heatmap(matrix.organoids, col = heatcolor.i, cluster_rows = F, cluster_columns = F,
                top_annotation = ha1,
                width = unit(5/23*6 *2, "cm"), column_title = "Retinal organoids",
                show_row_names = F, show_column_names = F, show_heatmap_legend = F)
  ht2 = Heatmap(matrix.C1, col = heatcolor.i, cluster_rows = F, cluster_columns = F,
                top_annotation = ha2,
                width = unit(5 /1.25, "cm"), column_title = "Subtype 1",
                show_row_names = F, show_column_names = F, show_heatmap_legend = F)
  ht3 = Heatmap(matrix.C2, col = heatcolor.i, cluster_rows = F, cluster_columns = F,
                top_annotation = ha3,
                width = unit(5/23*44 /1.25, "cm"), column_title = "Subtype 2",
                show_row_names = T, show_column_names = F, show_heatmap_legend = T, name = "Expression")

  split = rep(c("cone early", "cone late", "ganglion"), each = 4)
  ht_list = ht1 + ht2 + ht3
  draw(ht_list, split = split)

  # ht_final <- ha_row + ht1 + ht2 + ht3
  # draw(ht_final)
  # draw(ht_final, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

  dev.print(png, "heatmap_fig3.png", res=300, height=1200, width=3900)

  #### by heatmap.2

  # no center no scale
  # exp.input <- nano.exp.i

  # center not scale
  # exp.input <- nano.exp.i - rowMeans(nano.exp.i)

  # center and scale
  # exp.input <- t(scale(t(nano.exp.i)))
  # clipper <- 2
  # exp.input[exp.input < -clipper]<- -clipper
  # exp.input[exp.input > clipper] <- clipper

  # ColSideColors <- c(brewer.pal(9, "BuGn")[3:9], rep("goldenrod", length(samples.tumors.c1)), rep("cornflowerblue", length(samples.tumors.c2)))
  # RowSideColors <- c(rep("yellow1", 4), rep("yellow3", 4), rep("aquamarine1", 3), rep("aquamarine3", 4))
  #
  # heatmap.2(exp.input,
  #           #distfun = function(z) as.dist(z),
  #           dendrogram = "none",
  #           Colv = F,
  #           Rowv = F,
  #           col = myPalette(low = "blue", mid = "white", high = "red", k = 100),
  #           ColSideColors = ColSideColors,
  #           RowSideColors = RowSideColors,
  #           key = T,
  #           #density.info="none",
  #           #zlim = c(-2, 2),
  #           trace = "none")
  # dev.print(png, "heatmap2.png", res=200, width= 2400, height=1600)
  #
  # plot(0,0, frame=FALSE, xlim=c(0,110), ylim=c(0,20), axes=FALSE, cex=0.001, xlab="", ylab="")
  # legend("topleft", title="Legend", bty="n",
  #        fill=c(brewer.pal(9, "BuGn")[3:9], "goldenrod", "cornflowerblue", "yellow4", "yellow2", "aquamarine4", "aquamarine2"),
  #        legend=c(paste("Retinal organoids", c("D0", "D35", "D49", "D56", "D84", "D112", "D173"), sep = " "), "Subtype 1", "Subtype 2", "Early Cone","Late Cone","Early ganglion","Late ganglion")
  # )
  # dev.print(png, file="legend.png", res = 300, width=1800,height=3000)



  #
  # plot(0,0, frame=FALSE, xlim=c(0,110), ylim=c(0,20), axes=FALSE, cex=0.001, xlab="", ylab="")
  # legend("topleft",fill=c(subgroup,cell),legend=c("ConeDiff","Mixed","Cone genes","Cone Ganglion gene","Ganglion genes"),title="Legend",bty="n")
  # dev.print(png, file="heatmaplegend.png", res = 300, width=3000,height=3000)
  # dev.off()


  #### phylogenetic tree

  tmp <- nano.exp[genes_cone, c(samples.organoids[-c(1:3)], samples.tumors.c1.c2)]
  tmp <- tmp - rowMeans(tmp)
  d = dist(t(tmp))
  #d = as.dist(1 - cor(tmp))
  phylo = fastme.bal(d)
  plot(phylo, type = "u")


  group_phylo <- group103[phylo$tip.label, 2]
  group_phylo[!is.na(group_phylo)] <- paste("C", group_phylo[!is.na(group_phylo)], sep="")
  # group_phylo[grep("iPs_D0", phylo$tip.label)] <- "D0"
  group_phylo[grep("iPs_D35", phylo$tip.label)] <- "D35"
  group_phylo[grep("iPs_D49", phylo$tip.label)] <- "D49"
  group_phylo[grep("iPs_D56", phylo$tip.label)] <- "D56"
  group_phylo[grep("iPs_D84", phylo$tip.label)] <- "D84"
  group_phylo[grep("iPs_D112", phylo$tip.label)] <- "D112"
  group_phylo[grep("iPs_D173", phylo$tip.label)] <- "D173"
  # group_phylo <- factor(group_phylo, levels = c("D0", "D35", "D49", "D56", "D84", "D112", "D173", "C1", "C2"))
  group_phylo <- factor(group_phylo, levels = c("D35", "D49", "D56", "D84", "D112", "D173", "C1", "C2"))
  group_phylo
  # col = c(brewer.pal(9, "BuGn")[3:9], "goldenrod", "cornflowerblue")[group_phylo]
  # col = c(brewer.pal(9, "BuGn")[4:9], "goldenrod", "cornflowerblue")[group_phylo]
  col.organoids <- gray.colors(n = 16, start = 0.90, end = 0.20, gamma = 2.2, alpha = NULL)[c(1, 5, 9, 13, 15, 16)]
  col.organoids
  col = c(col.organoids, "goldenrod", "cornflowerblue")[group_phylo]


  shape = c(rep(22, 6), rep(21, 2) )[group_phylo]
  #
  phylo2 <- phylo
  phylo2$tip.label = rep("", length(phylo$tip.label))

  plotBreakLongEdges(phylo2, n=2, type="unrooted", no.margin=TRUE, lab4ut="axial", edge.width=2, tip.color = col, rotate.tree = 15)
  # tiplabels(pch=22, col="black", bg=col, cex=3)
  # tiplabels(pch=21, col="black", bg=col, cex=3)
  tiplabels(pch=shape, col="black", bg=col, cex=3)

  dev.print(png,filename = "phylo.organoids.cone.eucledian_break2_rotate15.png", res=300,width=2100*1.8,height=2100*1)



  #############
  p.ttest <- lapply(1:nrow(nano.exp.i), function(z) t.test(nano.exp.i[z, samples.tumors.c1], nano.exp.i[z, samples.tumors.c2])$p.value )
  p.adj <- p.adjust(p.ttest, "BH")
  significance <- rep("ns", length(p.adj))
  significance[which(p.adj < 0.05)] <- "*"
  significance[which(p.adj < 0.01)] <- "**"
  significance[which(p.adj < 0.001)] <- "***"
  significance[which(p.adj < 0.0001)] <- "****"
  logfc_c2vsc1 <- lapply(1:nrow(nano.exp.i), function(z) mean(nano.exp.i[z, samples.tumors.c2]) - mean(nano.exp.i[z, samples.tumors.c1]) )

  write.table(cbind(rownames(nano.exp.i), round(unlist(logfc_c2vsc1), 2), p.ttest, p.adj, significance), "t_test.csv", sep=";", dec=",", row.names = F )





  fc = logfc_c2vsc1 %>% unlist
  range(fc)
  ha_barplot = rowAnnotation(barplot1 = row_anno_barplot(fc, baseline = 0, bar_width = 0.6, axis = T, ylim = c(-4, 4), border = F,
                                                         gp = gpar(fill = ifelse(fc > 0, "cornflowerblue", "goldenrod"))),
                             width = unit(2, "cm"))
  ?anno_barplot


  ht1 = Heatmap(matrix.organoids, col = heatcolor.i, cluster_rows = F, cluster_columns = F,
                top_annotation = ha1,
                width = unit(5/23*6 *2, "cm"), column_title = "Retinal organoids",
                show_row_names = F, show_column_names = T, show_heatmap_legend = F)
  ht2 = Heatmap(matrix.C1, col = heatcolor.i, cluster_rows = F, cluster_columns = F,
                top_annotation = ha2,
                width = unit(5 /1.25, "cm"), column_title = "Subtype 1",
                show_row_names = F, show_column_names = T, show_heatmap_legend = F)
  ht3 = Heatmap(matrix.C2, col = heatcolor.i, cluster_rows = F, cluster_columns = F,
                top_annotation = ha3,
                width = unit(5/23*44 /1.25, "cm"), column_title = "Subtype 2",
                show_row_names = T, show_column_names = T, show_heatmap_legend = T, name = "Expression")
  #
  # split = rep(c("cone early", "cone late", "ganglion"), each = 4)
  ht_list = ht1 + ht2 + ht3
  draw(ha_barplot + ht_list, split = split)
  dev.print(png, "heatmap_fig3_addbar.png", res=300, height=1300, width=4500)



  ####

  matrix.organoids.output <- nano.exp[rownames(matrix.organoids), samples.organoids[-c(1:3)] ]
  # matrix.organoids.output <- nano.exp[rownames(matrix.C1), samples.organoids[-c(1:3)] ]
  supptable_nano <- cbind(matrix.organoids.output, matrix.C1, matrix.C2)
  write.table(supptable_nano, file="supptable_nano.txt", sep="\t", col.names = NA)
  ####
  setwd("/Users/jing/Downloads/")
  write.table(supptable_nano, file="supptable_nano_ganglion_20190902.txt", sep="\t", col.names = NA)


}


