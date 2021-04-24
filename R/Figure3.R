
#' Figure 3a limma analysis to generate top_table_1
Figure3_limma <- function(expRBRET){
  require(dplyr)
  require(limma)
  
  gp1 <- intersect(TableS1_FinalSubtype$sample_ID[which(TableS1_FinalSubtype$final_subtype == 1)], colnames(expRBRET))
  gp2 <- intersect(TableS1_FinalSubtype$sample_ID[which(TableS1_FinalSubtype$final_subtype == 2)], colnames(expRBRET))
  
  exp.i = expRBRET[, c(gp1, gp2)]
  
  design <- cbind(gp1=c(rep(1,length(gp1)),rep(0,length(gp2))),gp2=c(rep(0,length(gp1)),rep(1,length(gp2))))
  colnames(design) <- c("group1","group2")
  rownames(design) <- c(gp1,gp2)
  
  # data matrix
  data_i <- cbind(exp.i[, gp1], exp.i[, gp2])
  # limma
  fit <- lmFit(data_i,design)
  cont.matrix <- makeContrasts(contrasts="group2-group1", levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  # top table
  top_table_1 <- topTable(fit2, number=Inf, adjust="BH",p.value=1)
  
  return(top_table_1)
}
  

#' Plot Figure 3a volcano plot of differential expressed genes between 2 subtypes
#' @export
Figure3a_volcanoPlot <- function(expRBRET){
  require(dplyr)
  require(limma)

  top_table_1 <- Figure3_limma(expRBRET)
  
  # 
  gp1genes = c("EGF", "TPBG", "GUCA1C", "GUCA1B", "GUCA1A", "GNAT2", "GNGT2", "ARR3", "PDE6C", "PDE6H", "OPN1SW")
  gp2genes = c("TFF1", "CD24", "EBF3", "GAP43", "STMN2", "POU4F2", "SOX11", "EBF1", "DCX", "ROBO1", "PCDHB10")
  
  col <- rep("black", nrow(top_table_1))
  col[which(top_table_1$adj.P.Val<0.05 & top_table_1$logFC<0)]="goldenrod"
  col[which(top_table_1$adj.P.Val<0.05 & top_table_1$logFC>0)]="cornflowerblue"
  
  names(col) <- rownames(top_table_1)
  
  tiff(filename = "Figure3a_volcanoPlot.tiff",
       width = 720, height = 450, units = "px", pointsize = 12)
  
  par(family = "Arial")
  par(mar = c(5.5, 5, 1, 1) + 0.1)
  par(mgp = c(3.5,1,0))
  #par(cex.names = 1.2)
  #dev.new()
  plot(top_table_1$logFC, -log10(top_table_1$adj.P.Val), xlim=c(-7,7), ylim=c(0,13.5), 
       pch=20, cex=0.4, frame=FALSE, col=col, font.lab=2, cex.lab = 1.3, 
       ylab = "-log10 (adjusted p-value)", xlab = "mRNA fold change \n(log2 scale, mixed-type vs cone-like)")
  
  legends = c(gp1genes, gp2genes)
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  x = top_table_1[c(gp1genes, gp2genes), ]$logFC
  y = -log10(top_table_1[c(gp1genes, gp2genes), ]$adj.P.Val)
  names(x) <- c(gp1genes, gp2genes)
  names(y) <- c(gp1genes, gp2genes)
  
  y["GUCA1B"] = 6.6 ## GUCA1B
  y["GNAT2"] = 6.1 ## GNAT2
  y["ARR3"] = 5.7 ## ARR3
  y["GNGT2"] = 5.2 ## GNGT2
  y["STMN2"] = 4.25## STMN2
  #y[14] = 4.12 ## GUCA1A
  #x[14] = -1.8 ## GUCA1A
  
  y["PDE6C"] = 3.45 ## PDE6C
  y["PDE6H"] = 3.07 ## PDE6H
  y["OPN1SW"] = 2.69 ## OPN1SW
  
  y["POU4F2"] = 3.8 ## POU4F2
  
  names = c(gp1genes, gp2genes)
  positions = ifelse(col[c(gp1genes, gp2genes)] == "goldenrod", 2, 4)
  
  
  collabel = col
  col.legends = col[legends]
  col.legends = ifelse(col.legends == "goldenrod", "goldenrod4", "darkblue")
  col.legends
  col.legends["EGF"] = "black"
  col.legends["TPBG"] = "black"
  
  col.legends["TFF1"] = "black"
  col.legends["CD24"] = "black"
  
  for (i in legends){
    text(x = x[i], y = y[i], labels = i, col=col.legends[legends][i], cex = 1.2, pos=positions[i], font = 3, adj = 0.5)   ## color 
    #text(x = x[i], y = y[i], labels = i, col="black", cex = 1.2, pos=positions[i], font = 3, adj = 0.5)   ## color black
    # font = 3 italic; 2 bold; 1 normal
  }
  
  dev.off()

}


#' Plot Figure3b Heatmap
#' @export
Figure3b_heatmap <- function(expRBRET){
  require(dplyr)
  require(limma)
  require(ComplexHeatmap)
  
  top_table_1 <- Figure3_limma(expRBRET)
  sample.i_ordered <- c("RB28", "RB1", "RB49", "RB6", "RB2", "RB35", "RB4", "RB25", "RB40", "RB47", #C1
                        "RB216", "RB50", "RB213", "RB11", "RB221", "RB44", "RB33", "RB24", "RB217", "RB9", #C1
                        "RB219", "RB23", "RB34", "RB10", "RB63",
                        "RB52", #C1
                        "RB51", #C2
                        "RB42", "RB224", "RB14", "RB212", "RB222", "RB22", "RB46", "RB211", "RB21", "RB32", #C2
                        "RB43", "RB37", "RB41", "RB15", "RB31", "RB27", "RB38", "RB39", "RB7", "RB48", #C2
                        "RB45", "RB59", "RB54", "RB57", "RB220", "RB58", "RB62", "RB225", "RB223", "RB56" #C2
  )
  exp.i <- expRBRET[, sample.i_ordered]
  genes.i <- rownames(top_table_1)[which(top_table_1$adj.P.Val < 0.05)]
  
  exp.ic <- exp.i[genes.i, ] - rowMeans(exp.i[genes.i, ]) # centering
  
  
  # dendrogram
  distance = "pearson"
  linkage = "average"
  # 1-pearson correlation(pairs of genes) as distance
  d_gene = as.dist(1-cor(t(exp.ic), method=distance))
  # hierarchical clustering
  hc_gene = hclust(d_gene, method=linkage)
  hc_gene

  # 1-pearson correlation(pairs of genes) as distance
  d_sample = as.dist(1-cor(exp.ic, method=distance))
  # hierarchical clustering
  hc_sample = hclust(d_sample, method=linkage)
  hc_sample
  # plot(hc_sample)
  # 
  # annot column

  df_col <- data.frame( Subtype = paste("Subtype ",  TableS1_FinalSubtype[colnames(exp.ic), "final_subtype"] , sep=""))
  rownames(df_col) <- colnames(exp.ic)
  ha_col <- HeatmapAnnotation(df = df_col, show_legend = F,
                              col = list(Subtype = c("Subtype 1" = "goldenrod", "Subtype 2" = "cornflowerblue")),
                              show_annotation_name = F )
  
  
  table(cutree(hc_gene, k = 15))
  df_row_cluster <- as.data.frame(cutree(hc_gene, k = 15))
  df_row_cluster <- df_row_cluster[rownames(exp.ic), , drop=F]
  colnames(df_row_cluster) <- "Cluster"
  df_row_cluster$Cluster[which(df_row_cluster$Cluster == 2)] <- 1.1
  df_row_cluster$Cluster[which(df_row_cluster$Cluster == 3)] <- 1.2
  table(df_row_cluster)
  df_row_cluster[intersect(rownames(df_row_cluster), rownames(top_table_1)[which(top_table_1$logFC > 0)]), 1]<- 2
  df_row_cluster[which(!(df_row_cluster$Cluster %in% c(1.1, 1.2, 2))), 1] <- NA
  df_row_cluster$Cluster <- as.factor(df_row_cluster$Cluster)
  
  ha_row_cluster <- rowAnnotation(df = df_row_cluster,
                                  na_col = "white",
                                  show_annotation_name = T,
                                  width = unit( ncol(df_row_cluster) , "cm"))

  #
  H <- Heatmap(exp.ic, name="exp",
               # column_title = pname,
               # cluster_columns = hc_sample,
               show_row_dend = T,
               column_order = colnames(exp.ic),
               cluster_rows = hc_gene,
               # split = cluster_toplot,
               top_annotation = ha_col,
               show_row_names = FALSE
  )
  H + ha_row_cluster
  dev.print(png, "Figure2b_Heatmap.png", res = 300, height = 2000, width = 2700)
  
}







