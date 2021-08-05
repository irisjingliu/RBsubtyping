
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



#' Plot Figure3d Stemness
#' @export
Figure3d_stemness <- function(){

  ################################################
  #### plot stem indice

  # load stem indice data
  load("/Volumes/LaCie/Work/Retinoblastoma/8_Other/Stem_Indice/mRNAsi/predict_2020/results/centered/df.i_stem_indice_mRNA_centered_annot.RData")

  ## order sample by stem indice in each subtype
  df.i_ordered = df.i %>% arrange(subtype, stem_indice)
  sample.i_ordered = rownames(df.i_ordered)[which(df.i_ordered$subtype != "NA")]



  load("/Volumes/LaCie/Work/Retinoblastoma/7_SingleCellRNAseq/seurat/v3_20200309/data/list_collection_enricher_20200309.RData")
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
  levels(list_collection_enricher$h.all.v7.0$ont)

  ## HALLMARK pathway metascore by average expression

  ## HALLMARK pathway metascore by ssgsea
  library(GSVA)
  exp.i = exp[, grep("^RB", colnames(exp))]
  dim(exp.i)

  load("/Volumes/LaCie/Work/Retinoblastoma/7_SingleCellRNAseq/seurat/v3_20200309/data/list_collection_enricher_20200309.RData")
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
  paths = levels(list_collection_enricher$h.all.v7.0$ont)

  # gsva(expr = exp, gset.idx.list = list(genes_Miranda), method = "ssgsea")
  list_ssgseascore_hallmark = lapply(paths, function(path.i) gsva(expr = exp.i, gset.idx.list = list( intersect(list_collection_enricher$h.all.v7.0[which(list_collection_enricher$h.all.v7.0$ont == path.i), "gene"], rownames(exp.i) ) ), method = "ssgsea") )
  names(list_ssgseascore_hallmark) = paths
  matrix_ssgseascore_hallmark = do.call(rbind.data.frame, list_ssgseascore_hallmark) %>% as.matrix()
  df_ssgseascore_hallmark = t(matrix_ssgseascore_hallmark) %>% as.data.frame()
  colnames(df_ssgseascore_hallmark) = paste("ssgseascore", gsub("^HALLMARK_", "", paths), sep="_")
  rownames(df_ssgseascore_hallmark) = colnames(exp.i)


  ## HALLMARK pathway by average expression
  list_score_hallmark = lapply(paths, function(path.i) colMeans(exp.i[intersect(list_collection_enricher$h.all.v7.0[which(list_collection_enricher$h.all.v7.0$ont == path.i), "gene"], rownames(exp.i)), ]) )
  names(list_score_hallmark) = paths
  matrix_score_hallmark = do.call(rbind.data.frame, list_score_hallmark) %>% as.matrix()
  df_score_hallmark = t(matrix_score_hallmark) %>% as.data.frame()
  colnames(df_score_hallmark) = paste("metascore", gsub("^HALLMARK_", "", paths), sep="_")
  rownames(df_score_hallmark) = colnames(exp.i)

  ##
  df.ii = Reduce(cbind.data.frame, list(df.i, df_score_hallmark[rownames(df.i) , ], df_ssgseascore_hallmark[rownames(df.i) , ]))
  colnames(df_ssgseascore_hallmark)
  cor_stemindice_ssgseascore <- sapply(colnames(df_ssgseascore_hallmark), function(zz) cor(df.ii$stem_indice, df.ii[, zz], method = "spearman") )
  colnames(df_score_hallmark)
  cor_stemindice_score <- sapply(colnames(df_score_hallmark), function(zz) cor(df.ii$stem_indice, df.ii[, zz], method = "spearman") )

  df_cor = data.frame(cor_stemindice_score, cor_stemindice_ssgseascore)



  colnames(df.ii)
  dftoheatmap = df.ii[sample.i_ordered, c("stem_indice", # "stemness_score_ssgsea_Miranda",
                                          "metascore_E2F_TARGETS", "metascore_MYC_TARGETS_V2", "metascore_MYC_TARGETS_V1",
                                          "ssgseascore_E2F_TARGETS", "ssgseascore_MYC_TARGETS_V2", "ssgseascore_MYC_TARGETS_V1",
                                          "metascore_INTERFERON_ALPHA_RESPONSE", "metascore_INTERFERON_GAMMA_RESPONSE",
                                          "ssgseascore_INTERFERON_ALPHA_RESPONSE", "ssgseascore_INTERFERON_GAMMA_RESPONSE",
                                          "Monocytic_lineage", "B_lineage", "Cytotoxic_lymphocytes")]

  sapply(c("metascore_E2F_TARGETS", "metascore_MYC_TARGETS_V2", "metascore_MYC_TARGETS_V1",
           "metascore_INTERFERON_ALPHA_RESPONSE", "metascore_INTERFERON_GAMMA_RESPONSE",
           "Monocytic_lineage", "B_lineage", "Cytotoxic_lymphocytes"),
         function(zz) cor.test(df.ii[, "stem_indice"], df.ii[, zz])$p.value)



  colnames(dftoheatmap)
  rho = sapply(colnames(dftoheatmap),
               function(zz) cor(dftoheatmap[,"stem_indice"], dftoheatmap[,zz], method = "spearman")
  )


  rho.i = rho[c("metascore_E2F_TARGETS", "metascore_MYC_TARGETS_V2", "metascore_MYC_TARGETS_V1")]
  ha_row.metascoreC2 = HeatmapAnnotation(rho = anno_barplot(rho.i, baseline=0, ylim = c(-1, 1),
                                                            axis = F, border = F,
                                                            gp = gpar(fill = ifelse(rho.i>0, "red", "green"))),
                                         # signif = anno_text(ifelse(p.spearman_stemindice_mcpscore < 0.05, "*", "ns")),
                                         show_annotation_name = F,
                                         width = unit(2, "cm"),
                                         which = "row" )

  rho.i = rho[c("ssgseascore_E2F_TARGETS", "ssgseascore_MYC_TARGETS_V2", "ssgseascore_MYC_TARGETS_V1")]
  ha_row.ssgseascoreC2 = HeatmapAnnotation(rho = anno_barplot(rho.i, baseline=0, ylim = c(-1, 1),
                                                              axis = F, border = F,
                                                              gp = gpar(fill = ifelse(rho.i>0, "red", "green"))),
                                           # signif = anno_text(ifelse(p.spearman_stemindice_mcpscore < 0.05, "*", "ns")),
                                           show_annotation_name = F,
                                           width = unit(2, "cm"),
                                           which = "row" )

  rho.i = rho[c("metascore_INTERFERON_ALPHA_RESPONSE", "metascore_INTERFERON_GAMMA_RESPONSE")]
  ha_row.metascoreC1 = HeatmapAnnotation(rho = anno_barplot(rho.i, baseline=0, ylim = c(-1, 1),
                                                            axis = F, border = F,
                                                            gp = gpar(fill = ifelse(rho.i>0, "red", "green"))),
                                         # signif = anno_text(ifelse(p.spearman_stemindice_mcpscore < 0.05, "*", "ns")),
                                         show_annotation_name = F,
                                         width = unit(2, "cm"),
                                         which = "row" )

  rho.i = rho[c("ssgseascore_INTERFERON_ALPHA_RESPONSE", "ssgseascore_INTERFERON_GAMMA_RESPONSE")]
  ha_row.ssgseascoreC1 = HeatmapAnnotation(rho = anno_barplot(rho.i, baseline=0, ylim = c(-1, 1),
                                                              axis = F, border = F,
                                                              gp = gpar(fill = ifelse(rho.i>0, "red", "green"))),
                                           # signif = anno_text(ifelse(p.spearman_stemindice_mcpscore < 0.05, "*", "ns")),
                                           show_annotation_name = F,
                                           width = unit(2, "cm"),
                                           which = "row" )

  rho.i = rho[c("Monocytic_lineage", "B_lineage", "Cytotoxic_lymphocytes")]
  ha_row.mcpscore = HeatmapAnnotation(rho = anno_barplot(rho.i, baseline=0, ylim = c(-1, 1),
                                                         axis = T, border = F,
                                                         gp = gpar(fill = ifelse(rho.i>0, "red", "green"))),
                                      # signif = anno_text(ifelse(p.spearman_stemindice_mcpscore < 0.05, "*", "ns")),
                                      width = unit(2, "cm"),
                                      which = "row" )



  df_col <- data.frame( Subtype = paste("Subtype ", c(rep(1, 26), rep(2, 31)), sep=""))
  # rownames(df_col) <- colnames(exp.ic)
  ha_col <- HeatmapAnnotation(df = df_col, show_legend = F,
                              col = list(Subtype = c("Subtype 1" = "goldenrod", "Subtype 2" = "cornflowerblue")),
                              show_annotation_name = F,
                              border = T)


  matrix.cd24 = t(as.matrix(exp["CD24", rownames(dftoheatmap)]))
  rownames(matrix.cd24) = "CD24"

  matrix.stem_indice = t(as.matrix(dftoheatmap[, "stem_indice"]))
  rownames(matrix.stem_indice) = "Stemness Indice"

  # matrix.stemness_score = t(as.matrix(dftoheatmap[, "stemness_score_ssgsea_Miranda"]))
  # rownames(matrix.stemness_score) = "Stemness Score\nMiranda et al. 2019"

  matrix.toplot.ssgseascoreC1 = dftoheatmap[, c("ssgseascore_INTERFERON_ALPHA_RESPONSE", "ssgseascore_INTERFERON_GAMMA_RESPONSE")] %>% scale(center=T, scale=F) %>% t()
  rownames(matrix.toplot.ssgseascoreC1) = c("Interferon alpha response", "Interferon gamma response")

  matrix.toplot.metascoreC1 = dftoheatmap[, c("metascore_INTERFERON_ALPHA_RESPONSE", "metascore_INTERFERON_GAMMA_RESPONSE")] %>% scale(center=T, scale=F) %>% t()
  rownames(matrix.toplot.metascoreC1) = c("HM: Interferon alpha response", "HM: Interferon gamma response")

  matrix.toplot.ssgseascoreC2 = dftoheatmap[, c("ssgseascore_E2F_TARGETS", "ssgseascore_MYC_TARGETS_V2", "ssgseascore_MYC_TARGETS_V1")] %>% scale(center=T, scale=F) %>% t()
  rownames(matrix.toplot.ssgseascoreC2) = c("E2F targets", "MYC targets V2", "MYC targets V1")

  matrix.toplot.metascoreC2 = dftoheatmap[, c("metascore_E2F_TARGETS", "metascore_MYC_TARGETS_V2", "metascore_MYC_TARGETS_V1")] %>% scale(center=T, scale=F) %>% t()
  rownames(matrix.toplot.metascoreC2) = c("HM: E2F targets", "HM: MYC targets V2", "HM: MYC targets V1")

  matrix.mcpscore.selected = dftoheatmap[, c("Monocytic_lineage", "B_lineage", "Cytotoxic_lymphocytes")] %>% scale(center=T, scale=F) %>% t()
  rownames(matrix.mcpscore.selected) = c("MCP: Monocytic lineage", "MCP: B lineage", "MCP: Cytotoxic lymphocytes")

  # matrix.3 = dftoheatmap[, c("score_hm_notch", "score_hm_wnt", "score_hm_hedgehog", "score_hm_tnfa", "score_muller_glia")] %>% scale(center=T, scale=T) %>% t()
  # H3 = Heatmap(matrix.3, show_column_names = F, name = "score3", #right_annotation = ha_row1,
  #                   # row_names_side = "left",
  #                   cluster_columns = F, cluster_rows = F, border = T)
  # draw(H_cd24 %v% H_mRNAsi %v% H2 %v% H1 %v% Hmcp %v% H3, heatmap_legend_side = "bottom")
  # ggscatter(dftoheatmap, "score_muller_glia", "score_hm_notch")

  H_cd24 = Heatmap(matrix.cd24, show_column_names = F, name = "CD24",
                   cluster_columns = F, cluster_rows = F, border = T,
                   # row_names_side = "left",
                   top_annotation = ha_col )

  H_mRNAsi = Heatmap(matrix.stem_indice, show_column_names = F, name = "Stem Indice\nMalta et al, Cell, 2018", # row_names_side = "left",
                     top_annotation = ha_col, row_names_side = "left",
                     cluster_columns = F, cluster_rows = F, border = T)

  # H_ss.Miranda = Heatmap(matrix.stemness_score, show_column_names = F, name = "Stemness Score\n(Miranda et al. PNAS, 2019)",
  #                        # row_names_side = "left",
  #                        cluster_columns = F, cluster_rows = F, border = T)
  # top_annotation = ha_col )

  H.metascoreC1 = Heatmap(matrix.toplot.metascoreC1, show_column_names = F, name = "metascoreC1", #right_annotation = ha_row1, # row_names_side = "left",
                          right_annotation = ha_row.metascoreC1, row_names_side = "left",
                          cluster_columns = F, cluster_rows = F, border = T)

  H.metascoreC2 = Heatmap(matrix.toplot.metascoreC2, show_column_names = F, name = "metascoreC2", #right_annotation = ha_row1, # row_names_side = "left",
                          right_annotation = ha_row.metascoreC2, row_names_side = "left",
                          cluster_columns = F, cluster_rows = F, border = T)

  H.ssgseascoreC1 = Heatmap(matrix.toplot.ssgseascoreC1, show_column_names = F, name = "ssgseascoreC1", #right_annotation = ha_row1, # row_names_side = "left",
                            right_annotation = ha_row.ssgseascoreC1,
                            cluster_columns = F, cluster_rows = F, border = T)

  H.ssgseascoreC2 = Heatmap(matrix.toplot.ssgseascoreC2, show_column_names = F, name = "ssgseascoreC2", #right_annotation = ha_row1, # row_names_side = "left",
                            right_annotation = ha_row.ssgseascoreC2,
                            cluster_columns = F, cluster_rows = F, border = T)

  Hmcp = Heatmap(matrix.mcpscore.selected, show_column_names = F, name = "scoremcp", #right_annotation = ha_row2, # row_names_side = "left",
                 right_annotation = ha_row.mcpscore, row_names_side = "left",
                 cluster_columns = F, cluster_rows = F, border = T)


  # setwd(dir_res)
  setwd("/Users/jing/Downloads")
  draw(H_mRNAsi %v% H.metascoreC2 %v% H.metascoreC1 %v% Hmcp, heatmap_legend_side = "bottom")
  dev.print(png, "stemindice_metascore.png", res=300, width=1500, height=1500)

  draw(H_mRNAsi %v% H.ssgseascoreC2 %v% H.ssgseascoreC1 %v% Hmcp, heatmap_legend_side = "bottom")
  dev.print(png, "stemindice_ssgsea.png", res=300, width=1500, height=1500)


  draw(H_mRNAsi %v% H.metascoreC2 %v% H.metascoreC1 %v% Hmcp)
  dev.print(png, "stemindice_metascorewider2.png", res=300, width=2100, height=900)

  draw(H_mRNAsi %v% H.ssgseascoreC2 %v% H.ssgseascoreC1 %v% Hmcp)
  dev.print(png, "stemindice_ssgseawider.png", res=300, width=2000, height=1000)


  ## barplot
  dftoboxplot = df.i[which(df.i$subtype != "NA"), ]
  dftoboxplot$subtype = paste("Subtype", as.vector(dftoboxplot$subtype), sep="") %>% as.factor()
  dftoboxplot$subtypeMYCN = paste("Subtype", as.vector(dftoboxplot$subtypeMYCN), sep="")
  dftoboxplot$subtypeMYCN[which(dftoboxplot$subtypeMYCN == "SubtypeMYCNamp")] = "Subtype2MYCNamplified"
  dftoboxplot$subtypeMYCN = as.factor(dftoboxplot$subtypeMYCN)

  par(margin = c(4,4,4,4))
  ggboxplot(dftoboxplot, y = "stem_indice", ylab = "Stem Indice",
            x = "subtype", xlab = "",
            add = "jitter",
            color = "subtype", palette = c("goldenrod", "cornflowerblue")) + rremove("legend") # + stat_compare_means(vjust = -1) +  theme(plot.margin=unit(c(1,1,1,1),"cm"))

  dev.print(png, "stemindice_barplot.png", res=300, width=1000, height=900)

  ggboxplot(dftoboxplot, y = "stem_indice", ylab = "Stem Indice",
            x = "subtype", xlab = "",
            add = "jitter", add.params = list(color = c("goldenrod", "cornflowerblue", "darkred")[dftoboxplot$subtypeMYCN]),
            color = "subtype", palette = c("goldenrod", "cornflowerblue")
  ) +
    stat_compare_means() +
    rremove("legend")

  dev.print(png, "stemindice_barplot_3color.png", res=300, width=1200, height=900)




}



#' Plot Figure3e Heatmap - cone and ganglion genes
#' @export
Figure3e_heatmap_cone_ganglion <- function(){

  #### plot retinal marker heatmap


  # tumor and ret
  load("/Volumes/LaCie/Work/Retinoblastoma/1_Transcriptome/0_preprocess/Exp_data_RB_RET_CancerCell_20190531/results/annotation_exp_matrix_probe_combat_20408_features_85_samples_RB_RET_CancerCell_U133Plus2_RMA_BA_entrezgV23_20190531.RData")
  exp <- as.matrix(annotation_exp_matrix_probe_combat[, 9:ncol(annotation_exp_matrix_probe_combat)])
  rownames(exp) <- annotation_exp_matrix_probe_combat$symbol
  samples.i <- c("RET215_110 ", "RET1", "RET2", sample.i_ordered)
  samples.i <- c(sample.i_ordered)


  genes.cone <- c("OPN1SW", "CCDC136", "SYP", "XRCC4", "GNB3", "GNAT2", "PDE6H", "GNGT2", "GUCA1C", "ARR3",
                  "RORA", "RXRG", "THRB", "PDE6C", "PEX5L", "GRK7", "TTR", "SLC17A7", "RAMP1", "SALL3",
                  "ONECUT1", "RORB")

  genes.ganglion <- c("POU4F2", "EBF3", "ONECUT1", "SOX11", "ELAVL3", "STMN2", "GAP43", "EBF1")



  genes.i <- genes.cone
  celltype.i <- "Cone cell markers"

  genes.i <- genes.ganglion
  celltype.i <- "Ganglion cell markers"

  toptable$gene <- rownames(toptable)
  toptable.i <- toptable[genes.i, ] %>% arrange(logFC)
  if(celltype.i == "Ganglion cell markers") {toptable.i <- toptable[genes.i, ] %>% arrange(desc(logFC))}

  genes.i <- toptable.i$gene
  rownames(toptable.i) <- genes.i


  exp.i <- exp[genes.i, samples.i]
  exp.ic <- exp.i - rowMeans(exp.i)

  ## text annotation of significance
  pvalues.i <- toptable.i[genes.i, ]$adj.P.Val
  significance <- rep("ns", length(pvalues.i))
  significance[which(pvalues.i < 0.05)] <- "*"
  significance[which(pvalues.i < 0.01)] <- "**"
  significance[which(pvalues.i < 0.001)] <- "***"
  significance[which(pvalues.i < 0.0001)] <- "****"
  significance

  ha_sigif_text <- rowAnnotation(signif = row_anno_text(significance)) # , gp = gpar(fontsize = 1:12+4)))
  ha_rowname_text <- rowAnnotation(rowname = row_anno_text(rownames(exp.ic)) ) # , gp = gpar(fontsize = 1:12+4)))

  fc = toptable.i[genes.i,]$logFC %>% unlist
  range(fc)
  ha_fc_bar = rowAnnotation(barplot1 = row_anno_barplot(fc, baseline = 0, bar_width = 0.6, ylim = c(-3, 3),
                                                        border = F, axis = T,
                                                        gp = gpar(fill = ifelse(fc > 0, "cornflowerblue", "goldenrod"))),
                            show_annotation_name = F,
                            width = unit(2, "cm")
  )


  # tumor and retina
  df_col <- data.frame( Subtype = c(rep("Retina", 3), paste("Subtype ", c(rep(1, 26), rep(2, 31)), sep="")))
  # only tumor
  df_col <- data.frame( Subtype = paste("Subtype ", c(rep(1, 26), rep(2, 31)), sep=""))

  rownames(df_col) <- colnames(exp.ic)

  # tumor and retina
  ha_col <- HeatmapAnnotation(df = df_col, show_legend = F, border = T,
                              col = list(Subtype = c("Retina" = "pink", "Subtype 1" = "goldenrod", "Subtype 2" = "cornflowerblue")),
                              show_annotation_name = F )
  # only tumor
  ha_col <- HeatmapAnnotation(df = df_col, show_legend = F, border = T,
                              col = list(Subtype = c("Subtype 1" = "goldenrod", "Subtype 2" = "cornflowerblue")),
                              show_annotation_name = F )


  H <- Heatmap(exp.ic, name="exp", border = T, column_title = celltype.i,
               cluster_columns = F, cluster_rows = F,
               show_column_names = F, show_row_names = T, # column_title = "retina",
               top_annotation = ha_col,
               heatmap_legend_param = list(direction = "horizontal"))
  draw(H, heatmap_legend_side = "bottom")
  draw(H + ha_rowname_text + ha_sigif_text + ha_fc_bar, heatmap_legend_side = "bottom")
  setwd("/Users/jing/Downloads")
  if(celltype.i == "Cone cell markers"){dev.print(png, file = paste(celltype.i, ".png", sep=""), res = 300, height = 1800, width = 1500)}
  if(celltype.i == "Ganglion cell markers"){dev.print(png, file = paste(celltype.i, ".png", sep=""), res = 300, height = 920, width = 1500)}


  H <- Heatmap(exp.ic, name="exp", border = T, row_title = celltype.i,
               cluster_columns = F, cluster_rows = F,
               show_column_names = F, show_row_names = T, # column_title = "retina",
               top_annotation = ha_col,
               heatmap_legend_param = list(direction = "horizontal"))
  draw(H + ha_rowname_text + ha_sigif_text + ha_fc_bar, heatmap_legend_side = "bottom")

  if(celltype.i == "Cone cell markers"){dev.print(png, file = paste(celltype.i, "wider.png", sep=""), res = 300, height = 1700, width = 1800)}
  if(celltype.i == "Ganglion cell markers"){dev.print(png, file = paste(celltype.i, "wider.png", sep=""), res = 300, height = 920, width = 1800)}


}




