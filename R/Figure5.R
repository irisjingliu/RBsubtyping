
setwd("/Users/jing/Downloads")

#' Plot Figure 5
#' @export
Figure5 <- function(){
  #### Warning: Seurat V2 is required to reproduce the results!

  #### install Seurat 2.3.4 in R 4.0
  # SDMTools required: https://cran.r-project.org/src/contrib/Archive/SDMTools/SDMTools_1.1-221.2.tar.gz
  # install.packages("SDMTools_1.1-221.2.tar.gz", repos = NULL, type = "source")
  # source("https://z.umn.edu/archived-seurat") # install Seurat 2.3.4 hosted by Satija Lab

  #### Loading packages
  library(dplyr)
  library(Seurat)
  library(biomaRt)
  library(clustree)
  library(clusterProfiler)
  library(monocle)
  library(infercnv)

  #### Directories
  dir_data.singlecell <- "/Users/jing/Google Drive/Data/sc_expression_data_rawcount"
  dir_res <- "Fig3"
  dir.create(dir_res)


  #####################################
  #######  START: import data  ########
  #####################################

  #### Import data
  # RBSC11
  setwd(dir_data.singlecell)
  list.files()
  dir_RB_SC11 <- dir_data.singlecell
  RB_SC11 = Read10X(dir_RB_SC11)

  ## first look
  # define cutoff for min.cells (keep genes which expresed in more than n cells)
  cutoff_min.cells <- 1
  # define cutoff for min.genes (keep cells which expresed more than n genes)
  cutoff_min.genes <- 1

  srobj_RB_SC11 <- CreateSeuratObject(raw.data = RB_SC11, project = "RB_SC11", min.cells = cutoff_min.cells, min.genes = cutoff_min.genes)
  srobj_RB_SC11

  #####################################
  #######   END: import data   ########
  #####################################



  ##########################################
  #######  START: QC and Filtering  ########
  ##########################################

  ## first filter
  # define cutoff for min.cells (keep genes which expresed in more than n cells)
  cutoff_min.cells <- 3
  # define cutoff for min.genes (keep cells which expresed more than n genes)
  cutoff_min.genes <- 100

  srobj_RB_SC11 <- CreateSeuratObject(raw.data = RB_SC11, project = "RB_SC11", min.cells = cutoff_min.cells, min.genes = cutoff_min.genes)
  srobj_RB_SC11

  #### QC and filtering parameters
  setwd(dir_res)

  srobj <- srobj_RB_SC11
  name.srobj <- "RB_SC11"

  ## other filters
  cutoff_nGenes <- 500
  cutoff_prct.mito <- 0.05
  cutoff_prct.ribo <- 0.6

  #### QC
  ## Set up control genes
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = srobj@data), value = TRUE)
  percent.mito <- Matrix::colSums(srobj@raw.data[mito.genes, ])/Matrix::colSums(srobj@raw.data)
  table(percent.mito > cutoff_prct.mito)
  srobj <- AddMetaData(object = srobj, metadata = percent.mito, col.name = "percent.mito")

  ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = srobj@data), value = TRUE)
  # ribo.genes <- grep(pattern = "^RP", x = rownames(x = srobj@data), value = TRUE)
  percent.ribo <- Matrix::colSums(srobj@raw.data[ribo.genes, ])/Matrix::colSums(srobj@raw.data)
  table(percent.ribo > cutoff_prct.ribo)
  srobj <- AddMetaData(object = srobj, metadata = percent.ribo, col.name = "percent.ribo")

  table(srobj@meta.data$nGene < cutoff_nGenes)

  #
  summary(srobj@meta.data)
  capture.output(summary(srobj@meta.data), file="0_summary_all.txt")

  ##
  pdf(file="0_QC_all.pdf", paper="a4r", width=12, height=7)

  quantile(srobj@meta.data$nUMI)
  par(mfrow=c(2,2), mar=c(2,3,6,1))
  barplot(sort(srobj@meta.data$nUMI),
          main=paste("Number of UMI/Cells\n\n", capture.output(quantile(srobj@meta.data$nUMI))[1], "\n", capture.output(quantile(srobj@meta.data$nUMI))[2])
          #main="Number of UMI/Cells\n",
          #xlab=paste(capture.output(quantile(srobj@meta.data$nGene))[1], "\n", capture.output(quantile(srobj@meta.data$nGene))[2])
  )

  quantile(srobj@meta.data$nGene)
  barplot(sort(srobj@meta.data$nGene),
          main=paste("Number of genes detected/Cells\n\n", capture.output(quantile(srobj@meta.data$nGene))[1], "\n", capture.output(quantile(srobj@meta.data$nGene))[2])
          #main="Number of genes detected/Cells",
          #xlab=paste(capture.output(quantile(srobj@meta.data$nGene))[1], "\n", capture.output(quantile(srobj@meta.data$nGene))[2])
  )

  quantile(srobj@meta.data$percent.mito)
  barplot(sort(srobj@meta.data$percent.mito),
          main=paste("Percent of mitocondria genes\n\n", capture.output(quantile(srobj@meta.data$percent.mito))[1], "\n", capture.output(quantile(srobj@meta.data$percent.mito))[2])
          #main="Percent of mitocondria genes",
          #xlab=paste(capture.output(quantile(srobj@meta.data$percent.mito))[1], "\n", capture.output(quantile(srobj@meta.data$percent.mito))[2])
  )

  quantile(srobj@meta.data$percent.ribo)
  barplot(sort(srobj@meta.data$percent.ribo),
          main=paste("Percent of ribosome genes\n\n", capture.output(quantile(srobj@meta.data$percent.ribo))[1], "\n", capture.output(quantile(srobj@meta.data$percent.ribo))[2])
          #main="Percent of ribosome genes",
          #xlab=paste(capture.output(quantile(srobj@meta.data$percent.ribo))[1], "\n", capture.output(quantile(srobj@meta.data$percent.ribo))[2])
  )


  par(mfrow = c(1, 1))
  VlnPlot(object = srobj, features.plot = c("nUMI", "nGene", "percent.mito", "percent.ribo"), nCol = 4)

  par(mfrow = c(2, 3), mar=c(5,5,5,1))
  GenePlot(object = srobj, gene1 = "nUMI", gene2 = "nGene")
  GenePlot(object = srobj, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = srobj, gene1 = "nUMI", gene2 = "percent.ribo")
  GenePlot(object = srobj, gene1 = "nGene", gene2 = "percent.mito")
  GenePlot(object = srobj, gene1 = "nGene", gene2 = "percent.ribo")
  GenePlot(object = srobj, gene1 = "percent.mito", gene2 = "percent.ribo")

  dev.off()


  #### filter cells
  srobj <- FilterCells(object = srobj, subset.names = c("nGene", "percent.mito", "percent.ribo"),
                       low.thresholds = c(cutoff_nGenes, -Inf, -Inf), high.thresholds = c(Inf, cutoff_prct.mito, cutoff_prct.ribo))

  # filtering criteria
  capture.output(print("cutoff_nGenes"), cutoff_nGenes,
                 print("cutoff_prct.mito"), cutoff_prct.mito,
                 print("cutoff_prct.ribo"), cutoff_prct.ribo,
                 print("srobj@meta.data$nGene < cutoff_nGenes"), rownames(srobj@meta.data)[which(srobj@meta.data$nGene < cutoff_nGenes)],
                 print("percent.mito > cutoff_prct.mito"), names(which(percent.mito > cutoff_prct.mito)),
                 print("percent.ribo > cutoff_prct.ribo"), which(percent.ribo > cutoff_prct.ribo),
                 file="0_summary_filtering_criteria.txt")

  # summary after filtering
  summary(srobj@meta.data)
  capture.output(summary(srobj@meta.data), file="0_summary_filtered.txt")

  ##
  pdf(file="0_QC_filtered.pdf", paper="a4r", width=12, height=7)

  quantile(srobj@meta.data$nUMI)
  par(mfrow=c(1,2), mar=c(2,3,6,1))
  barplot(sort(srobj@meta.data$nUMI),
          main=paste("After filtering: Number of UMI/Cells\n\n", capture.output(quantile(srobj@meta.data$nUMI))[1], "\n", capture.output(quantile(srobj@meta.data$nUMI))[2])
          #main="After filtering: Number of UMI/Cells\n",
          #xlab=paste(capture.output(quantile(srobj@meta.data$nGene))[1], "\n", capture.output(quantile(srobj@meta.data$nGene))[2])
  )

  quantile(srobj@meta.data$nGene)
  barplot(sort(srobj@meta.data$nGene),
          main=paste("After filtering: Number of genes detected/Cells\n\n", capture.output(quantile(srobj@meta.data$nGene))[1], "\n", capture.output(quantile(srobj@meta.data$nGene))[2])
          #main="After filtering: Number of genes detected/Cells",
          #xlab=paste(capture.output(quantile(srobj@meta.data$nGene))[1], "\n", capture.output(quantile(srobj@meta.data$nGene))[2])
  )

  quantile(srobj@meta.data$percent.mito)
  barplot(sort(srobj@meta.data$percent.mito),
          main=paste("After filtering: Percent of mitocondria genes\n\n", capture.output(quantile(srobj@meta.data$percent.mito))[1], "\n", capture.output(quantile(srobj@meta.data$percent.mito))[2])
          #main="After filtering: Percent of mitocondria genes",
          #xlab=paste(capture.output(quantile(srobj@meta.data$percent.mito))[1], "\n", capture.output(quantile(srobj@meta.data$percent.mito))[2])
  )

  quantile(srobj@meta.data$percent.ribo)
  barplot(sort(srobj@meta.data$percent.ribo),
          main=paste("After filtering: Percent of ribosome genes\n\n", capture.output(quantile(srobj@meta.data$percent.ribo))[1], "\n", capture.output(quantile(srobj@meta.data$percent.ribo))[2])
          #main="Percent of ribosome genes",
          #xlab=paste(capture.output(quantile(srobj@meta.data$percent.ribo))[1], "\n", capture.output(quantile(srobj@meta.data$percent.ribo))[2])
  )

  par(mfrow = c(1, 1))
  VlnPlot(object = srobj, features.plot = c("nGene", "nUMI", "percent.mito", "percent.ribo"), nCol = 4)

  par(mfrow = c(2, 3), mar=c(5,5,5,1))
  GenePlot(object = srobj, gene1 = "nUMI", gene2 = "nGene")
  GenePlot(object = srobj, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = srobj, gene1 = "nUMI", gene2 = "percent.ribo")
  GenePlot(object = srobj, gene1 = "nGene", gene2 = "percent.mito")
  GenePlot(object = srobj, gene1 = "nGene", gene2 = "percent.ribo")
  GenePlot(object = srobj, gene1 = "percent.mito", gene2 = "percent.ribo")

  dev.off()

  ##########################################
  #######   END: QC and Filtering   ########
  ##########################################




  #####################################################
  #######  START: Normalization, Scaling, PCA  ########
  #####################################################

  #### Normalization using Seurat
  median(srobj@meta.data$nUMI)
  # srobj <- NormalizeData(object = srobj, normalization.method = "LogNormalize", scale.factor = 10000) default is 10000
  srobj <- NormalizeData(srobj, scale.factor = median(srobj@meta.data$nUMI))
  #boxplot(as.matrix(srobj@data),las=2,cex.axis=0.5)


  #### Detection of variable genes across the single cells
  par(mfrow = c(1, 1))
  srobj <- FindVariableGenes(object = srobj, mean.function = ExpMean, dispersion.function = LogVMR,
                             x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 0.5)

  length(x = srobj@var.genes)
  dev.print(png, paste("0_var.genes_", length(x = srobj@var.genes), ".png", sep=""), res=300, height=3000, width=3000)

  #### Scale data
  srobj <- ScaleData(object = srobj, vars.to.regress = "nUMI")

  #### PCA
  srobj <- RunPCA(object = srobj, pc.genes = srobj@var.genes, pcs.compute = 100,
                  do.print = TRUE, pcs.print = 1:5, genes.print = 5)

  dev.new()
  VizPCA(object = srobj, pcs.use = 1:21)
  dev.print(png, "1_PCA_gene.png", res=300, height=10000, width=2100)

  PCAPlot(object = srobj, dim.1 = 1, dim.2 = 2)
  dev.print(png, "1_PCA_plot.png", res=300, height=2100, width=2400)

  #srobj <- ProjectPCA(object = srobj, do.print = FALSE)
  #PCHeatmap(object = srobj, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
  PCHeatmap(object = srobj, pc.use = 1:30, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
  dev.print(png, "1_PCA_heatmap.png", res=300, height=9000, width=2400)

  #
  PCA_table <- GetGeneLoadings(object = srobj, reduction.type = "pca", dims.use = 1:100) %>% as.data.frame()
  head(PCA_table)

  #
  mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  filter_gene <- rownames(PCA_table)
  #attributes_biomart <- c("external_gene_name", "description", "chromosome_name", "band", "gene_biotype", "phenotype_description")
  attributes_biomart <- c("external_gene_name", "description", "chromosome_name", "band")
  out_biomart <- getBM(attributes = attributes_biomart,
                       filters = "external_gene_name",
                       values = filter_gene,
                       mart = mart)

  genes.duplicated <- unique(out_biomart$external_gene_name[which(duplicated(out_biomart$external_gene_name))])
  out_biomart.nd <- out_biomart[which(!(out_biomart$external_gene_name %in% genes.duplicated)), ]

  out_biomart.d <- out_biomart[which(out_biomart$external_gene_name %in% genes.duplicated), ]
  out_biomart.d_rm <- out_biomart.d[-grep("^CHR_", out_biomart.d$chromosome_name), ]

  dim(out_biomart.nd)
  dim(out_biomart.d_rm)

  out_biomart.final <- rbind.data.frame(out_biomart.nd, out_biomart.d_rm)

  #
  out_PCA_table <- merge(x=out_biomart.final, y=PCA_table, by.x="external_gene_name", by.y = "row.names", all.y=T)
  colnames(out_PCA_table)
  out_PCA_table <- out_PCA_table %>% arrange(desc(abs(PC1)), desc(abs(PC2)))
  write.table(out_PCA_table, file="1_PCA_table_gene.csv", sep=";", dec=",", row.names = F)

  #
  PCElbowPlot(object = srobj, num.pc = 100)
  dev.print(png, "1_PCA_plot_elbowplot.png", res=300, height=1500, width=3000)

  #####################################################
  #######   END: Normalization, Scaling, PCA   ########
  #####################################################



  #####################################
  ########  START: Clustering  ########
  #####################################

  #### Clustering
  #
  nPC <- 20

  #### FindClusters by different resolutions
  resolutions.i <- seq(0.4, 1.4, 0.1)
  for (res.i in resolutions.i) {
    srobj <- FindClusters(object = srobj, reduction.type = "pca", dims.use = 1:nPC,
                          resolution = res.i, print.output = 0, save.SNN = TRUE, force.recalc = T)
  }

  #### Count nCells for each Cluster
  list_nCells_Cluster <- lapply(1:length(resolutions.i), function(zz) table(srobj@meta.data[, paste("res.", resolutions.i[zz], sep="")]) )
  names(list_nCells_Cluster) <- paste("res.", resolutions.i, sep="")
  list_nCells_Cluster
  capture.output( list_nCells_Cluster, file="2_nCells_Cluster.txt" )

  ## prepare clustreeobj for clustree use
  clustreeobj <- srobj@meta.data

  #### Do ClusTree plot
  clustree(clustreeobj, prefix = "res.")
  dev.print(png, paste("2_clustree_treeplot.png", sep=""), res=300, height=2700, width=2700)


  #### Calulate T-SNE dimention reduction
  srobj <- RunTSNE(object = srobj, dims.use = 1:nPC, do.fast = TRUE)
  PrintFindClustersParams(object = srobj)
  PrintTSNEParams(srobj)

  colnames(srobj@meta.data)
  TSNEPlot(object = srobj, group.by = "res.0.6")
  TSNEPlot(object = srobj, group.by = "res.0.9")
  TSNEPlot(object = srobj, group.by = "res.1.2")

  #
  clustreeobj <- cbind.data.frame(clustreeobj, srobj@dr$tsne@cell.embeddings)
  colnames(clustreeobj)
  clustree_overlay(clustreeobj, prefix = "res.", x_value = "tSNE_1", y_value = "tSNE_2", label_nodes = TRUE)
  dev.print(png, "2_clustree_tSNE_overlay.png", res=300, height=3000, width=3600)

  resolutions.i <- seq(0.4, 1.4, 0.1)
  pdf("2_tSNE_allres.pdf", width=8.5, height=8)
  for (res.i in resolutions.i) {
    TSNEPlot(object = srobj, group.by = paste("res.", res.i, sep=""), plot.title = paste("res.", res.i, sep="") )
    #dev.print(png, paste("2_tSNE_res.", res.i, ".png", sep=""), res=300, height=2100, width=2400 )
  }
  dev.off()

  ####
  res.i = 0.6
  srobj <- FindClusters(object = srobj, reduction.type = "pca", dims.use = 1:nPC, plot.SNN = F,
                        resolution = res.i, print.output = 0, save.SNN = TRUE, force.recalc = T)
  srobj@ident

  #### save RDS
  save(srobj, file = paste(name.srobj, "_nPC", nPC,".RData", sep="") )

  ####
  TSNEPlot(object = srobj)
  dev.print(png, paste("2_tSNE_res.", res.i, ".png", sep=""), res=300, height=2100, width=2400 )

  ####
  srobj.markers <- FindAllMarkers(object = srobj, test.use = "wilcox",
                                  only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  save(srobj.markers, file = paste("srobj.markers_", name.srobj, "_res.", res.i, ".RData", sep="") )


  ####
  cluster.markers <- srobj.markers # %>% group_by(cluster)

  #
  mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  filter_gene <- rownames(cluster.markers)
  attributes_biomart <- c("external_gene_name", "description", "chromosome_name", "band")
  out_biomart <- getBM(attributes = attributes_biomart,
                       filters = "external_gene_name",
                       values = filter_gene,
                       mart = mart)

  # remove duplicated lines and chr start with CHR_
  genes.duplicated <- unique(out_biomart$external_gene_name[which(duplicated(out_biomart$external_gene_name))])
  out_biomart.nd <- out_biomart[which(!(out_biomart$external_gene_name %in% genes.duplicated)), ]

  out_biomart.d <- out_biomart[which(out_biomart$external_gene_name %in% genes.duplicated), ]
  out_biomart.d_rm <- out_biomart.d[-grep("^CHR_", out_biomart.d$chromosome_name), ]

  dim(out_biomart.nd)
  dim(out_biomart.d_rm)

  out_biomart.final <- rbind.data.frame(out_biomart.nd, out_biomart.d_rm)

  #
  out_cluster.markers <- merge(x=out_biomart.final, y=cluster.markers, by.x="external_gene_name", by.y = "gene", all.y=T) %>%
    arrange(cluster, desc(avg_logFC))
  colnames(out_cluster.markers)[c(1, 3)] <- c("gene", "chr")
  write.table(out_cluster.markers, file = paste("3_Cluster.markers_Table_res.", res.i, ".csv", sep=""),
              sep=";", row.names = F, na="")


  #### ClusterProfiler
  # MSigDB Collection
  setwd(dir_msigdb)
  list.files(path = dir_msigdb)
  msigdb = read.gmt("msigdb.v6.2.symbols.gmt.txt")

  msigdb.c2.cp = read.gmt("c2.cp.v6.2.symbols.gmt.txt")
  names_msigdb.c2.cp <- unique(as.vector(msigdb.c2.cp$ont))
  # names_msigdb.c2.cp

  ##
  setwd(dir_res)
  srobj.markers_split <- srobj.markers %>% group_split(cluster)
  enricher_input <- lapply(srobj.markers_split, function(zz) zz$gene)
  enrich_list <- lapply(enricher_input, function(zz) enricher(zz, TERM2GENE=msigdb, minGSSize = 3, universe = rownames(srobj@data) ) )

  ## All collection in MSigDB
  res_enrich_list <- lapply(1:length(enrich_list), function(zz) {
    if(!is.null(enrich_list[[zz]])){
      return(enrich_list[[zz]]@result)
    }
    else if (is.null(enrich_list[[zz]])){
      return(NULL)} } )
  # tmp <- res_enrich_list[[2]]

  ## GO
  res_enrich_GO_list <- lapply(1:length(res_enrich_list), function(zz) {
    if(!is.null(res_enrich_list[[zz]])){
      res_enrich_list[[zz]][grep("^GO_", rownames(res_enrich_list[[zz]])), ] %>% na.omit%>% arrange(p.adjust, desc(Count))
    }
    else if (is.null(res_enrich_list[[zz]])){
      return(NULL)} } )
  # tmp <- res_enrich_GO_list[[1]]

  ## HALLMARK
  res_enrich_HALLMARK_list <- lapply(1:length(res_enrich_list), function(zz) {
    if(!is.null(res_enrich_list[[zz]])){
      res_enrich_list[[zz]][grep("^HALLMARK_", rownames(res_enrich_list[[zz]])), ] %>% na.omit%>% arrange(p.adjust, desc(Count))
    }
    else if (is.null(res_enrich_list[[zz]])){
      return(NULL)} } )
  # tmp <- res_enrich_HALLMARK_list[[1]]

  ## C2CP
  res_enrich_C2CP_list <- lapply(1:length(res_enrich_list), function(zz) {
    if(!is.null(res_enrich_list[[zz]])){
      res_enrich_list[[zz]][which(rownames(res_enrich_list[[zz]]) %in% names_msigdb.c2.cp), ] %>% na.omit%>% arrange(p.adjust, desc(Count))
    }
    else if (is.null(res_enrich_list[[zz]])){
      return(NULL)} } )
  # tmp <- res_enrich_C2CP_list[[1]]

  names_enricher_list <- paste("C", levels(srobj.markers$cluster), sep="")

  #### save
  setwd(dir_res)
  dir.create("3_ClusterProfiler")
  setwd("3_ClusterProfiler")
  dir.create( paste("res.", res.i, sep="") )
  setwd( paste("res.", res.i, sep="") )

  lapply(1:length(res_enrich_list),
         function(zz) {
           if (!is.null(res_enrich_list[[zz]])){
             write.table(res_enrich_list[[zz]][, -grep("Description", colnames(res_enrich_list[[zz]]))],
                         file = paste("Enrichment_", names_enricher_list[zz], "__msigdb_All.csv", sep=""), sep=";", row.names = F)
           } } )

  lapply(1:length(res_enrich_GO_list),
         function(zz) {
           if (!is.null(res_enrich_GO_list[[zz]])){
             write.table(res_enrich_GO_list[[zz]][, -grep("Description", colnames(res_enrich_GO_list[[zz]]))],
                         file = paste("Enrichment_", names_enricher_list[zz], "__msigdb_GO.csv", sep=""), sep=";", row.names = F)
           } } )

  lapply(1:length(res_enrich_HALLMARK_list),
         function(zz) {
           if (!is.null(res_enrich_HALLMARK_list[[zz]])){
             write.table(res_enrich_HALLMARK_list[[zz]][, -grep("Description", colnames(res_enrich_HALLMARK_list[[zz]]))],
                         file = paste("Enrichment_", names_enricher_list[zz], "__msigdb_HALLMARK.csv", sep=""), sep=";", row.names = F)
           } } )

  lapply(1:length(res_enrich_C2CP_list),
         function(zz) {
           if (!is.null(res_enrich_C2CP_list[[zz]])){
             write.table(res_enrich_C2CP_list[[zz]][, -grep("Description", colnames(res_enrich_C2CP_list[[zz]]))],
                         file = paste("Enrichment_", names_enricher_list[zz], "__msigdb_C2CP.csv", sep=""), sep=";", row.names = F)
           } } )


  #### Retinal markers
  setwd(dir_data)

  ####
  retinal_markers_DO <- read.table("retinal_markers.csv", sep=";", header = T, stringsAsFactors = F)
  retinal_markers_H <- read.table("retinal.markers_Hoshino.csv", sep=";", header = T, stringsAsFactors = F)

  retinal_markers_DO <- retinal_markers_DO[-which(retinal_markers_DO$Consensus == ""), ]

  unique(c(retinal_markers_H$Cells, retinal_markers_DO$Consensus))
  unique(retinal_markers_H$Cells)
  retinal_markers_DO$Consensus <- paste(gsub(" ", "_", retinal_markers_DO$Consensus), "DO", sep="_")
  retinal_markers_H$Cells <- paste(retinal_markers_H$Cells, "H", sep="_")

  colnames(retinal_markers_DO)
  colnames(retinal_markers_DO)[3] <- "Cells"
  colnames(retinal_markers_H)
  retinal_markers <- rbind.data.frame(retinal_markers_DO[, c("Cells", "Genes")], retinal_markers_H[, c("Cells", "Genes")])
  colnames(retinal_markers) <- c("ont", "gene")
  retinal_markers$ont <- as.factor(retinal_markers$ont)
  retinal_markers$gene <- toupper(retinal_markers$gene)
  retinal_markers

  ####
  # tmp = enricher(A_O, TERM2GENE=retinal_markers, universe = N)
  # tmp@result

  enrich_retinal_list <- lapply(enricher_input, function(zz) enricher(zz, TERM2GENE=retinal_markers, minGSSize = 2, universe = rownames(srobj@data)) )
  # tmp <- !is.null(res_enrich_retinal_list[[4]])
  res_enrich_retinal_list <- lapply(1:length(enrich_retinal_list), function(zz) {
    if(!is.null(enrich_retinal_list[[zz]])){
      return(enrich_retinal_list[[zz]]@result)
    }
    else if (is.null(enrich_retinal_list[[zz]])){
      return(NULL)} })

  #### save
  setwd(dir_res)
  setwd("3_ClusterProfiler")
  setwd( paste("res.", res.i, sep="") )

  lapply(1:length(res_enrich_retinal_list),
         function(zz) {
           if (!is.null(res_enrich_retinal_list[[zz]])){
             write.table(res_enrich_retinal_list[[zz]][, -grep("Description", colnames(res_enrich_retinal_list[[zz]]))],
                         file = paste("Enrichment_", names_enricher_list[zz], "__RetinalMarkers.csv", sep=""), sep=";", row.names = F)
           } } )



  #### Heatmap
  setwd(dir_res)
  DoHeatmap(object = srobj, genes.use = cluster.markers$gene, slim.col.label = TRUE, remove.key = TRUE)
  dev.print(png, paste("4_Heatmap_cluster.markers_res.", res.i, ".png", sep=""),
            res=300, width=2400, height = nrow(cluster.markers)/10 * 300 + 300)

  top10 <- srobj.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
  DoHeatmap(object = srobj, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
  dev.print(png, paste("4_Heatmap_cluster.markers_res.", res.i, "_top10.png", sep=""),
            res=300, width=2400, height = nrow(top10)/10 * 300 + 300)



  #####################################
  ########   END: Clustering   ########
  #####################################




  ###################################
  ########  START: inferCNV  ########
  ###################################

  setwd(dir_res)
  dir.create("infercnv")
  setwd("infercnv")


  ####
  setwd(dir_res)
  load("RB_SC11_postPCA.RData")

  nPC <- 20

  res.i = 0.6
  FindClusters(object = srobj.merged, reduction.type = "pca", dims.use = 1:nPC,
               resolution = res.i, print.output = 0, save.SNN = TRUE, force.recalc = T)


  ####
  #### save data for infercnv use
  dir.create("infercnv_input")
  infercnv_annotation <- cbind(rownames(srobj@meta.data), paste0("C", srobj@meta.data[, paste("res.", res.i, sep="")], sep="") )
  write.table(infercnv_annotation, file="infercnv_input/annotations_file.txt", quote=F, sep="\t", col.names = F, row.names = F)


  infercnv_refgroup <- c("C5", "C6") # immune cells in this tumor sample as reference

  srobj@raw.data[1:50, 1:50]

  range(srobj@raw.data)
  range(srobj@data)

  #### data matrix of organoids (data public)
  setwd(dir_RO_Public)
  matrix_RO_public <- read.table("GSE119893_rawCounts.csv", sep=",", header = T, row.names = 1)
  matrix_RO_public[1:50, 1:50]
  matrix_RO_public. <- as.sparse(matrix_RO_public)
  range(matrix_RO_public)
  tmp = unlist(matrix_RO_public)
  summary(tmp)

  matrix_RO_public_log <- log(matrix_RO_public + 1)
  matrix_RO_public_log[1:50, 1:50]
  range(matrix_RO_public_log)

  tmp =  as.sparse(matrix_RO_public)
  srobj_RO_public <- CreateSeuratObject(raw.data = matrix_RO_public, min.cells = 3, min.genes = 100)
  srobj_RO_public <- NormalizeData(srobj_RO_public, scale.factor = median(srobj_RO_public@meta.data$nUMI) )
  range(srobj_RO_public@data)

  srobj.infercnv <- AddSamples(srobj, srobj_RO_public@data, add.cell.id = "RO_public", do.normalize = F)
  which(srobj.infercnv@data["ARR3",] != 0)
  head(srobj.infercnv@meta.data)
  tail(srobj.infercnv@meta.data)

  srobj.infercnv@meta.data$res.0.6[grep("^RO", rownames(srobj.infercnv@meta.data))] <- gsub(".*day", "RO_day", rownames(srobj.infercnv@meta.data)[grep("^RO", rownames(srobj.infercnv@meta.data))] )
  srobj.infercnv@meta.data$res.0.6[which(!grepl("^RO", srobj.infercnv@meta.data$res.0.6))] <- paste("RB_C", srobj.infercnv@meta.data$res.0.6[which(!grepl("^RO", srobj.infercnv@meta.data$res.0.6))], sep="")

  #
  levels(as.factor(infercnv_annotation[,2]))
  infercnv_refgroup <- c("RO_day60", "RO_day90", "RO_day200")


  #### datamatrix from retinal ganglion cells

  dir_RO_ganglion <- "/Volumes/LaCie/Work/SingleCellSeq/Publicdata/2018_RetinalOrganoidsGanglionCells_SData/Preprocess/1_cellranger/v0_180626/results/Aggregate/outs/filtered_gene_bc_matrices_mex/GRCh38/"
  RO_ganglion = Read10X(dir_RO_ganglion)

  srobj_RO_ganglion <- CreateSeuratObject(raw.data = RO_ganglion, project = "RO_ganglion", min.cells = 3, min.genes = 100)

  median(srobj_RO_ganglion@meta.data$nUMI)
  srobj_RO_ganglion <- NormalizeData(srobj_RO_ganglion, scale.factor = median(srobj_RO_ganglion@meta.data$nUMI))
  range(srobj_RO_ganglion@data)

  #
  srobj.infercnv <- AddSamples(srobj, srobj_RO_ganglion@data, add.cell.id = "RO_ganglion", do.normalize = F)
  head(srobj.infercnv@meta.data)
  tail(srobj.infercnv@meta.data)


  srobj.infercnv@meta.data$res.0.6[grep("^RO", rownames(srobj.infercnv@meta.data))] <- gsub("^RO_.*-", "RO_ganglion_", rownames(srobj.infercnv@meta.data)[grep("^RO", rownames(srobj.infercnv@meta.data))] )
  srobj.infercnv@meta.data$res.0.6[which(!grepl("^RO", srobj.infercnv@meta.data$res.0.6))] <- paste("RB_C", srobj.infercnv@meta.data$res.0.6[which(!grepl("^RO", srobj.infercnv@meta.data$res.0.6))], sep="")

  #### save data for infercnv use
  setwd("~/Downloads/RO_ganglion")
  dir.create("infercnv_input", recursive = T)
  infercnv_annotation <- cbind(rownames(srobj.infercnv@meta.data), srobj.infercnv@meta.data$res.0.6 )
  write.table(infercnv_annotation, file="infercnv_input/annotations_file.txt", quote=F, sep="\t", col.names = F, row.names = F)

  #
  levels(as.factor(infercnv_annotation[,2]))
  infercnv_refgroup <- c("RO_ganglion_1", "RO_ganglion_2")
  dir_infercnv_output <- "infercnv_output"

  infercnv_refgroup <- c("RO_ganglion_2")
  dir_infercnv_output <- "infercnv_output_ref_THY1neg"



  #### datamatrix PBMC 4k from 10x genomics

  dir_PBMC_4k <- "/Volumes/LaCie/Work/SingleCellSeq/Publicdata/10xGenomics/Chromium_v2_Chemistry/PBMC_4k/filtered_gene_bc_matrices/GRCh38/"
  PBMC_4k = Read10X(dir_PBMC_4k)

  srobj_PBMC_4k <- CreateSeuratObject(raw.data = PBMC_4k, project = "PBMC_4k", min.cells = 3, min.genes = 100)

  median(srobj_PBMC_4k@meta.data$nUMI)
  srobj_PBMC_4k <- NormalizeData(srobj_PBMC_4k, scale.factor = median(srobj_PBMC_4k@meta.data$nUMI))
  range(srobj_PBMC_4k@data)

  #
  srobj.infercnv <- AddSamples(srobj, srobj_PBMC_4k@data, add.cell.id = "PBMC_4k", do.normalize = F)
  head(srobj.infercnv@meta.data)
  tail(srobj.infercnv@meta.data)

  srobj.infercnv@meta.data$res.0.6[which(!grepl("^PBMC", srobj.infercnv@meta.data$res.0.6))] <- paste("RB_C", srobj.infercnv@meta.data$res.0.6[which(!grepl("^RO", srobj.infercnv@meta.data$res.0.6))], sep="")

  # save data for infercnv use

  setwd("~/Downloads/")
  dir.create("PBMC_4k")
  setwd("PBMC_4k")
  dir.create("infercnv_input", recursive = T)
  infercnv_annotation <- cbind(rownames(srobj.infercnv@meta.data), srobj.infercnv@meta.data$res.0.6 )
  write.table(infercnv_annotation, file="infercnv_input/annotations_file.txt", quote=F, sep="\t", col.names = F, row.names = F)

  #
  levels(as.factor(infercnv_annotation[,2]))
  infercnv_refgroup <- c("PBMC_4k")
  dir_infercnv_output <- "infercnv_output_ref_PBMC_4k"




  #### create infercnv object
  cnvobj = CreateInfercnvObject(raw_counts_matrix = as.matrix(srobj.infercnv@data),
                                annotations_file = "infercnv_input/annotations_file.txt",
                                delim = "\t",
                                gene_order_file = paste(dir_data.infercnv, "gencode_v19_gene_pos.txt", sep=""),
                                ref_group_names = infercnv_refgroup)

  # perform infercnv operations to reveal cnv signal
  cnvobj = run(cnvobj,
               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
               out_dir=dir_infercnv_output,  # dir is auto-created for storing outputs
               cluster_by_groups=T,   # cluster
               include.spike=F)

  save(cnvobj, file= "cnvobj.RData")
  # load("cnvobj.RData")

  ###################################
  ########   END: inferCNV   ########
  ###################################


  ########################
  ########################
  # setwd(dir_res)
  capture.output(srobj@calc.params, file="calc.params.txt")
  capture.output(Sys.time(), print("sessionInfo"), sessionInfo(), file = "sessionInfo.txt")




}


