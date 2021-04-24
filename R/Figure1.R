
#########
#' Plot Figure 1a Upper Panel left, Consensus Clustering mRNA
#' @export
Figure1a_ConsensusClustering_mRNA <- function(input_exp = expRBRET){
  require(cit.utils)
  require(cit.unsupervised)

  projectName <- "mRNA"

  expRBRET <- input_exp
  expRB <- expRBRET[, grep("^RB", colnames(expRBRET))]
  exp.i <- expRB - rowMeans( expRB )

  # params for ICA
  nbcomp_ICA = 3
  cutoff_ICA = 2.5

  resJade <- JADE::JADE(exp.i, n.comp = nbcomp_ICA, maxiter = 10000)
  resS <- resJade$S
  # find the component related to normal cell infiltration (Rod and Muller glia markers)
  # resS[c("RHO", "SAG", "RLBP1", "HES1"), ]

  IC.i = names(which.max(abs(resS["RHO", ])))
  if (resS["RHO", IC.i] > 0) {
    genes.i <- rownames(resS)[which(resS[,IC.i] > cutoff_ICA)]
  } else if (resS["RHO", IC.i] < 0) {
    genes.i <- rownames(resS)[which(resS[,IC.i] < -1*cutoff_ICA)]
  }
  genes.i
  resS[genes.i,IC.i][order(resS[genes.i,IC.i])][1:20]

  exp.i <- expRB[setdiff(rownames(expRB), genes.i), ]

  clustering = cit.clusteringAnalysis(projectName, data = as.data.frame(exp.i),
                                      doTests = F, K.clustering=6, meth.Clustering=c("complete", "average", "ward"), rowCentering=TRUE)
  mp <- clustering$multipartition
  mpK2 <- mp[,grep("K=2",names(mp))]

  ch <- cit.getConsHier(multpart =mpK2,K=2)
  conspart <- ch$consensusPartition   # partition consensus
  conshier <- ch$consensusHierarchy    # hierarchie consensus
  coclas <- as.data.frame(cit.countPairs(mpK2))  # matrice de co-classement
  rownames(coclas) <- names(coclas) <- rownames(mpK2)

  palette <- cit.myPalette(low = "white", high = "blue", k = 11)
  persoc2 <- colorRampPalette( c("#A7CEE4","#1F78B4") )
  dev.new(width = 10, height = 10)
  tut=cit.heatmap(as.matrix(coclas)  ,colclust.hclust =conshier,rowclust.hclust =conshier,heatmapcolors = palette,
                  colpart = conspart,  colpart.bgcol =  c("#A7CEE4","#1F78B4") ,colclust.heightlabel = "",
                  collabels.cex = 0.001,collabels.ymax = 0.01,
                  colpart.labels = "",
                  ## suppression des lables samples et label de groupes
                  lw = c(2, 5, 1, 20),left.margin = 0, top.margin = 0,lh = c(7, 1,2, 20),title = projectName)
  print(tut)
  dev.print(png,filename = "Figure1a_upperpanel_ConsensusClustering_mRNA.png", res=360,width=2400,height=2400)



}



#' Plot Figure 1a Upper Panel middle, Consensus Clustering DNA methylation
#' @export
Figure1a_ConsensusClustering_meth <- function(input_meth = meth){
  require(cit.utils)
  require(cit.unsupervised)

  projectName <- "DNA methylation"
  meth.i <- input_meth

  clustering <- cit.clusteringAnalysis(projectName, data = as.data.frame(meth.i), rowCentering=TRUE,
                                       doTests = F, K.clustering=6, meth.Clustering=c("complete", "average", "ward") )
  mp <- clustering$multipartition
  # rm(mp)
  # mp <- get(load( paste(projectName, "-multPart-centered.RData", sep="") ) )$multipartition

  mpK2 <- mp[,grep("K=2",names(mp))]
  ch <- cit.getConsHier(multpart =mpK2,K=2)
  conspart <- ch$consensusPartition   # partition consensus
  conshier <- ch$consensusHierarchy    # hierarchie consensus
  coclas <- as.data.frame(cit.countPairs(mpK2))  # matrice de co-classement
  rownames(coclas) <- names(coclas) <- rownames(mpK2)

  palette <- cit.myPalette(low = "white", high = "blue", k = 11)
  persoc2 <- colorRampPalette( c("#A7CEE4","#1F78B4"))
  dev.new(width = 10, height = 10)
  tut=cit.heatmap(as.matrix(coclas), colclust.hclust = conshier, rowclust.hclust = conshier, colclust.colors = "black", heatmapcolors = palette,
                  colpart = conspart,  colpart.bgcol =  c("#A7CEE4","#1F78B4") ,colclust.heightlabel = "",
                  collabels.cex = 0.001,collabels.ymax = 0.01,
                  colpart.labels = "",
                  ## suppression des lables samples et label de groupes
                  lw = c(2, 5, 1, 20),left.margin = 0, top.margin = 0,lh = c(7, 1,2, 20),title = projectName)
  print(tut)
  dev.print(png,filename = "Figure1a_upperpanel_ConsensusClustering_meth.png", res=360,width=2400,height=2400)

}



#' Plot Figure 1a Upper Panel right, Consensus Clustering CNV # bug
#' @export
Figure1a_ConsensusClustering_CNV <- function(input_CNV){
  require(cit.utils)
  require(cit.unsupervised)

  #### Copy number data ####
  # sample of interest
  load(file="/Users/jing/Google Drive/Data/expression_data.RData")
  load(file="/Users/jing/Google Drive/Data/methylation_data.RData")

  samps <- unique(c( grep("^RB", colnames(expRBRET), value = T) , colnames(meth)))
  samps

  # copy
  load(file="/Users/jing/Google Drive/Data/copynum.RData")

  # samps<-unique(c(colnames(exp), colnames(meth)))
  # setwd("/Users/jing/Downloads")
  lesions <- read.table("/Users/jing/Google Drive/Data/try3.all_lesions.conf_75.txt", header=TRUE, sep="\t", skip=0)
  # setwd(inputFolderWork)
  regions <- lesions[,c(1,5)]
  regions <- regions[1:(nrow(regions)/2),]
  regions <- cbind(regions, kind= unlist(strsplit(as.vector(regions[,1]), " Peak"))[seq(1, nrow(regions)*2, by=2)] )
  regions <- regions[,-1]
  regions <- unique(regions)

  chrom <-  unlist(strsplit(as.vector(regions[,1]), ":"))[seq(1, nrow(regions)*3, by=3)]
  chrom <-   unlist(strsplit(chrom, "chr"))[seq(2, nrow(regions)*2, by=2)]

  start <- unlist(strsplit(as.vector(regions[,1]), ":"))[seq(2, nrow(regions)*3, by=3)]
  start <- as.numeric(unlist(strsplit(start, "-"))[seq(1, nrow(regions)*2, by=2)])

  end <- unlist(strsplit(as.vector(regions[,1]), ":"))[seq(2, nrow(regions)*3, by=3)]
  end <- unlist(strsplit(end, "-"))[seq(2, nrow(regions)*2, by=2)]
  end <- unlist(strsplit(end, "probes"))[seq(1, nrow(regions)*2, by=2)]
  end <- as.numeric(substr(end, 1, nchar(end)-1))

  regions <- cbind(regions, chrom=chrom, start=start, end=end)
  regions[,1] <- paste(regions[,1], regions[,2])
  rownames(regions) <- as.vector(regions[,1])


  listmeans <- list()
  for ( i in 1:nrow(regions)){

    temp  <-copynum[which(copynum$chrom==as.vector(regions$chrom[i])),]
    w1 <- which(temp$start>=regions$end[i])
    if(length(w1)>0) {temp <- temp[-w1,] }
    w2 <- which(temp$end<=regions$start[i])
    if(length(w2)>0) {temp <- temp[-w2,] }

    temp$start[which(temp$start<regions$start[i])] <- regions$start[i]
    temp$end[which(temp$end>regions$end[i])] <- regions$end[i]

    meanvals <- vector()


    for( j in 1:length(samps)){
      tempsub <- temp[which(temp[,1]==samps[j]),]
      tempsub <- cbind(tempsub, length=tempsub$end-tempsub$start)
      meanvals[j] <- sum(tempsub$length* tempsub$value) /   sum(tempsub$length)
    }
    listmeans[[i]] <- meanvals
  }

  alltog <- listmeans[[1]]

  for(i in 2:length(listmeans)){ alltog <- rbind(alltog, listmeans[[i]])}

  colnames(alltog) <- samps
  rownames(alltog) <- rownames(regions)
  ########## remove NA sample
  alltog <- alltog[,-which(colnames(alltog)=="RB64")]
  ############## tetraploid sample
  alltog[,which(colnames(alltog)=="RB57")] <-  alltog[,which(colnames(alltog)=="RB57")]-2

  regions <- cbind(regions, alltog)
  temp=apply(alltog,2,function(z){-(z<2)+(z>2)})

  freqs <- regions[,1:5]
  perc <- vector()
  for( i in 1:nrow(freqs)){
    if(freqs[,2][i]=="Amplification"){ perc[i]=length(which(temp[i,]==1))/ncol(temp)*100}
    if(freqs[,2][i]=="Deletion"){ perc[i]=length(which(temp[i,]==-1))/ncol(temp)*100}
  }
  freqs<- cbind(freqs, perc=perc)

  getwd()
  save(alltog, freqs, regions, file="GNL-gistic.RData")

  gnl<- alltog
  gnl[which(gnl>4)]=4



  ###### consensus clustering copy number

  projectName <- "CNV"

  getwd()
  clustering <- cit.clusteringAnalysis(projectName, data = as.data.frame(gnl), doTests = F,
                                       q.UnsupSelection = c( 0.5, 0.6, 0.7, 0.8),
                                       K.clustering=6 , rowCentering=TRUE,
                                       meth.Clustering= c("ward","complete", "average")
                                       )
  mp <- clustering$multipartition

  # if doTests = T, function is not working.
  #   Error in png(paste(projectTitle, "asso.clust.annot.bestpvals-", ifelse(rowCentering,  :
  #   invalid quartz() device size
  # if doTests = F, the function can be run, but q.UnsupSelection is not working.
  # Therefore, I used the pre-calculated RData instead:
  mp <- get(load("/Users/jing/Google Drive/Data/retino-07-2017-cnv-multPart-centered.RData"))$multipartition

  # rm(mp)
  # mp <- get(load( paste(projectName, "-multPart-centered.RData", sep="") ) )$multipartition
  # mp <- get(load("retino-cnv-multPart-centered.RData"))$multipartition

  mpK2 <- mp[,grep("K=5",names(mp))]
  ch <- cit.getConsHier(multpart =mpK2,K=5)
  conspart <- ch$consensusPartition   # partition consensus
  conshier <- ch$consensusHierarchy    # hierarchie consensus
  coclas <- as.data.frame(cit.countPairs(mpK2))  # matrice de co-classement
  rownames(coclas) <- names(coclas) <- rownames(mpK2)

  # cnv_group<- conspartc
  palette <- cit.myPalette(low = "white", high = "blue", k = 11)
  persoc5 <- colorRampPalette( c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",  "#FB9A99"))
  dev.new(width = 10, height = 10)
  tut=cit.heatmap(as.matrix(coclas), colclust.hclust = conshier, rowclust.hclust = conshier,
                  colclust.colors = "black",
                  heatmapcolors = palette,
                  colpart = conspart,
                  # colpart.bgcol =  persoc5,
                  # colpart.bgcol =  c("#A7CEE4","#1F78B4") ,
                  colclust.heightlabel = "",
                  collabels.cex = 0.001,collabels.ymax = 0.01,
                  colpart.labels = "",
                  ## suppression des lables samples et label de groupes
                  lw = c(2, 5, 1, 20),left.margin = 0, top.margin = 0,lh = c(7, 1,2, 20),title = projectName)
  print(tut)
  dev.print(png,filename = "Figure1a_upperpanel_ConsensusClustering_CNV.png", res=360,width=2400,height=2400)

  # cit.clusteringAnalysis("retino-cnv",data=gnl, K.clustering=6 , meth.Clustering= c("ward","complete", "average"), q.UnsupSelection=c( 0.5,0.6,0.7, 0.8))
  #
  # mp <- get(load("retino-cnv-multPart-centered.RData"))$multipartition
  # # rm(mp)
  # mpK2 <- mp[,grep("K=5",names(mp))]
  # ch <- cit.getConsHier(multpart=mpK2,K=5)
  # # partition consensus:
  # conspartc <- ch$consensusPartition
  # # hiérarchie consensus:
  # conshierc <- ch$consensusHierarchy
  # # matrice de co-classement
  # coclasc <- as.data.frame(cit.countPairs(mpK2))
  # names(coclasc) <- rownames(mpK2)
  # rownames(coclasc) <- rownames(mpK2)
  #
  # cnv_group<- conspartc


  # dev.new()
  # EMA::clustering.plot(as.matrix(coclas), tree =conshier, tree.sup =conshier, heatcol = palette, trim.heatmap=0.96, lab=conspart, legend=FALSE, palette="persoc5", scale="none")
  # dev.print(png,filename = "figure_1_A_coclasc.png", res=360,width=2400,height=2400)


}





#' Plot Figure 1a Middle and Bottom Panel, clusters of clusters and centroid clustering plot # plot only
#' @export
Figure1a_Plot_ClustersofClusters_CentroidClustering <- function(){

  load(file="/Users/jing/Google Drive/Data/cluster_centroid_colors.RData")
  alldatacols_part_C <- cbind(framecol[,c(2,3,4,6)],cor_centroids_col[,c(3,4,5)])


  dev.new()

  plot(0,0, frame=FALSE, xlim=c(0,110), ylim=c(0,20), axes=FALSE, cex=0.001, xlab="", ylab="")
  abs <-15
  poscol <- 1:nrow(alldatacols_part_C)
  hauteurs <- rep(0.7, ncol(alldatacols_part_C))
  gaps <- rep(0, length(hauteurs))
  gaps[3]=0.5
  gaps[4]=5
  gaps[6]=0.5
  for(j in 1:(ncol(alldatacols_part_C))){
    for(i in 1:nrow(alldatacols_part_C)){
      if(alldatacols_part_C[,j][i] == "lightgrey") {
        alldatacols_part_C[,j][i] <- "grey"
      }
      rect(poscol[i],abs, poscol[i]+1,abs-hauteurs[j],  col= alldatacols_part_C[,j][i], border="grey80", lwd=0.5)
    }
    abs <- abs - hauteurs[j] - gaps[j]
  }

  dev.print(png,filename = "Figure1a_middlepanel_ClustersofClusters_bottompanel_CentroidClustering.png", res=360,width=2400,height=2400)
  dev.off()
}




########
#' Plot Figure 1b left Panel, heatmap of 9 CpGs # to be finished
#'
#' @export
Figure1b_Heatmap_9CpG_metharray <- function(){
  require(gplots)

  load("/Users/jing/Google Drive/Data/methylation_data.RData")
  # load("/Volumes/LaCie/Work/Database/Retinoblastoma/DataCancerCell/group102.RData")
  # group.i = group102[which(!is.na(group102$FinalGroup)),]

  group.i = TableS1_FinalSubtype[which(TableS1_FinalSubtype$final_subtype !=3),]

  s = intersect(rownames(group.i), colnames(meth))
  gp = group.i[s, "final_subtype"]
  cpg = colnames(pyro.sig)[2:10]

  dev.new()
  heatmap.2(meth[cpg,s],
            distfun = function(x) d = as.dist(1 - cor(t(x))),
            hclustfun = function(d) hclust(d, method = "complete"),
            Rowv = F,
            scale = "none",
            dendrogram = "none",
            margins = c(5, 10),
            col = colorRampPalette(c("blue", "black", "yellow")),
            ColSideColors = c("goldenrod", "cornflowerblue")[gp],
            trace = "none")
  dev.print(png, "Fig1b_leftpanel_Heatmap_9CpG.png", res = 100, height = 800, width = 800)

  dev.off()


}
Figure1b_CentroidPredictor_9CpG <- function(){

  meth <- get(load("/Users/jing/Google Drive/Data/methylation_data.RData"))
  dim(meth)

  methc <- sweep(meth, 1, rowMeans(meth), "-")

  signsm <- vector()
  pvalsm <- vector()

  methc_group1 <- methc[,intersect(rownames(frame)[which(frame[,5]==1)], colnames(methc))]
  methc_group2 <- methc[,intersect(rownames(frame)[which(frame[,5]==2)], colnames(methc))]
  methc_unclass <- methc[,intersect(rownames(frame)[which(frame[,5]==3)], colnames(methc))]
  if(F){
    for (i in 1:nrow(methc)){
      pvalsm[i] <- wilcox.test( methc_group1[i,], methc_group2[i,])$p.value
      signsm[i] <- sign( mean(methc_group1[i,]) - mean(methc_group2[i,] ))
    }
    save(pvalsm,signsm,file="C:/CIT/projets_Autres/P_Meriem_Sefta/pvalsm_signsm.RData")
  }else{
    load("C:/CIT/projets_Autres/P_Meriem_Sefta/pvalsm_signsm.RData")
  }
  sig_cpgs <- c(rownames(methc)[which(signsm==1)][order(pvalsm[which(signsm==1)])][1:5000],
                rownames(methc)[which(signsm==-1)][order(pvalsm[which(signsm==-1)])][1:5000])

  centroid1 <- rowMeans(methc_group1[sig_cpgs,])
  centroid2 <- rowMeans(methc_group2[sig_cpgs,])

  curie.methy.centroids=data.frame(centroid1,centroid2,row.names=sig_cpgs)                                                # added by AdR
  save(curie.methy.centroids,file="C:/CIT/projets Autres/Meriem Sefta/curie.methy.centroids.RData")                       # added by AdR


  S<-intersect(cit.which(frame,5,1:2),colnames(meth)); z= cit.supervised::my.test(meth[,S],frame[S,5],type.test="wilcox") # added by AdR
  zd=z[which(z[,6]<0),]                                                                                                   # added by AdR
  zd[intersect(cit.which(zd,1,20,"lowest"),cit.which(zd,6,-.4,"<")),c(1,4:6)]                                             # added by AdR
  #           wilcox.test p.values      GM.1      GM.2    FC.1vs2
  #cg12750745         8.859330e-14 0.1285418 0.5541586 -0.4256168
  #cg08091439         1.157760e-13 0.2400174 0.6874420 -0.4474246    # <<---
  #cg23877497         1.505587e-13 0.2331429 0.6819304 -0.4487875
  #cg11324957         1.944947e-13 0.3973842 0.8604477 -0.4630635
  #cg19787076         1.944947e-13 0.2322002 0.7062312 -0.4740310
  zu=z[which(z[,6]>0),]                                                                                                   # added by AdR
  zu[intersect(cit.which(zu,1,20,"lowest"),cit.which(zu,6,.4,">")),c(1,4:6)]                                              # added by AdR
  #           wilcox.test p.values      GM.1      GM.2   FC.1vs2
  #cg20641531         2.218993e-16 0.6682888 0.1552796 0.5130092
  #cg03670369         1.664245e-15 0.8413417 0.2757081 0.5656335
  #cg10007051         2.496367e-15 0.5610808 0.1567623 0.4043185
  #cg07857792         3.716813e-15 0.7318207 0.2670876 0.4647331
  #cg17341366         3.716813e-15 0.6596688 0.1988143 0.4608545
  sig_cpgs <- c( intersect(cit.which(zd,1,20,"lowest"),cit.which(zd,6,-.4,"<")),                                          # added by AdR
                 intersect(cit.which(zu,1,20,"lowest"),cit.which(zu,6,.4,">")) )                                          # added by AdR
  reducedCpGSig=annotCpGs[sig_cpgs,]                                                                                      # added by AdR
  write.table(reducedCpGSig,file="C:/CIT/projets Autres/Meriem Sefta/signature-centroide-methy-reduiteV2.txt",sep="\t")     # added by AdR


  zd=z[which(z[,6]<0),]                                                                                                   # added by AdR
  zd[intersect(cit.which(zd,1,50,"lowest"),cit.which(zd,6,-.38,"<")),c(1,4:6)]                                             # added by AdR
  #           wilcox.test p.values      GM.1      GM.2    FC.1vs2
  #cg12750745         8.859330e-14 0.1285418 0.5541586 -0.4256168
  #cg08091439         1.157760e-13 0.2400174 0.6874420 -0.4474246
  #cg23877497         1.505587e-13 0.2331429 0.6819304 -0.4487875
  #cg11324957         1.944947e-13 0.3973842 0.8604477 -0.4630635
  #cg19787076         1.944947e-13 0.2322002 0.7062312 -0.4740310
  #cg09470010         4.070743e-13 0.1612397 0.5420055 -0.3807658     <- NEW
  #cg18693822         1.963587e-12 0.3243872 0.7083503 -0.3839631     <- NEW

  zu=z[which(z[,6]>0),]                                                                                                   # added by AdR
  zu[intersect(cit.which(zu,1,30,"lowest"),cit.which(zu,6,.4,">")),c(1,4:6)]                                              # added by AdR
  #           wilcox.test p.values      GM.1      GM.2   FC.1vs2
  #cg20641531         2.218993e-16 0.6682888 0.1552796 0.5130092
  #cg03670369         1.664245e-15 0.8413417 0.2757081 0.5656335
  #cg10007051         2.496367e-15 0.5610808 0.1567623 0.4043185
  #cg07857792         3.716813e-15 0.7318207 0.2670876 0.4647331
  #cg17341366         3.716813e-15 0.6596688 0.1988143 0.4608545
  #cg21214455         5.381058e-15 0.6511162 0.1647284 0.4863878       <- NEW
  #cg10316527         5.381058e-15 0.5646209 0.1613257 0.4032952       <- NEW

  sig_cpgs <- c( intersect(cit.which(zd,1,50,"lowest"),cit.which(zd,6,-.38,"<")),                                          # added by AdR
                 intersect(cit.which(zu,1,30,"lowest"),cit.which(zu,6,.4,">")) )                                          # added by AdR
  reducedCpGSig=cbind(z[sig_cpgs,c(1,4:6)],annotCpGs[sig_cpgs,] )                                                                                     # added by AdR
  write.table(reducedCpGSig,file="C:/CIT/projets Autres/Meriem Sefta/signature-centroide-methy-reduiteV3.txt",sep="\t")     # added by AdR

  sig_cpgs <- setdiff(sig_cpgs , c("cg11324957","cg03670369"))

  centroid1 <- rowMeans(methc_group1[sig_cpgs,])
  centroid2 <- rowMeans(methc_group2[sig_cpgs,])


  cor_centroid1 <- vector()
  cor_centroid2 <- vector()

  for(i in 1:ncol(methc)){
    cor_centroid1[i] <-  cor.test(centroid1, methc[sig_cpgs,][,i])$estimate
    cor_centroid2[i] <-  cor.test(centroid2, methc[sig_cpgs,][,i])$estimate
  }


  cor_meth_centroids <- data.frame(names=colnames(methc), cor_meth_centroid1=cor_centroid1, cor_meth_centroid2=cor_centroid2)

  names=rownames(frame)[which(frame[,5]==3)]

  table(apply( cor_meth_centroids [,2:3],1,which.max),frame[as.character(cor_meth_centroids [,1]),5])
  #
  #     1  2  3
  #  1 24  0  5
  #  2  0 36  1


  TableS1_FinalSubtype


  m = as.matrix(pyro.sig[,3:11])
  rownames(m) = pyro.sig$Sample
  md = apply(m, 2, median)
  m.i = t(apply(m, 2, "-", md))
  d.i = as.dist(1- cor(t(m.i), method = "pearson"), )
  heatmap(m.i,
          distfun = function(x) as.dist(1- cor(t(x), method = "pearson"), ),
          hclustfun = function(d, method = "centroid") hclust(d,method),
          scale = "none",
          ColSideColors = ifelse(group102[intersect(rownames(m), rownames(group102)),1] == 1, "goldenrod", "cornflowerblue"))



} # analysis


#' Plot Figure 1b middle Panel, scatter of 9 CpGs correlating CpGarray and pyroseq
#' @export
Figure1b_Scatter_CpGarray_pyroseq <- function(){
  require(tidyr)
  require(dplyr)
  require(RColorBrewer)
  require(colorRamps)

  s = intersect(rownames(array.sig), rownames(pyro.sig))

  # plot(array.sig[,3], pyro.sig[,3])
  colnames(array.sig)
  df.array <- array.sig[s, ] %>% gather(cpg, methlevel, cg12750745:cg10316527)
  df.pyro <- pyro.sig[s, ] %>% gather(cpg, methlevel, cg12750745:cg10316527)

  df <- merge(x=df.array, y=df.pyro, by=c("sample", "group", "cpg"))
  df$cpg <- factor(df$cpg, levels=c("cg12750745", "cg08091439", "cg23877497", "cg20641531", "cg10007051", "cg07857792", "cg17341366", "cg21214455", "cg10316527"))
  colnames(df)[4:5] <- c("array", "pyro")


  dev.new()
  par(mfrow = c(1,1))
  par(mar=c(5,5,3,1))
  par(mgp=c(3,1,0))
  #plot(df$array, df$pyro, col = col.var, main = paste("Correlation of methylome and pyroseq\nr = ", signif(cor(df$array, df$pyro),2), sep = ""),
  #     ylab = "Methylation level by pyroseq", xlab = "Methylation level by Illumina 450K array", pch = 20, cex = 1.5, cex.lab = 2, cex.main = 2, cex.axis = 2,
  #     frame=FALSE, ylim = c(0,1.15), xlim = c(0,1))
  plot(df$array*100, df$pyro*100, col = primary.colors(11)[df$cpg],
       ylab = "Methylation level by pyrosequencing (%)       ", xlab = "Methylation level by Illumina 450K array (%)",
       pch = 20, cex = 1.5, cex.lab = 2, cex.main = 2, cex.axis = 2,
       frame=FALSE, ylim = c(0,115), xlim = c(0,100))
  legend("topleft", pch = 20, legend = levels(df$cpg),col = primary.colors(11), bty = "n", cex = 1.75)
  dev.print(png,filename = "Figure1b_middlepanel_cor_pyroseq_metharray.png", res=360,width=2760,height=3120)

  res_cor.test <- cor.test(df$array*100, df$pyro*100)
  capture.output(res_cor.test, file = "correlation p value.txt")

}


#' Plot Figure 1b right Panel, classification of retinoblastoma by omics/9-CpG predictor
#' @export
Figure1b_Classification_CpGpredictor <- function(){

  ### centroid / pyro Validation Set

  centroid <- c(rep("goldenrod",10),rep("cornflowerblue",8))
  dev.new()
  plot(0,0, frame=FALSE, xlim=c(0,50), ylim=c(0,10), axes=FALSE, cex=0.001, xlab="", ylab="")
  abs <-8
  poscol <- 1:length(centroid)
  hauteurs <- 0.7
  for(i in 1:length(centroid)){
    rect(poscol[i],abs, poscol[i]+1,abs-hauteurs,  col= centroid[i], border="grey80", lwd=0.5)
  }

  dev.print(png,filename = "Figure1b_rightpanel_InitialSeries_omicspredictor.png", res=360,width=2400,height=2400)


  ######## initial set avec NA

  ### centroid / pyro Validation Set

  centroid <- c(rep("goldenrod",11),rep("cornflowerblue",9))
  centroid_2 <- c(rep("goldenrod",10),"grey",rep("cornflowerblue",8),"grey")
  dev.new()
  plot(0,0, frame=FALSE, xlim=c(0,50), ylim=c(0,10), axes=FALSE, cex=0.001, xlab="", ylab="")
  abs <-8
  poscol <- 1:length(centroid)
  hauteurs <- 0.7
  for(i in 1:length(centroid)){
    rect(poscol[i],abs, poscol[i]+1,abs-hauteurs,  col= centroid[i], border="grey80", lwd=0.5)
  }
  abs <- 7
  for(i in 1:length(centroid)){
    rect(poscol[i],abs, poscol[i]+1,abs-hauteurs,  col= centroid_2[i], border="grey80", lwd=0.5)
  }

  dev.print(png,filename = "Figure1b_rightpanel_InitialSeries_9CpGpredictor.png", res=360,width=2400,height=2400)


  ###################################################
  ### Additional Samples

  pyro<- c(rep("goldenrod",7),rep("cornflowerblue",20),rep("grey",3))
  dev.new()
  plot(0,0, frame=FALSE, xlim=c(0,50), ylim=c(0,10), axes=FALSE, cex=0.001, xlab="", ylab="")
  abs <-8
  poscol <- 1:length(pyro)
  hauteurs <- 0.7
  for(i in 1:length(pyro)){
    rect(poscol[i],abs, poscol[i]+1,abs-hauteurs,  col= pyro[i], border="grey80", lwd=0.5)
  }

  dev.print(png,filename = "Figure1b_rightpanel_AdditionalSamples_9CpGpredictor.png", res=360,width=2400,height=2400)


}





########
#' Plot Figure 1c, Heatmap of different clinical features in two subtypes of retinoblastoma
#' @export
Figure1c_Plot_Clinical <- function(input){

  load("/Users/jing/Google Drive/Data/Data_i_figure1_C.RData")

  ### on met tout en couleur ....
  library(EMA)

  alldatacols <- data_i
  nacol <- "white"

  # ## Groupe Cone-like / Mixed-type
  alldatacols[,"groupe_final"] <- as.vector(alldatacols[,"groupe_final"])
  groupcols <- rep(nacol, nrow(alldatacols))
  groupcols[which(alldatacols$groupe_final==1)]="goldenrod"
  groupcols[which(alldatacols$groupe_final==2)]="cornflowerblue"
  groupcols[which(alldatacols$groupe_final=="NA")]="grey"
  alldatacols$groupe_final <- groupcols

  ### RB germline
  alldatacols[,"RB1_germ"] <- as.vector(alldatacols[,"RB1_germ"])
  alldatacols[,"RB1_germ"][which(alldatacols[,"RB1_germ"]=="no")] <- "grey"
  alldatacols[,"RB1_germ"][which(alldatacols[,"RB1_germ"]!="grey")] <-"black"
  alldatacols[,"RB1_germ"][which(is.na(alldatacols[,"RB1_germ"])==TRUE)] <- "white"

  ### AGE
  agecol <- rep(nacol, nrow(alldatacols))
  agecol[which(alldatacols$Age<19)]="navajowhite1"
  agecol[which(alldatacols$Age>=19 & alldatacols$Age/30.5<35)]="orangered"
  agecol[which(alldatacols$Age>=35)]="red4"
  alldatacols$Age_diagnosis<-agecol

  ## Laterality
  latcols <- rep(nacol, nrow(alldatacols))
  latcols[which(alldatacols$Laterality=="unilateral")]="skyblue1"
  latcols[which(alldatacols$Laterality=="bilateral")]="darkblue"
  alldatacols$Laterality <- latcols

  ## Growth pattern
  growthcols <- rep(nacol, nrow(alldatacols))
  growthcols[which(alldatacols$growth=="mixed")]="seagreen"
  growthcols[which(alldatacols$growth=="exophytic")]="purple"
  growthcols[which(alldatacols$growth=="endophytic")]="yellow"
  alldatacols$growth <- growthcols

  ## Diameter tumor
  diamcols <- myPalette(low = "white", high = "navyblue", k = 7)[2:7]
  diamcols <- diamcols[as.numeric(cut(alldatacols$Diameter_mm, breaks=6))]
  diamcols[which(is.na(diamcols))]=nacol
  alldatacols$Diameter_mm <- diamcols

  ### Necrosis
  neccols <- rep(nacol, nrow(alldatacols))
  neccols[which(alldatacols$necrosis=="necrosis_present")]="sienna"
  neccols[which(alldatacols$necrosis=="none")]="gray82"
  alldatacols$necrosis <- neccols

  ### Optic
  optncols<- myPalette(low = "paleturquoise1", high = "cyan4", k = 4)
  alldatacols$optic_nerve_invasion <- as.vector(alldatacols$optic_nerve_invasion)
  alldatacols$optic_nerve_invasion[which(alldatacols$optic_nerve_invasion=="none")]=optncols[1]
  alldatacols$optic_nerve_invasion[which(alldatacols$optic_nerve_invasion=="absent")]=optncols[1]
  alldatacols$optic_nerve_invasion[which(alldatacols$optic_nerve_invasion=="prelaminar")]=optncols[2]
  alldatacols$optic_nerve_invasion[which(alldatacols$optic_nerve_invasion=="intralaminar")]=optncols[3]
  alldatacols$optic_nerve_invasion[which(alldatacols$optic_nerve_invasion=="postlaminar")]=optncols[4]
  alldatacols$optic_nerve_invasion[which(is.na(alldatacols$optic_nerve_invasion))]=nacol

  ##Choroid
  chsccols<- myPalette(low = "lightpink", high = "deeppink3", k = 5)
  alldatacols$choroid_and_sclera_invasion <- as.vector(alldatacols$choroid_and_sclera_invasion)
  alldatacols$choroid_and_sclera_invasion[which(alldatacols$choroid_and_sclera_invasion=="none")]=chsccols[1]
  alldatacols$choroid_and_sclera_invasion[which(alldatacols$choroid_and_sclera_invasion=="choroid_minimal")]=chsccols[2]
  alldatacols$choroid_and_sclera_invasion[which(alldatacols$choroid_and_sclera_invasion=="choroid_deep")]=chsccols[3]
  alldatacols$choroid_and_sclera_invasion[which(alldatacols$choroid_and_sclera_invasion=="choroid_extended")]=chsccols[4]
  alldatacols$choroid_and_sclera_invasion[which(alldatacols$choroid_and_sclera_invasion=="choroid_extended_and_sclera")]=chsccols[5]
  alldatacols$choroid_and_sclera_invasion[which(is.na(alldatacols$choroid_and_sclera_invasion))]=nacol

  rownames(alldatacols) <- alldatacols[,"ID"]

  # save(alldatacols,file="all_data_color_figure1_C_juillet_2017.RData")

  #### plot figure
  ### nouvel ordre des données
  new_order_col <- c("ID","groupe_final","RB1_germ","Laterality","Age_diagnosis","growth","Diameter_mm","necrosis","optic_nerve_invasion","choroid_and_sclera_invasion")

  alldatacols <- alldatacols[,new_order_col]
  # save(alldatacols,file="all_data_color_figure1_C_juillet_2017.RData")
  # poscol <- c(1:39, 41:98, 100:105)
  poscol <- c(1:38, 42:99, 102:107)
  dev.new()

  plot(0,0, frame=FALSE, xlim=c(0,110), ylim=c(0,20), axes=FALSE, cex=0.001, xlab="", ylab="")
  abs <-18
  #poscol <- 1:nrow(alldatacols)

  hauteurs <- rep(0.7, ncol(alldatacols))
  gaps <- rep(0, length(hauteurs))
  gaps[c(5,6)]=0.5
  gaps[2]=0.7
  for(j in 2:(ncol(alldatacols))){
    print(colnames(alldatacols)[j])
    for(i in 1:nrow(alldatacols)){
      rect(poscol[i],abs, poscol[i]+1,abs-hauteurs[j],  col= alldatacols[,j][i], border="grey80", lwd=0.5)
    }
    abs <- abs - hauteurs[j] - gaps[j]
  }

  dev.print(png,filename = "Figure1c_ClinicalPlot.png", res=360,width=2400,height=2400)



}


