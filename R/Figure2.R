# load("/Users/jing/Google Drive/Data/methylation_data.RData")
# dim(meth)
# input <- meth
#
# load("/Users/jing/Google Drive/Data/TableS1_FinalSubtype.rda")
# load("/Users/jing/Google Drive/Data/TableS2_DataFig2D.rda")

# library(RBsubtyping)



#' Plot Figure 2a-d plot of mutations #### to update
#' @export
Figure2abcd <- function(){

  ##########################################
  #### Part A

  ### Loading data GNL, LOH, po (probe_description) and  annot

  GNL.file <- "GNL_Figure2.RData"
  LOH.file <- "LOH_Figure2.RData"
  po.file <- "probe_description.RData"
  annot.file <- "annot.txt"
  resdir <- "Figure_2/Part_A/"

  load(GNL.file);dim(GNL)
  load(LOH.file);dim(LOH)
  load(po.file);dim(po)
  po$pos <- po$Start
  annot <- read.delim(annot.file,as.is=T);dim(annot)


  ### Functions
  geco.annotToCol <- function(annotS=NULL,
                              annotT=NULL, #table d'annotation
                              missing=c("",NA),
                              anotype=NULL,
                              maxnumcateg=2,
                              categCol=NULL,
                              quantitCol=NULL,
                              plotLegend=T,
                              plotLegendFile=NULL
  )
  {
    # Case with only one annotCol
    if(is.null(ncol(annotS))){
      annotS<-data.frame(annotS)
      colnames(annotS)=annotCol
      rownames(annotS)=rownames(annotT)
    }
    # Set all missing values to NA
    for(j in 1:ncol(annotS))	annotS[which(annotS[,j] %in% missing),j] <- NA

    # Get the type of each annotation column
    if(is.null(anotype)){
      anotype <- rep("categ",ncol(annotS));names(anotype) <- colnames(annotS)
      classes <- sapply(1:ncol(annotS),function(j) class(annotS[,j]))
      nmodal <- sapply(1:ncol(annotS),function(j) length(unique(setdiff(annotS[,j],NA))))
      anotype[which(classes %in% c("integer","numeric") & nmodal > maxnumcateg)] <- "quantit"
      anotype[which(nmodal==2)] <- "binary"
    }

    # Convert annotations to colors
    anocol <- annotS

    if(plotLegend)	pdf(plotLegendFile)

    if(is.null(categCol))	categCol <- c("royalblue", "palevioletred1", "red", "palegreen4", "skyblue", "sienna2", "slateblue3", "pink2", "slategray", "black", "orange", "turquoise4", "yellow3", "orangered4", "orchid", "palegreen2", "orchid4", "red4", "peru", "orangered", "palevioletred4", "purple", "sienna4", "turquoise1")

    k <- 1
    for(j in which(anotype=="categ")){
      tmp <- as.factor(anocol[,j])
      classes <- as.character(levels(tmp))
      ncat <- length(levels(tmp))
      if(k+ncat > length(categCol))	categCol <- c(categCol,categCol)
      levels(tmp) <- categCol[k:(k+ncat-1)]
      fill <- as.character(levels(tmp))
      anocol[,j] <- as.character(tmp)
      k <- k+ncat
      if(plotLegend){
        par(mar=c(0,0,0,0))
        plot(-10,axes=F,xlim=c(0,5),ylim=c(0,5),xlab="",ylab="")
        legend(1,5,legend=classes,fill=fill,title=colnames(anocol)[j],xjust=0.5,yjust=1)
      }
    }

    memcol <- c()
    for(j in which(anotype=="binary")){
      new <- setdiff(anocol[,j],c(NA,memcol))
      if(length(new)==2){memcol <- c(memcol,c("dodgerblue4","firebrick"));names(memcol)[(length(memcol)-1):length(memcol)] <- sort(new)}
      if(length(new)==1){memcol <- c(memcol,setdiff(c("dodgerblue4","firebrick"),memcol[setdiff(anocol[,j],c(NA,new))]));names(memcol)[length(memcol)] <- new}
      anocol[,j] <-  as.character(anocol[,j])
      for (z in 1:length(memcol)){
        anocol[which(anocol[,j]==names(memcol)[z]),j] <- memcol[z]
      }

      if(plotLegend){
        par(mar=c(0,0,0,0))
        plot(-10,axes=F,xlim=c(0,5),ylim=c(0,5),xlab="",ylab="")
        classes <- intersect(names(memcol),annotS[,j]);fill <- memcol[classes]
        legend(1,5,legend=classes,fill=fill,title=colnames(anocol)[j],xjust=0.5,yjust=1)
      }
    }

    if(is.null(quantitCol))	quantitCol <- c("orangered1","darkgreen","darkblue","darkred","darkgoldenrod4","darkorchid4","darkolivegreen4","darkorange4","darkslategray")
    k <- 1
    for(j in which(anotype=="quantit")){
      colrange <- colorRampPalette(c("white",quantitCol[k]))(100)
      anocol[,j] <- colrange[round(geco.changeRange(anocol[,j],newmin=1,newmax=100))]
      if(k < length(quantitCol)){k <- k+1}else{k <- 1}
      if(plotLegend){
        par(mar=c(8,2,5,1))
        lims <- seq(-1,1,length.out=200)
        image(matrix(lims,nc=1),col= colrange,axes=F,xlab=colnames(anocol)[j])
      }
    }

    if(plotLegend)	dev.off()

    for(j in 1:ncol(anocol))	anocol[which(is.na(anocol[,j])),j] <- "white"
    as.matrix(anocol)
  }

  geco.imageCol <- function(matcol=NULL,
                            strat=NULL,
                            xlab.cex=0.5,
                            ylab.cex=0.5,
                            drawLines=c("none","h","v","b")[1],
                            ...
  )
  {
    if(is.null(ncol(matcol))){
      matcol<-data.frame(matcol)
      colnames(matcol)=colnames(anocol)
    }
    matcol <- matcol[,ncol(matcol):1]
    if(is.null(ncol(matcol))){
      matcol<-data.frame(matcol)
      colnames(matcol)=colnames(anocol)
    }
    csc <- matcol
    csc.colors <- matrix()
    csc.names <- names(table(csc))
    csc.i <- 1
    for(csc.name in csc.names){
      csc.colors[csc.i] <- csc.name
      csc[csc == csc.name] <- csc.i
      csc.i <- csc.i + 1
    }

    if(dim(csc)[2]==1){
      csc<-matrix(as.numeric(unlist(csc)), nrow = dim(csc)[1])
    }	else {
      csc <- matrix(as.numeric(csc), nrow = dim(csc)[1])
    }

    image(csc, col = as.vector(csc.colors), axes = FALSE, ...)
    if(xlab.cex!=0){
      axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1),colnames(matcol),las = 2,tick = FALSE,cex.axis=xlab.cex, ...)
    }
    if(ylab.cex!=0){
      axis(3, 0:(dim(csc)[1] - 1)/(dim(csc)[1] - 1),rownames(matcol),las = 2,tick = FALSE,cex.axis=ylab.cex, ...)
    }
    if(drawLines %in% c("h","b"))	abline(h=-0.5:(dim(csc)[2] - 1)/(dim(csc)[2] - 1));box()
    if(drawLines %in% c("v","b"))	abline(v=0.5:(dim(csc)[1] - 1)/(dim(csc)[1] - 1));box()
    if(!is.null(strat)){
      z <- factor(matcol[,strat]);levels(z) <- 1:length(levels(z))
      z <- geco.vectorToSegments(as.numeric(z))
      abline(v=geco.changeRange(c(0.5,z$Ind_K+0.5)/max(z$Ind_K),newmin=par()$usr[1],newmax=par()$usr[2]),lwd=2,lty=2)
    }
  }

  geco.changeRange <- function (v, newmin = 1, newmax = 10)
  {
    oldmin <- min(v, na.rm = TRUE)
    oldmax <- max(v, na.rm = TRUE)
    newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
  }




  ### 3. Compute and plot aberration frequencies (panel A)
  indAut <- which(po$Chr %in% 1:22);length(indAut)
  GNL <- GNL[indAut,]
  LOH <- LOH[indAut,]
  po <- po[indAut,]

  all(colnames(GNL) == annot$ID)
  #ind <- which(annot$Platform %in% c("370k","Affy2","WES_Curie","WES_IGfirst"));length(ind) # 87 tumors with LOH data
  nLOH <- sum(annot$Platform %in% c("370k","Affy2","WES_Curie","WES_IGfirst"));nLOH
  fgain <- apply(GNL,1,function(z){mean(z > 0,na.rm=T)})
  floss <- apply(GNL,1,function(z){mean(z < 0,na.rm=T)})
  fampl <- apply(GNL,1,function(z){mean(z==2,na.rm=T)})
  fhdel <- apply(GNL,1,function(z){mean(z==-2,na.rm=T)})
  floh <- apply(LOH,1,mean,na.rm=T)*ncol(LOH)/nLOH+floss # we correct LOH frequency that is only studied in nLOH samples (with SNP array or WES)
  png(file.path(resdir,"All_tumors.png"),width=1400,height=400,res=150)
  par(mar=rep(2,4))
  plot(fgain,type="h",col="indianred1",ylim=c(-1,1),lwd=0.2,las=1,axes=F,xlab="",ylab="")
  points(fampl*2,type="h",col="darkred",ylim=c(-1,1))
  points(-floh,type="h",col="gray75",ylim=c(-1,1))
  points(-floss,type="h",col="skyblue3",ylim=c(-1,1))
  points(-fhdel*2,type="h",col="darkblue",ylim=c(-1,1))
  chrom <- po$Chr
  tmp <- 1;mycol <- rep(c("white","grey20"),20)
  for(C in unique(chrom))
  {
    d <- min(which(chrom==C))
    f <- max(which(chrom==C))
    segments(x0=d,y0=0.9,x1=d,y1=1)
    #rect(xleft=d,ybottom=-0.88,xright=f,ytop=-0.82,col=mycol[match(C,unique(chrom))])
    #if(C==22){segments(x0=f,y0=0.9,x1=f,y1=1)}
    text(x=(d+f)/2,y=tmp,labels=C,cex=0.8);if(tmp==1){tmp <- 0.9}else{tmp <- 1.}
    #segments(f,-0.82,f,0.6,lty=3)
  }
  tmarks <- seq(-0.8,0.8,by=0.4)
  axis(side=2,at=tmarks,labels=abs(tmarks)*100,las=1,pos=-1e4,cex.axis=0.8)
  axis(side=4,at=tmarks,labels=abs(tmarks)*50,las=1,pos=nrow(GNL)+1e4,cex.axis=0.8)
  #mtext("Alteration frequency (%)",side=2,at=0,line=1,cex=0.8)
  dev.off()

  gp1 <- intersect(annot[which(annot$groupe_final==1),"ID"],colnames(GNL));length(gp1) # n=38
  nLOH <- sum(annot$groupe_final==1 & annot$Platform %in% c("370k","Affy2","WES_Curie","WES_IGfirst"),na.rm=T);nLOH # 30 with LOH data
  fgain <- apply(GNL[,gp1],1,function(z){mean(z > 0,na.rm=T)})
  floss <- apply(GNL[,gp1],1,function(z){mean(z < 0,na.rm=T)})
  fampl <- apply(GNL[,gp1],1,function(z){mean(z==2,na.rm=T)})
  fhdel <- apply(GNL[,gp1],1,function(z){mean(z==-2,na.rm=T)})
  floh <- apply(LOH[,gp1],1,mean,na.rm=T)*length(gp1)/nLOH+floss
  png(file.path(resdir,"Groupe1.png"),width=1400,height=400,res=150)
  par(mar=rep(2,4))
  plot(fgain,type="h",col="indianred1",ylim=c(-1,1),lwd=0.2,las=1,axes=F,xlab="",ylab="")
  points(fampl*2,type="h",col="darkred",ylim=c(-1,1))
  points(-floh,type="h",col="gray75",ylim=c(-1,1))
  points(-floss,type="h",col="skyblue3",ylim=c(-1,1))
  points(-fhdel*2,type="h",col="darkblue",ylim=c(-1,1))
  chrom <- po$Chr
  tmp <- 1;mycol <- rep(c("white","grey20"),20)
  for(C in unique(chrom))
  {
    d <- min(which(chrom==C))
    f <- max(which(chrom==C))
    segments(x0=d,y0=0.9,x1=d,y1=1)
    #rect(xleft=d,ybottom=-0.88,xright=f,ytop=-0.82,col=mycol[match(C,unique(chrom))])
    #if(C==22){segments(x0=f,y0=0.9,x1=f,y1=1)}
    text(x=(d+f)/2,y=tmp,labels=C,cex=0.8);if(tmp==1){tmp <- 0.9}else{tmp <- 1.}
    #segments(f,-0.82,f,0.6,lty=3)
  }
  tmarks <- seq(-0.8,0.8,by=0.4)
  axis(side=2,at=tmarks,labels=abs(tmarks)*100,las=1,pos=-1e4,cex.axis=0.8)
  axis(side=4,at=tmarks,labels=abs(tmarks)*50,las=1,pos=nrow(GNL)+1e4,cex.axis=0.8)
  #mtext("Alteration frequency (%)",side=2,at=0,line=1,cex=0.8)
  dev.off()

  gp2 <- intersect(annot[which(annot$groupe_final==2),"ID"],colnames(GNL));length(gp2) # n=58
  nLOH <- sum(annot$groupe_final==2 & annot$Platform %in% c("370k","Affy2","WES_Curie","WES_IGfirst"),na.rm=T);nLOH #51 with LOH data
  fgain <- apply(GNL[,gp2],1,function(z){mean(z > 0,na.rm=T)})
  floss <- apply(GNL[,gp2],1,function(z){mean(z < 0,na.rm=T)})
  fampl <- apply(GNL[,gp2],1,function(z){mean(z==2,na.rm=T)})
  fhdel <- apply(GNL[,gp2],1,function(z){mean(z==-2,na.rm=T)})
  floh <- apply(LOH[,gp2],1,mean,na.rm=T)*length(gp2)/nLOH+floss
  png(file.path(resdir,"Groupe2.png"),width=1400,height=400,res=150)
  par(mar=rep(2,4))
  plot(fgain,type="h",col="indianred1",ylim=c(-1,1),lwd=0.2,las=1,axes=F,xlab="",ylab="")
  points(fampl*2,type="h",col="darkred",ylim=c(-1,1))
  points(-floh,type="h",col="gray75",ylim=c(-1,1))
  points(-floss,type="h",col="skyblue3",ylim=c(-1,1))
  points(-fhdel*2,type="h",col="darkblue",ylim=c(-1,1))
  chrom <- po$Chr
  tmp <- 1;mycol <- rep(c("white","grey20"),20)
  for(C in unique(chrom))
  {
    d <- min(which(chrom==C))
    f <- max(which(chrom==C))
    segments(x0=d,y0=0.9,x1=d,y1=1)
    #rect(xleft=d,ybottom=-0.88,xright=f,ytop=-0.82,col=mycol[match(C,unique(chrom))])
    #if(C==22){segments(x0=f,y0=0.9,x1=f,y1=1)}
    text(x=(d+f)/2,y=tmp,labels=C,cex=0.8);if(tmp==1){tmp <- 0.9}else{tmp <- 1.}
    #segments(f,-0.82,f,0.6,lty=3)
  }
  tmarks <- seq(-0.8,0.8,by=0.4)
  axis(side=2,at=tmarks,labels=abs(tmarks)*100,las=1,pos=-1e4,cex.axis=0.8)
  axis(side=4,at=tmarks,labels=abs(tmarks)*50,las=1,pos=nrow(GNL)+1e4,cex.axis=0.8)
  #mtext("Alteration frequency (%)",side=2,at=0,line=1,cex=0.8)
  dev.off()


  ### 5. Estimating fraction of altered genome and comparison between groups (panel B)
  FAG <- sapply(1:ncol(GNL),function(j)	mean(GNL[,j] !=0,na.rm=T) + mean(LOH[,j] !=0,na.rm=T))
  names(FAG) <- colnames(GNL)
  annot$FAG <- FAG[annot$ID]

  wilcox.test(annot$FAG ~ annot$groupe_final) #p-value = 3.294e-07

  ### 6. MYCN, 1q and 16q in each group
  indMYCN <- which(po$Chr==2 & po$pos > 16.080e6 & po$pos < 16.088e6);length(indMYCN)
  ind1q <- which(po$Chr==1 & po$pos > 126e6);length(ind1q)
  ind16q <- which(po$Chr==16 & po$pos > 37e6);length(ind16q)

  GNL[indMYCN,]
  annot$MYCN <- as.numeric(GNL[indMYCN,][annot$ID]==2)
  tt <- table(annot$MYCN,annot$groupe_final);tt;fisher.test(tt) #p-value = 0.005506
  # annot[which(annot$MYCN==1 & annot$groupe_final==1),]

  tmp <- apply(GNL[ind1q,],2,function(z)	mean(z > 0,na.rm=T));sort(tmp)
  annot$gain1q <- as.numeric(tmp[annot$ID] > 0.4)
  tt <- table(annot$gain1q,annot$groupe_final);tt;fisher.test(tt) #p-value = 5.502e-11

  tmp <- apply(GNL[ind16q,],2,function(z)	mean(z < 0,na.rm=T));sort(tmp)
  annot$loss16q <- as.numeric(tmp[annot$ID] > 0.4)
  tt <- table(annot$loss16q,annot$groupe_final);tt;fisher.test(tt) #p-value = 1.8e-07


  ind2p <- which(po$Chr==2 & po$pos < 95391419)
  GNL[ind2p,]
  tmp <- apply(GNL[ind2p,],2,function(z)	mean(z > 0,na.rm=T));sort(tmp)
  annot$Gain2p <- as.numeric(tmp[annot$ID] > 0.4)
  tt <- table(annot$Gain2p,annot$groupe_final);tt;fisher.test(tt)


  anocol <- geco.annotToCol(annotS=annot[,c("groupe_final","MYCN","gain1q","loss16q")],annotT=annot,plotLegendFile=NULL)
  anocol[,1] <- c("goldenrod","cornflowerblue","grey80")[match(annot[,"groupe_final"],c(1,2,NA))]
  anocol[,2] <- c("grey80","darkred","white")[match(annot[,"MYCN"],c(0,1,NA))]
  anocol[,3] <- c("grey80","indianred1","white")[match(annot[,"gain1q"],c(0,1,NA))]
  anocol[,4] <- c("grey80","skyblue3","white")[match(annot[,"loss16q"],c(0,1,NA))]

  rownames(anocol) <- annot[,"ID"]
  save(anocol,file="annot_col_Fig2_partA.RData")
  save(annot,file="Data_Fig2_partA.Rdata")

  ##############################
  #### PartB

  load(file ="annot_fig.RData")

  groupe_final_MYCN <- annot$groupe_final
  names(groupe_final_MYCN)<- annot$ID

  samples_MYCN_Amp <- c("RB222","RB22","RB14","RB224","RB215","RBsjd2","RBsjd7","RBsjd3","RB659","RB15","RB13")
  annot$groupe_final_MYCN <- groupe_final_MYCN
  annot[samples_MYCN_Amp,"groupe_final_MYCN"] <- "2_MYCN_Amp"

  ### Statistiques ....
  library(ggpubr)

  ### 2 gp
  compare_means(FAG ~ groupe_final, data = annot)
  # # A tibble: 1 x 8
  # .y.   group1 group2           p       p.adj p.format p.signif method
  # <chr> <chr>  <chr>        <dbl>       <dbl> <chr>    <chr>    <chr>
  #   1 FAG   1      2      0.000000329 0.000000329 3.3e-07  ****     Wilcoxon
  # ### 3 gp
  compare_means(FAG ~ groupe_final_MYCN, data = annot)
  # A tibble: 3 x 8
  # .y.   group1     group2                p       p.adj p.format p.signif method
  # <chr> <chr>      <chr>             <dbl>       <dbl> <chr>    <chr>    <chr>
  #   1 FAG   1          2_MYCN_Amp 0.149        0.15        0.149    ns    Wilcoxon
  # 2 FAG   1          2          0.0000000231 0.000000069 2.3e-08  ****     Wilcoxon
  # 3 FAG   2_MYCN_Amp 2          0.0123       0.025       0.012    *        Wilcoxon
  #

  write.table(annot[c("FAG","groupe_final_MYCN")],file="instability_genomic_groupes.csv",sep = ";", col.names = NA)

  dev.new()
  par(mfrow=c(2,2))

  boxplot(1:5, xlim = c(0, 6), ylim=c(0,0.9),
          horizontal=TRUE,outline=F,axes=F) #invisible boxes

  boxplot(annot$FAG ~ annot$groupe_final,frame=FALSE, pch=20, col=c("goldenrod","cornflowerblue"),
          horizontal=TRUE, ylim=c(0,0.45),cex=1.5,names=c("C1","C2"),outline=F,axes=F,at=c(5,4), add = TRUE)

  boxplot(annot[which(annot$groupe_final_MYCN %in% c("2","2_MYCN_Amp")),"FAG"] ~ annot[which(annot$groupe_final_MYCN  %in% c("2","2_MYCN_Amp")),"groupe_final_MYCN"],
          frame=FALSE, pch=20, col=c("cornflowerblue","darkred"),
          horizontal=TRUE, ylim=c(0,0.45),cex=1.5,outline=F,axes=T,at=c(2,1), add = TRUE)

  stripchart(as.numeric(as.vector(annot$FAG)) ~ as.vector(annot$groupe_final), data = annot,
             vertical = FALSE, method = "jitter", jitter = 0.3,
             pch = 20, col = "black", cex= 0.5,at=c(5,4),
             add = TRUE)

  stripchart(annot[which(annot$groupe_final_MYCN  %in% c("2","2_MYCN_Amp")),"FAG"] ~ annot[which(annot$groupe_final_MYCN  %in% c("2","2_MYCN_Amp")),"groupe_final_MYCN"], data = annot,
             vertical = FALSE, method = "jitter", jitter = 0.3,
             pch = 20, col = "black", cex= 0.5,at=c(2,1),
             add = TRUE)

  dev.print(png,filename = "boxplot_instability_genomic_3gp.png",res=360,width=2400,height=2400)


  ##############################
  #### PartC
  ############# NUMBER OF MUTATIONS #######################
  load(file ="annot_fig2_Part_C.RData")

  ### Statistiques ....
  # library(ggpubr)
  # ### 2 gp
  # compare_means(nb_mutations ~ groupe_final, data = number_mutations_groupe)
  # # A tibble: 1 x 8
  # .y.          group1 group2           p       p.adj p.format p.signif method
  # <chr>        <chr>  <chr>        <dbl>       <dbl> <chr>    <chr>    <chr>
  #   1 nb_mutations 2      1      0.000000810 0.000000810 8.1e-07  ****     Wilcoxon
  # ### 3 gp
  # compare_means(nb_mutations ~ groupe_final_MYCN, data = number_mutations_groupe)
  # # A tibble: 3 x 8
  # .y.          group1 group2          p     p.adj p.format p.signif method
  # <chr>        <chr>  <chr>       <dbl>     <dbl> <chr>    <chr>    <chr>
  #   1 nb_mutations 2      1      0.00000354 0.0000106 3.5e-06  ****     Wilcoxon
  # 2 nb_mutations 2      3      0.775      0.775     0.775    ns       Wilcoxon
  # 3 nb_mutations 1      3      0.00105    0.00210   0.001    **       Wilcoxon


  nb_c1 <- length(number_mutations_groupe[which(number_mutations_groupe$groupe_final == "1"),"nb_mutations"] )
  nb_c2 <- length(number_mutations_groupe[which(number_mutations_groupe$groupe_final == "2"),"nb_mutations"] )
  nb_c2_b <- length(number_mutations_groupe[which(number_mutations_groupe$groupe_final_MYCN == "2"),"nb_mutations"] )
  nb_c3 <- length(number_mutations_groupe[which(number_mutations_groupe$groupe_final_MYCN == "3"),"nb_mutations"] )


  number_mutations_groupe[which(number_mutations_groupe$groupe_final_MYCN %in% c("2","3")),]
  write.table(number_mutations_groupe,file="number_mutations_groupes.csv",sep = ";",col.names = NA)

  dev.new()
  par(mfrow=c(2,2))

  boxplot(x=1:5,xlim=c(0,6),ylim=c(0,8),
          horizontal=TRUE,outline=F,axes=F,boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1)) #invisible boxes

  boxplot(number_mutations_groupe[,"nb_mutations"] ~ as.numeric(number_mutations_groupe[,"groupe_final"]),
          col=c("goldenrod","cornflowerblue"),whiskcol="black",staplecol="black",frame=FALSE,
          horizontal=TRUE,cex=0.7,names=c("C1","C2"),col.axis="black",las=2,outline=F,axes=T,at=c(5,4), add=TRUE)

  stripchart(number_mutations_groupe[,"nb_mutations"] ~ as.numeric(number_mutations_groupe[,"groupe_final"]),
             vertical = FALSE, method = "jitter", jitter = 0.5,
             pch = 20, col="black", at=c(5,4),cex= 0.7,
             add = TRUE)

  boxplot(number_mutations_groupe[which(number_mutations_groupe$groupe_final_MYCN %in% c("2","3")),"nb_mutations"] ~ number_mutations_groupe[which(number_mutations_groupe$groupe_final_MYCN  %in% c("2","3")),"groupe_final_MYCN"],
          col=c("cornflowerblue","darkred"),whiskcol=c("black","black"),staplecol=c("black","black"),boxcol=c("black","black"),outcol=c("black","black"),outbg=c("black","black"),
          frame=FALSE,xaxt="n",
          horizontal=TRUE,cex=1,names=c("C2_no_M","C2_MYCN_Amp"),outline=F,axes=T,at=c(2,1),las=2, add=TRUE)

  stripchart(number_mutations_groupe[which(number_mutations_groupe$groupe_final_MYCN %in% c("2","3")),"nb_mutations"] ~ number_mutations_groupe[which(number_mutations_groupe$groupe_final_MYCN  %in% c("2","3")),"groupe_final_MYCN"],
             vertical = FALSE, method = "jitter", jitter = 0.5,
             pch = 20, col="black",  at=c(2,1),cex= 0.7,
             add = TRUE)

  setwd(inputFolderWork)
  dev.print(png,filename = "fig_2_C_3gp.png",res=360,width=2400,height=2400)


  ###############################
  #### Part D
  ###############################
  dir_output <- "./Downloads"
  #### data
  DataFig2d <- read.table("/Users/jing/Downloads/TableS2_DataFig2d.csv", sep = ";", header = T)

  str(DataFig2d)

  #### types of mutations for hit 1:
  unique(DataFig2d$RB1.hit.1)
  # [1] "nonsense"             "missense"             "frameshift indel"     "splice site"          "deletion exons 7-17"
  # [6] "deletion 1 allele"    NA                     "in_frame_indel"       "promoter methylation" "deletion exons 3-27"
  # [11] "deletion exons 7-12"  "no"

  # deletion 1 allele/some exons --> deletion
  DataFig2d[grep("deletion", DataFig2d$RB1.hit.1), ]$RB1.hit.1 <- "deletion"
  unique(DataFig2d$RB1.hit.1)

  # in_frame_indel --> in-frame indel
  DataFig2d[grep("in_frame_indel", DataFig2d$RB1.hit.1), ]$RB1.hit.1 <- "in-frame indel"
  unique(DataFig2d$RB1.hit.1)

  # NA --> "NA"
  DataFig2d[is.na(DataFig2d$RB1.hit.1), ]$RB1.hit.1 <- "NA"
  unique(DataFig2d$RB1.hit.1)

  #### hit 1 germline:
  # NA --> "NA"
  DataFig2d[is.na(DataFig2d$RB1.hit.1.germline), ]$RB1.hit.1.germline <- "NA"
  unique(DataFig2d$RB1.hit.1.germline)

  #### types of mutations for hit 2:
  unique(DataFig2d$RB1.hit.2)
  # [1] "copy neutral LOH"     "nonsense"             "frameshift indel"     NA                     "promoter methylation"
  # [6] "deletion 1 allele"    "splice site"          "deletion exons 18-27" "no"

  # deletion 1 allele/some exons --> deletion
  DataFig2d[grep("deletion", DataFig2d$RB1.hit.2), ]$RB1.hit.2 <- "deletion"
  unique(DataFig2d$RB1.hit.2)

  # NA --> "NA"
  DataFig2d[is.na(DataFig2d$RB1.hit.2), ]$RB1.hit.2 <- "NA"
  unique(DataFig2d$RB1.hit.2)

  #### BCOR.mutation:
  unique(DataFig2d$BCOR.mutation)

  # frameshift_indel --> frameshift indel
  DataFig2d[which(DataFig2d$BCOR.mutation == "frameshift_indel"), ]$BCOR.mutation <- "frameshift indel"

  # NA --> "NA"
  DataFig2d[is.na(DataFig2d$BCOR.mutation), ]$BCOR.mutation <- "NA"
  unique(DataFig2d$BCOR.mutation)

  #### ARID1A.mutation:
  unique(DataFig2d$ARID1A.mutation)

  # frameshift_indel --> frameshift indel
  DataFig2d[which(DataFig2d$ARID1A.mutation == "frameshift_indel"), ]$ARID1A.mutation <- "frameshift indel"

  # NA --> "NA"
  DataFig2d[is.na(DataFig2d$ARID1A.mutation), ]$ARID1A.mutation <- "NA"
  unique(DataFig2d$ARID1A.mutation)



  #### Use ComplexeHatmap to plot
  library(ComplexHeatmap)
  library(circlize)

  str(DataFig2d)

  m <- DataFig2d[, c("MYCN.amplification", "X1q.gain", "X16q.loss")]
  m$MYCN.amplification <- ifelse(m$MYCN.amplification == "yes", 2, 0)
  m$X1q.gain <- ifelse(m$X1q.gain == "yes", 1, 0)
  m$X16q.loss <- ifelse(m$X16q.loss == "yes", -1, 0)

  df.annot <-  DataFig2d[, c("subtype", "RB1.hit.1", "RB1.hit.1.germline", "RB1.hit.2", "BCOR.mutation", "ARID1A.mutation")]
  df.annot$subtype <- paste("Subtype", df.annot$subtype, sep = " ")
  ha <- HeatmapAnnotation(# df = df.annot,
    Subtype = df.annot$subtype,
    Mutation = as.matrix(df.annot[, c("RB1.hit.1", "RB1.hit.1.germline", "RB1.hit.2", "BCOR.mutation", "ARID1A.mutation")]),

    col = list(Subtype = c("Subtype 1" = "goldenrod", "Subtype 2" = "cornflowerblue"),

               Mutation = c("yes" = "black", "no" = "grey82", "NA" = "white",
                            "nonsense" = "salmon", "missense" = "green", "frameshift indel" = "darkslategray1",
                            "splice site" = "yellow", "deletion" = "blue", "in-frame indel" = "plum",
                            "promoter methylation" = "maroon4", "copy neutral LOH" = "purple" )
    ),

    gp = gpar(col = "grey80"),
    gap = unit(c(1,0,0,0,0,0), "mm"),
    simple_anno_size = unit(5, "mm"),
    border = TRUE
  )

  col_fun = colorRamp2(c(-2, -1, 0, 1, 2), c("darkblue", "skyblue3", "grey", "indianred1", "darkred"))
  h <- Heatmap(t(m), cluster_rows = F, cluster_columns = F,
               column_split = df.annot$subtype, column_gap = unit(2, "mm"),
               height = unit(15, "mm"), width = unit(18, "cm"),
               col = col_fun,
               rect_gp = gpar(col = "grey80"),
               border = T,
               top_annotation = ha
  )
  print(h)
  setwd(dir_output)
  dev.print(png, "Fig2d_ComplexeHeatmap.png", res = 300, height = 1500, width = 4000)






}


#' Plot Figure 2efg heatmap of differential methylated CpGs
#' @export
Figure2efg_Meth_FigureS2_Meth <- function(){
  require(gplots)
  require(ComplexHeatmap)
  require(minfi)
  require(GLAD)
  require(ggpubr)
  require(gridExtra)
  require(dplyr)
  require(tidyr)

  #### import data:
  load("/Users/jing/Google Drive/Data/annotCpGs_ilmn12.hg19_minfi1.28.4_Manifest0.4.0.RData")
  load("/Users/jing/Google Drive/Data/methylation_data.RData")

  input = meth


  ## diff gene expression
  df_diffgene <- read.table("/Users/jing/Google Drive/Data/top_table_1.txt", header = T, row.names = 1, sep = "\t", dec=".", quote="")
  df_diffgene.sig <- df_diffgene[which(df_diffgene$adj.P.Val < 0.05), ]
  dim(df_diffgene.sig)


  # prepare samples of interest
  samples_RBc1 <- intersect(rownames(TableS1_FinalSubtype[which(TableS1_FinalSubtype$final_subtype == 1), ] ), colnames(input))
  samples_RBc2 <- intersect(rownames(TableS1_FinalSubtype[which(TableS1_FinalSubtype$final_subtype == 2), ] ), colnames(input))
  samples_RBc2mycn <- intersect(samples_RBc2, rownames(TableS2_DataFig2D)[which(TableS2_DataFig2D$MYCN_amplification == "yes")] )
  samples_RBc2nomycn <- setdiff(samples_RBc2, samples_RBc2mycn)

  #### BEGIN: Prepare statinfo for use later ####
  #### Prepare Variance info
  ## across RBc1 and RBc2
  meth.i <- input[ , c(samples_RBc1, samples_RBc2)]
  dim(meth.i)
  var_meth_RBc1c2 <- apply(meth.i, 1, var) # %>% sort(decreasing=TRUE)

  var_meth_RBc1c2[1:10]

  #### Prepare Diff meth info
  ## C2 vs C1
  meth.i1 <- meth[, samples_RBc1]
  meth.i2 <- meth[, samples_RBc2]
  dim(meth.i1)
  dim(meth.i2)
  p_wilcox_meth_RBc1c2 <- lapply(1:nrow(meth.i1), function(zz) wilcox.test(meth.i1[zz, ], meth.i2[zz, ])$p.value )
  names(p_wilcox_meth_RBc1c2) <- rownames(meth.i1)
  p_wilcox_meth_RBc1c2 <- unlist(p_wilcox_meth_RBc1c2)
  p_wilcox_meth_RBc1c2[1:10]

  p.adj_wilcox_meth_RBc1c2 <- p.adjust(p_wilcox_meth_RBc1c2, method = "BH")
  p.adj_wilcox_meth_RBc1c2[1:10]

  avg_meth_RBc1 <- rowMeans(meth.i1)
  avg_meth_RBc2 <- rowMeans(meth.i2)

  diff_meth_RBc1c2 <- avg_meth_RBc1 - avg_meth_RBc2
  diff_meth_RBc1c2[1:10]

  table(p_wilcox_meth_RBc1c2 < 0.05)
  table(p.adj_wilcox_meth_RBc1c2 < 0.05)
  table(p_wilcox_meth_RBc1c2 < 0.05, diff_meth_RBc1c2 > 0.1)
  table(p.adj_wilcox_meth_RBc1c2 < 0.05, diff_meth_RBc1c2 > 0.1)

  statinfo_meth_df <- cbind.data.frame(var_meth_RBc1c2,
                                       p_wilcox_meth_RBc1c2, p.adj_wilcox_meth_RBc1c2, diff_meth_RBc1c2,
                                       avg_meth_RBc1, avg_meth_RBc2)
  setwd("/Users/jing/Downloads")
  save(statinfo_meth_df, file="statinfo_meth_df_BH.RData")

  statinfo_meth_list <- list(var_meth_RBc1c2 = var_meth_RBc1c2,
                             p_wilcox_meth_RBc1c2 = p_wilcox_meth_RBc1c2, p.adj_wilcox_meth_RBc1c2 = p.adj_wilcox_meth_RBc1c2,
                             diff_meth_RBc1c2 = diff_meth_RBc1c2, avg_meth_RBc1 = avg_meth_RBc1, avg_meth_RBc2 = avg_meth_RBc2)
  names(statinfo_meth_list)
  save(statinfo_meth_list, file="statinfo_meth_list_BH.RData")

  #### END: Prepare statinfo for use later ####



  #### Sample of interest
  samples.i <- c(samples_RBc1, samples_RBc2nomycn, samples_RBc2mycn)
  type_samples.i <- c(
    rep("RBC1", length(samples_RBc1)),
    rep("RBC2woMYCN", length(samples_RBc2nomycn)), rep("RBC2MYCN", length(samples_RBc2mycn)) )
  type_samples.i <- factor(type_samples.i, levels = c("RBC1", "RBC2woMYCN", "RBC2MYCN"))
  color_samples.i <- c("goldenrod", "cornflowerblue", "darkred")[type_samples.i]


  #### cpg of interest to plot
  cpgs.i <- names(p.adj_wilcox_meth_RBc1c2)[which(p.adj_wilcox_meth_RBc1c2 < 0.05 & abs(diff_meth_RBc1c2) > 0.2)]
  length(cpgs.i)


  #### plot heatmap
  dev.new()
  heatmap.2(meth.i[cpgs.i, samples.i],
            Colv = F,
            distfun = function(x) d = as.dist(1 - cor(t(x))),
            hclustfun = function(d) hclust(d, method = "complete"),
            col = colorRampPalette(c("blue", "black", "yellow")),
            ColSideColors = color_samples.i,
            trace = "none")
  dev.print(png, file="heatmap_2020.png", res=300, height=4000, width=3000)
  dev.off()



  #### Select CpGs ####
  annotCpGs.i <- annotCpGs[rownames(input), ]

  levels(factor(annotCpGs$Relation_to_Island))
  # [1] "Island"  "N_Shelf" "N_Shore" "OpenSea" "S_Shelf" "S_Shore"

  #### By Islands
  ## Island
  CpGs_island <- rownames(annotCpGs)[grepl("Island", annotCpGs$Relation_to_Island)]
  length(CpGs_island)
  ## CpGs outside island
  CpGs_outsideisland <- rownames(annotCpGs)[!grepl("Island", annotCpGs$Relation_to_Island)]
  length(CpGs_outsideisland)
  #### CpGs shore
  CpGs_shore <- rownames(annotCpGs)[grepl("Shore", annotCpGs$Relation_to_Island)]
  length(CpGs_shore)
  CpGs_Nshore <- rownames(annotCpGs)[grepl("N_Shore", annotCpGs$Relation_to_Island)]
  length(CpGs_Nshore)
  CpGs_Sshore <- rownames(annotCpGs)[grepl("S_Shore", annotCpGs$Relation_to_Island)]
  length(CpGs_Sshore)
  #### CpGs shore
  CpGs_shelf <- rownames(annotCpGs)[grepl("Shelf", annotCpGs$Relation_to_Island)]
  length(CpGs_shelf)
  CpGs_Nshelf <- rownames(annotCpGs)[grepl("N_Shelf", annotCpGs$Relation_to_Island)]
  length(CpGs_Nshelf)
  CpGs_Sshelf <- rownames(annotCpGs)[grepl("S_Shelf", annotCpGs$Relation_to_Island)]
  length(CpGs_Sshelf)
  #### CpGs Sea
  CpGs_sea <- rownames(annotCpGs)[grepl("Sea", annotCpGs$Relation_to_Island)]
  length(CpGs_sea)
  #### CpGs All
  CpGs_all <- c(CpGs_island, CpGs_outsideisland)


  #### piechart
  plot_piechart = function(zz, input_table, input_table_name) {
    if (grepl(pattern = "island", x = input_table_name) ){
      color.i <- c("forestgreen", "khaki1", "khaki3", "tan1", "tan3", "dodgerblue")
    }
    if (grepl(pattern = "gene", x = input_table_name) ){
      color.i <- c("skyblue4", "skyblue3", "skyblue1", "orangered", "yellow1", "lightgreen")
    }
    df = data.frame(percentage = input_table[,zz], relation_island = rownames(input_table))
    df$relation_island = factor(df$relation_island, levels = rownames(df))
    ## Barplot
    bp <- ggplot(df, aes(x = "", y = percentage, fill = relation_island)) +
      ggtitle( colnames(input_table)[zz] ) + xlab("") +
      geom_bar(width = 1, stat = "identity") +
      scale_fill_manual(values =  color.i ) +
      coord_flip() +
      theme_minimal()
    pc <- bp + coord_polar("y", start=0) +
      theme_minimal()
    return(pc)
  }


  pcs <- lapply(1:ncol(prop_gene), plot_piechart, prop_gene, "prop_gene")


  cpgs.i.up <- names(p.adj_wilcox_meth_RBc1c2)[which(p.adj_wilcox_meth_RBc1c2 < 0.05 & diff_meth_RBc1c2 > 0.2)]
  cpgs.i.dn <- names(p.adj_wilcox_meth_RBc1c2)[which(p.adj_wilcox_meth_RBc1c2 < 0.05 & diff_meth_RBc1c2 < -0.2)]

  annotIsland.i.up <- annotCpGs.i[cpgs.i.up, "Relation_to_Island"]
  length(annotIsland.i.up)
  annotIsland.i.dn <- annotCpGs.i[cpgs.i.dn, "Relation_to_Island"]
  length(annotIsland.i.dn)

  table(annotIsland.i.up)
  table(annotIsland.i.dn)

  prop.table(table(annotIsland.i.up))
  prop.table(table(annotIsland.i.dn))
  prop.table(table(annotCpGs[rownames(beta), ]$Relation_to_Island))

  prop_gene <- cbind(dninC1 = prop.table(table(annotIsland.i.dn)),
                     upinC1 = prop.table(table(annotIsland.i.up)),
                     background = prop.table(table(annotCpGs[rownames(beta), ]$Relation_to_Island))
  )
  prop_gene <- prop_gene[c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"), ]
  prop_gene
  write.table(prop_gene, file="prop_gene.csv", sep=";", col.names = NA)

  table_gene <- cbind(dninC1 = table(annotIsland.i.dn),
                      upinC1 = table(annotIsland.i.up),
                      background = table(annotCpGs[rownames(input), ]$Relation_to_Island)
  )
  table_gene <- table_gene[c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"), ]
  table_gene
  write.table(table_gene, file="table_gene.csv", sep=";", col.names = NA)


  #
  color.i <- c("forestgreen", "khaki1", "khaki3", "tan1", "tan3", "dodgerblue")

  ## up
  df = data.frame(prop.table(table(annotIsland.i.up)))
  colnames(df) = c("relation_island", "percentage")
  df$relation_island = factor(df$relation_island, levels = c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"))

  ## Barplot
  bp <- ggplot(df, aes(x = "", y = percentage, fill = relation_island)) +
    ggtitle( "Hypo-methylated CpGs in Subtype2" ) + xlab("") +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values =  color.i ) +
    coord_flip() +
    theme_minimal()
  pc <- bp + coord_polar("y", start=0) +
    theme_minimal()
  pc

  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )

  pc + blank_theme
  dev.print(png, file="Fig2f_piechart_hypoinC2.png", res=300, height=1200, width=1200)


  ## dn
  df = data.frame(prop.table(table(annotIsland.i.dn)))
  colnames(df) = c("relation_island", "percentage")
  df$relation_island = factor(df$relation_island, levels = c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"))

  ## plot
  bp <- ggplot(df, aes(x = "", y = percentage, fill = relation_island)) +
    ggtitle( "Hyper-methylated CpGs in Subtype2" ) + xlab("") +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values =  color.i ) +
    coord_flip() +
    theme_minimal()
  pc <- bp + coord_polar("y", start=0) +
    theme_minimal()
  pc

  pc+blank_theme
  dev.print(png, file="Fig2f_piechart_hyperinC2.png", res=300, height=1200, width=1200)

  ## background
  df = data.frame(prop.table(table(annotCpGs[rownames(beta), ]$Relation_to_Island)))
  colnames(df) = c("relation_island", "percentage")
  df$relation_island = factor(df$relation_island, levels = c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"))

  ## plot
  bp <- ggplot(df, aes(x = "", y = percentage, fill = relation_island)) +
    ggtitle( "background" ) + xlab("") +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values =  color.i ) +
    coord_flip() +
    theme_minimal()
  pc <- bp + coord_polar("y", start=0) +
    theme_minimal()
  pc
  pc+blank_theme
  dev.print(png, file="Fig2f_piechart_background.png", res=300, height=1200, width=1200)



  #### distribution
  df_diffgene.sig <- df_diffgene[which(df_diffgene$adj.P.Val < 0.05), 4:9]

  cutoff_diff <- 0.2
  relatedtogene <- F
  fun_stat_island <- function(cutoff_diff = 0, relatedtogene = F){
    annot.hyperC1 <- annotCpGs.i[which(p.adj_wilcox_meth_RBc1c2 < 0.05 & diff_meth_RBc1c2 > cutoff_diff ), ]
    if (relatedtogene) {annot.hyperC1 <- annot.hyperC1[which(annot.hyperC1$UCSC_RefGene_Name != ""), ]}
    count.hyperC1 <- table(annot.hyperC1$Relation_to_Island)
    prct.hyperC1  <- round((count.hyperC1 / sum(count.hyperC1)) * 100, 1)

    annot.hyperC2 <- annotCpGs.i[which(p.adj_wilcox_meth_RBc1c2 < 0.05 & diff_meth_RBc1c2 < -cutoff_diff ), ]
    if (relatedtogene) {annot.hyperC2 <- annot.hyperC2[which(annot.hyperC2$UCSC_RefGene_Name != ""), ]}
    count.hyperC2 <- table(annot.hyperC2$Relation_to_Island)
    prct.hyperC2 <- round((count.hyperC2 / sum(count.hyperC2)) * 100, 1)

    stat_island_6cat <- cbind(count.hyperC1, prct.hyperC1, count.hyperC2, prct.hyperC2)
    stat_island_4cat <- rbind(Island = stat_island_6cat["Island", ],
                              Shore = colSums(stat_island_6cat[grep("Shore", rownames(stat_island_6cat)), ]),
                              Shelf = colSums(stat_island_6cat[grep("Shelf", rownames(stat_island_6cat)), ]),
                              Sea = stat_island_6cat["OpenSea", ]  )
    stat_island_2cat <- rbind(Island = stat_island_6cat["Island", ],
                              OutsideIsland = colSums(stat_island_6cat[!grepl("Island", rownames(stat_island_6cat)), ])  )
    chisq.test(stat_island_2cat[, c(1, 3)])$p.value
    return(list(stat = stat_island_2cat, chisq.p = chisq.test(stat_island_2cat[, c(1, 3)])$p.value))
  }


  setwd("/Users/jing/Downloads/")
  fun_plot_meth <- function(CpGs.i, name_CpGs.i, nbyvar = NA, cutoff_p.adj = 0.2, cutoff_diff = 0.15, cutoff_hypo = NA, cutoff_hyper = NA,
                            samples.i = samples_all, type_samples.i = type_samples_all,
                            relatedtogene = FALSE,
                            plot = TRUE) {
    var_meth_RBc1c2.i <- var_meth_RBc1c2[CpGs.i] %>% sort(decreasing = T)
    # var_meth_RBc1c2.i[1:10]
    if (is.numeric(nbyvar)){
      CpGs.input <- names(var_meth_RBc1c2.i[1 : nbyvar])
      title_param <- paste(name_CpGs.i, ", byvar", ", n = ", length(CpGs.input), sep="")
      filename_param <- gsub("=| ", "", gsub(", ", "_", title_param))
    } else if (is.character(nbyvar)) {
      CpGs.input <- names(var_meth_RBc1c2.i[1 : floor(length(var_meth_RBc1c2.i) * as.numeric(gsub("%", "", nbyvar)) / 100)])
      title_param <- paste(name_CpGs.i, ", byvar ", nbyvar, ", n = ", length(CpGs.input), sep="")
      filename_param <- gsub("%", "prct", gsub("=| ", "", gsub(", ", "_", title_param)))
    } else if (is.na(nbyvar)) {
      if (is.na(cutoff_hypo) & is.na(cutoff_hyper)) {
        CpGs.i_hyperinC1 <- intersect(names(which(p.adj_wilcox_meth_RBc1c2 < cutoff_p.adj &
                                                    diff_meth_RBc1c2 > cutoff_diff)), CpGs.i)
        CpGs.i_hyperinC2 <- intersect(names(which(p.adj_wilcox_meth_RBc1c2 < cutoff_p.adj &
                                                    diff_meth_RBc1c2 < -cutoff_diff)), CpGs.i)
      } else if (is.na(cutoff_hypo)) {
        CpGs.i_hyperinC1 <- intersect(names(which(p.adj_wilcox_meth_RBc1c2 < cutoff_p.adj &
                                                    diff_meth_RBc1c2 > cutoff_diff &
                                                    avg_meth_RBc1 > cutoff_hyper)), CpGs.i)
        CpGs.i_hyperinC2 <- intersect(names(which(p.adj_wilcox_meth_RBc1c2 < cutoff_p.adj &
                                                    diff_meth_RBc1c2 < -cutoff_diff &
                                                    avg_meth_RBc2 > cutoff_hyper)), CpGs.i)
      } else if (is.na(cutoff_hyper)){
        CpGs.i_hyperinC1 <- intersect(names(which(p.adj_wilcox_meth_RBc1c2 < cutoff_p.adj &
                                                    diff_meth_RBc1c2 > cutoff_diff &
                                                    avg_meth_RBc2 < cutoff_hypo)), CpGs.i)
        CpGs.i_hyperinC2 <- intersect(names(which(p.adj_wilcox_meth_RBc1c2 < cutoff_p.adj &
                                                    diff_meth_RBc1c2 < -cutoff_diff &
                                                    avg_meth_RBc1 < cutoff_hypo)), CpGs.i)
      }
      length(CpGs.i_hyperinC1)
      length(CpGs.i_hyperinC2)
      CpGs.input <- c(CpGs.i_hyperinC1, CpGs.i_hyperinC2)

      title_param <- paste(name_CpGs.i, ", ", "n = ", length(CpGs.input), "\n",
                           "p.adj", "=", cutoff_p.adj, ", diff", "=", cutoff_diff, ", hypo", "=", cutoff_hypo, ", hyper", "=", cutoff_hyper, sep="")
      filename_param <- gsub("=| ", "", gsub("\n|, ", "_", title_param))
    }

    # CpGs.input
    if (length(CpGs.input) == 0) {
      stop ("No probe meets the criteria")
    } else {
      ## annot CpGs and gene
      # colnames(as.data.frame(annotCpGs))
      annotCpGs.input <- as.data.frame(annotCpGs)[CpGs.input, c("chr", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group") ]
      annotCpGs.input <- cbind.data.frame(annotCpGs.input, statinfo_meth_df[CpGs.input, ])

      if (relatedtogene) {
        CpGs.input <- rownames(annotCpGs.input)[which(annotCpGs.input$UCSC_RefGene_Name != "")]
        annotCpGs.input <- as.data.frame(annotCpGs)[CpGs.input, c("chr", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group") ]
        annotCpGs.input <- cbind.data.frame(annotCpGs.input, statinfo_meth_df[CpGs.input, ])

        title_param <- paste(name_CpGs.i, ", related to gene, ", "n = ", length(CpGs.input), "\n",
                             "p.adj", "=", cutoff_p.adj, ", diff", "=", cutoff_diff, ", hypo", "=", cutoff_hypo, ", hyper", "=", cutoff_hyper, sep="")
        filename_param <- gsub("=| ", "", gsub("\n|, ", "_", title_param))
      }

      dir.create(paste(getwd(), "/", filename_param, sep=""))
      if (length(which(annotCpGs.input$diff_meth_RBc1c2 < 0)) != 0) {
        ## annot CpG
        annotCpGs.input_hyperinC2 <- annotCpGs.input[which(annotCpGs.input$diff_meth_RBc1c2 < 0), ]
        write.table(annotCpGs.input_hyperinC2, file=paste(getwd(), "/", filename_param,  "/annotCpG_hyperinC2_", filename_param, ".csv", sep=""), sep=";", col.names = NA, na = "")
        ## annot gene
        gene <- unlist(lapply(as.list(annotCpGs.input_hyperinC2$UCSC_RefGene_Name), function(zz) unique(unlist(strsplit(zz, ";"))) ))
        gene_hyperinC2 <- as.data.frame(table(gene))
        gene_hyperinC2 <- cbind.data.frame(gene_hyperinC2, df_diffgene.sig[as.vector(gene_hyperinC2$gene), c(1:2, 5)] ) %>% arrange(logFC)
        write.table(gene_hyperinC2, file=paste(getwd(), "/", filename_param,  "/gene_hyperinC2_", filename_param, ".csv", sep=""), sep=";", row.names = F, na = "")
      }
      if (length(which(annotCpGs.input$diff_meth_RBc1c2 > 0)) != 0) {
        ## annot CpG
        annotCpGs.input_hyperinC1 <- annotCpGs.input[which(annotCpGs.input$diff_meth_RBc1c2 > 0), ]
        write.table(annotCpGs.input_hyperinC1, file=paste(getwd(), "/", filename_param,  "/annotCpG_hyperinC1_", filename_param, ".csv", sep=""), sep=";", col.names = NA, na = "")
        ## annot gene
        gene <- unlist(lapply(as.list(annotCpGs.input_hyperinC1$UCSC_RefGene_Name), function(zz) unique(unlist(strsplit(zz, ";"))) ))
        gene_hyperinC1 <- as.data.frame(table(gene))
        gene_hyperinC1 <- cbind.data.frame(gene_hyperinC1, df_diffgene.sig[as.vector(gene_hyperinC1$gene), c(1:2, 5)] ) %>% arrange(desc(logFC))
        write.table(gene_hyperinC1, file=paste(getwd(), "/", filename_param, "/gene_hyperinC1_", filename_param, ".csv", sep=""), sep=";", row.names = F, na = "")
      }

      ##
      if (plot) {
        input_meth <- meth[CpGs.input, samples.i]

        pdf( paste("plot_", filename_param, ".pdf", sep=""), width=12, height=8, paper="a4r")

        ##
        type_samples_factor <- factor(type_samples.i, levels = c("RBC1", "RBC2woMYCN", "RBC2MYCN"))
        color_type_samples <- c("goldenrod", "cornflowerblue", "darkred")[type_samples_factor]
        palette.i <- c("goldenrod", "cornflowerblue", "darkred")
        ColSideColor.i <- palette.i[type_samples_factor]

        HeatColor.i <- myPalette(low="blue", high="yellow", mid="white")
        # HeatColor.i <- colorRamps::blue2yellow(100)


        ## densityPlot of meth value
        densityPlot(input_meth,
                    main= title_param,
                    pal = palette.i,
                    sampGroups = type_samples_factor,
                    legend = F)

        ## Boxplot - nCpGs in each 0.1 interval
        df_percentage_stat <- matrix(NA, nrow=11, ncol=ncol(input_meth))
        colnames(df_percentage_stat) <- colnames(input_meth)

        df_percentage_stat_0.2 <- matrix(NA, nrow=6, ncol=ncol(input_meth))
        colnames(df_percentage_stat_0.2) <- colnames(input_meth)

        for (column.i in colnames(input_meth) ) {
          df_percentage_stat[1, column.i] <- length(which(input_meth[, column.i] >= 0 & input_meth[, column.i] < 0.1))
          df_percentage_stat[2, column.i] <- length(which(input_meth[, column.i] >= 0.1 & input_meth[, column.i] < 0.2))
          df_percentage_stat[3, column.i] <- length(which(input_meth[, column.i] >= 0.2 & input_meth[, column.i] < 0.3))
          df_percentage_stat[4, column.i] <- length(which(input_meth[, column.i] >= 0.3 & input_meth[, column.i] < 0.4))
          df_percentage_stat[5, column.i] <- length(which(input_meth[, column.i] >= 0.4 & input_meth[, column.i] < 0.5))
          df_percentage_stat[6, column.i] <- length(which(input_meth[, column.i] >= 0.5 & input_meth[, column.i] < 0.6))
          df_percentage_stat[7, column.i] <- length(which(input_meth[, column.i] >= 0.6 & input_meth[, column.i] < 0.7))
          df_percentage_stat[8, column.i] <- length(which(input_meth[, column.i] >= 0.7 & input_meth[, column.i] < 0.8))
          df_percentage_stat[9, column.i] <- length(which(input_meth[, column.i] >= 0.8 & input_meth[, column.i] < 0.9))
          df_percentage_stat[10, column.i] <- length(which(input_meth[, column.i] >= 0.9 & input_meth[, column.i] <= 1))
          df_percentage_stat[11, column.i] <- nrow(input_meth)

          df_percentage_stat_0.2[1, column.i] <- length(which(input_meth[, column.i] >= 0 & input_meth[, column.i] < 0.2))
          df_percentage_stat_0.2[2, column.i] <- length(which(input_meth[, column.i] >= 0.2 & input_meth[, column.i] < 0.4))
          df_percentage_stat_0.2[3, column.i] <- length(which(input_meth[, column.i] >= 0.4 & input_meth[, column.i] < 0.6))
          df_percentage_stat_0.2[4, column.i] <- length(which(input_meth[, column.i] >= 0.6 & input_meth[, column.i] < 0.8))
          df_percentage_stat_0.2[5, column.i] <- length(which(input_meth[, column.i] >= 0.8 & input_meth[, column.i] <= 1))
          df_percentage_stat_0.2[6, column.i] <- nrow(input_meth)

        }

        #
        df_percentage_stat <- t(df_percentage_stat)
        df_percentage_stat <- cbind.data.frame(CpGsType = rep(name_CpGs.i, nrow(df_percentage_stat)),
                                               SampleID = samples.i,
                                               SampleType = type_samples_factor,
                                               df_percentage_stat)
        colnames(df_percentage_stat) <- c("CpGsType", "SampleID", "SampleType",
                                          "meth_0_0.1", "meth_0.1_0.2", "meth_0.2_0.3", "meth_0.3_0.4", "meth_0.4_0.5",
                                          "meth_0.5_0.6", "meth_0.6_0.7", "meth_0.7_0.8", "meth_0.8_0.9", "meth_0.9_1",
                                          "total")
        long_df_percentage_stat <- df_percentage_stat %>% gather(interval, n_CpGs, meth_0_0.1:meth_0.9_1)

        ggbox_percent <- ggboxplot(long_df_percentage_stat, x = "SampleType", y = "n_CpGs",
                                   xlab = "", ylab = "",
                                   title = paste("nCpGs in each 0.1 interval: ", gsub("\n", ", ", title_param), sep=""),
                                   color = "SampleType", palette = palette.i,
                                   x.text.angle = 30,
                                   # facet.by = "interval",
                                   add = "jitter")
        print(ggbox_percent + facet_wrap(. ~ interval, nrow=2, scales="free") + stat_compare_means(comparisons = list(c("RBC1", "RBC2woMYCN"))))

        #
        df_percentage_stat_0.2 <- t(df_percentage_stat_0.2)
        df_percentage_stat_0.2 <- cbind.data.frame(CpGsType = rep(name_CpGs.i, nrow(df_percentage_stat_0.2)),
                                                   SampleID = samples.i,
                                                   SampleType = type_samples_factor,
                                                   df_percentage_stat_0.2)
        colnames(df_percentage_stat_0.2) <- c("CpGsType", "SampleID", "SampleType",
                                              "meth_0_0.2", "meth_0.2_0.4", "meth_0.4_0.6", "meth_0.6_0.8", "meth_0.8_1",
                                              "total")
        long_df_percentage_stat_0.2 <- df_percentage_stat_0.2 %>% gather(interval, n_CpGs, meth_0_0.2:meth_0.8_1)

        ggbox_percent_0.2 <- ggboxplot(long_df_percentage_stat_0.2, x = "SampleType", y = "n_CpGs",
                                       xlab = "", ylab = "",
                                       title = paste("nCpGs in each 0.2 interval: ", gsub("\n", ", ", title_param), sep=""),
                                       color = "SampleType", palette = palette.i,
                                       x.text.angle = 30,
                                       # facet.by = "interval",
                                       add = "jitter")
        print(ggbox_percent_0.2 + facet_wrap(. ~ interval, nrow=1, scales="free") + stat_compare_means(comparisons = list(c("RBC1", "RBC2woMYCN"))))

        ## statistical test for nCpGs in 0.1 interval
        colnames(df_percentage_stat)

        colMeans(df_percentage_stat[samples_RBc1, 4:13])
        colMeans(df_percentage_stat[samples_RBc2, 4:13])
        pvalues.t <- sapply(4:13, function(z) t.test(df_percentage_stat[samples_RBc1, z], df_percentage_stat[samples_RBc2, z])$p.value )
        pvalues.wilcox <- sapply(4:13, function(z) wilcox.test(df_percentage_stat[samples_RBc1, z], df_percentage_stat[samples_RBc2, z])$p.value )

        stat_df_percentage_stat <- cbind(NA, NA, NA,
                                         rbind.data.frame(round(colMeans(df_percentage_stat[samples_RBc1, 4:13]), 1),
                                                          round(colMeans(df_percentage_stat[samples_RBc2, 4:13]), 1),
                                                          round(pvalues.t, 3),
                                                          round(pvalues.wilcox, 3) ),
                                         NA)
        colnames(stat_df_percentage_stat) <- colnames(df_percentage_stat)
        rownames(stat_df_percentage_stat) <- c("mean_C1", "mean_C2", "t.test", "wilcox.test")
        output_df_interval <- rbind.data.frame(df_percentage_stat, stat_df_percentage_stat)

        write.table(output_df_interval, file = paste(getwd(), "/", filename_param, "/table_CpG_number_interval0.1_,", filename_param, ".csv", sep=""), col.names = NA, sep=";")

        ## Boxplot of average meth value
        input_boxplot <- cbind.data.frame(avg_meth = colMeans(input_meth), type_samples = type_samples_factor)
        ggbox <- ggboxplot(input_boxplot, y="avg_meth", x="type_samples",
                           title = gsub("\n", ", ", title_param), ylab = "overall meth value", xlab = "",
                           #facet.by = "type.CpGs",
                           color = "type_samples", palette = palette.i,
                           orientation = "horizontal",
                           add = "jitter") + stat_compare_means(comparisons = list(c("RBC1", "RBC2woMYCN")))
        print(ggbox)

        ## Heatmap
        if (nrow(input_meth) <= 25000 ){
          heatmap.2(input_meth,
                    # main = nrow(input_heatmap),
                    main = title_param,
                    distfun = function(z) as.dist(1-cor(t(z))),
                    hclustfun = function(z) hclust(z, method = "complete"),
                    ColSideColors = ColSideColor.i,
                    col = HeatColor.i,
                    trace="none")
        }

        ##
        dev.off()
      }
    }
  }


  relatedtogene = F
  ##
  # dir.i <- paste(dir_res, "/padj0.05_diff0.2", sep="")
  # dir.create(dir.i)
  # setwd(dir.i)
  # setwd(dir_res)
  list_CpGs.i <- list(CpGs_island, CpGs_shore, CpGs_shelf, CpGs_sea, CpGs_outsideisland, CpGs_all)
  list_Name_CpGs.i <- list("CpGs_island", "CpGs_shore", "CpGs_shelf", "CpGs_sea", "CpGs_outsideisland", "CpGs_all")

  lapply(1:length(list_CpGs.i), function(zz) { fun_plot_meth(list_CpGs.i[[zz]], list_Name_CpGs.i[[zz]], cutoff_p.adj = 0.05, cutoff_diff = 0.2, samples.i = samples.i, type_samples.i = type_samples.i, relatedtogene = relatedtogene) })



}


