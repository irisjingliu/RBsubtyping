#' Plot Figure 6
#' @export
Figure6c <- function(){
  require(ggpubr)

  colnames(tff1_initial)
  colnames(tff1_metastasis)


  tff1.1 <- tff1_initial[, c("ID", "subtype", "TFF1_QS")]
  tff1.2 <- tff1_metastasis[, c("ID", "metastatic", "TFF1_QS_tumor")]
  colnames(tff1.1) <- colnames(tff1.2 ) <- c("ID", "annot", "TFF1_QS")
  tff1 <- rbind.data.frame(tff1.1, tff1.2)
  tff1[which(tff1$annot == "yes"), "annot"] <- "Metastatic"
  tff1[which(tff1$annot == "No"), "annot"] <- "Non metastatic"
  tff1[which(tff1$annot == "1"), "annot"] <- "Subtype 1"
  tff1[which(tff1$annot == "2"), "annot"] <- "Subtype 2"

  tff1$series <- c(rep("Intial series", nrow(tff1.1)), rep("HPRF series", nrow(tff1.2)) )

  dev.new()
  p1 = ggboxplot(tff1, x="annot", y="TFF1_QS", facet.by = "series",
            color="annot", palette = c("goldenrod", "cornflowerblue", "black", "grey"),
            add = "jitter" )

  p2 = ggboxplot(tff1, x="annot", y="TFF1_QS",
            xlab = "",
            legend.title = "",
            x.text.angle = 30,
            color="annot", palette = c("goldenrod", "cornflowerblue", "black", "grey"),
            add = "jitter") +
    geom_hline(yintercept = 50,
               linetype = 2)

  print(p2)
  dev.print(png, "Figure6c_TFF1_QS_boxplot_merged.png", res=300, height=1200, width=1500)
  dev.off()

  write.table(tff1, file = "data_for_figure6.csv", sep=";", row.names = F)

}


