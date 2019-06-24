boxED<-function(ed, alld){
  library(forcats)
  ##specific mut
  msel=alld[which(alld$mut %in% as.character(ed)),]
  colnames(msel)<-c("mut","vaf")
  mselagg=aggregate(vaf ~ mut, data = msel, median)

  colfunc <- colorRampPalette(c("dodgerblue4", "lightcyan"))
  cols <- colfunc(length(unique(msel$mut)))[order(ed)]

  so.summary <- function(y) {
    return(data.frame(y=median(y), size=3))
  }

  a1=ggplot(msel, aes(x=fct_reorder(mut, vaf, .fun=median, .desc =TRUE), y=vaf, fill=mut)) +
    xlab("Mutation")+
    ylab("VAF")+
    #geom_boxplot(alpha=0.7)+
    geom_violin(trim=FALSE)+
    coord_flip() +
    #ggtitle(desc[1]) +
    scale_fill_manual(values=cols) +
    geom_hline(yintercept=0.5, linetype="dashed",
               color = "honeydew3", size=1)+
    geom_hline(yintercept=0.25, linetype="dashed",
               color = "honeydew3", size=0.5)+
    geom_hline(yintercept=0.75, linetype="dashed",
               color = "honeydew3", size=0.5)+
    theme_classic()+
    theme(axis.text.x = element_text(size=6, angle = 60, hjust = 1))+
    stat_summary(fun.data=so.summary, geom="point", color="lightsalmon2")

  return(a1)
}
