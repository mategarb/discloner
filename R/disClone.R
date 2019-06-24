#############
# DISCLONER #
#############
#
#   Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

disClone<-function(vaf, method="GMM"){


  clstr=c()
  clone_num=c()

  ### gmm EM
  if(method=="GMM"){
    ## mclust
    ## Density estimation via Gaussian finite mixture modeling
    mod <- densityMclust(as.numeric(vaf))

    clstr=mod$classification
    clone_num=mod$G
  }



  ## distance based
  if(method=="distpamk"){
    vaf=sort(vaf[sample(length(vaf))]) #sort(aVAF)
    mat=outer(vaf, vaf, "-")

    # option 1 pamk
    pamk.best <- pamk(mat)
    predNum2=pamk.best$nc

    clstr=pamk.best$pamobject$clustering
    clone_num=pamk.best$nc
  }

  if(method=="distap"){
    #option 2 apcluster
    vaf=sort(vaf[sample(length(vaf))]) #sort(aVAF)
    mat=outer(vaf, vaf, "-")
    apres <- apcluster(negDistMat(r=2), mat, details=TRUE, q=0.1)
    predNum=length(apres@clusters)

    lst=unlist(lapply(apres@clusters,length))
    cc=c()
    for(l in 1:length(lst)){
      cc=c(cc,rep(l, lst[l]))
    }

    clstr=cc
    clone_num=predNum
  }

  ### kmeans 1D
  if(method=="kmeans"){
    res <- Ckmeans.1d.dp(vaf, k=c(2,200))
    kopt <- length(res$centers)

    clstr=res$cluster
    clone_num=kopt
  }

  ## kmedians 1D
  if(method=="kmedians"){
    vaf=as.numeric(as.matrix(vcfData2$`Variant_allele_ratio%`))
    res <- Ckmedian.1d.dp(vaf, k=c(2,20))
    kopt <- length(res$centers)

    clstr=res$cluster
    clone_num=kopt
  }

  aa=data.frame(vaf, paste0("Clone ",as.character(unlist(clstr))))
  colnames(aa)<-c("vaf","Clones")

  #density plot
  plot_dens<-ggplot(aa, aes(vaf, fill = Clones)) +
    geom_density(alpha = 0.8) +
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_fill_brewer(palette="Set3")

  #histogram
  plot_hist<-ggplot(aa, aes(vaf, fill = Clones)) +
    geom_histogram(alpha = 0.8, aes(y = ..density..),binwidth=0.01, position = 'identity')+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_fill_brewer(palette="Set3")

  return(list(clones=clone_num, clusters=clstr, vaf=vaf, hist=plot_hist, dens=plot_dens))

}
