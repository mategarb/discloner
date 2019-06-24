## input ##
## mutation_nam, vaf, patient_id

#earlyDrivers<-function(vaf, mut, pat_id, dist_method="euclidean", tm="ward.D", test="barnard", signQuartile=0.05, info=FALSE, infoVec, minMut=10, maxMut=500){
earlyDrivers<-function(vaf, mut, pat_id, VAFthr = 0.65, info=FALSE, infoVec){

  if(info){ ##just in case of adding aditional info about mutation e.g. chromosome or unique name
    mut=paste0(mut,"-",infoVec)
  }

  tcgaID <- unique(as.character(pat_id)) ##unique patients
  # remove patients with only 1 mutation
  if(length(which(tcgaID %in% names(which(table(as.character(pat_id))==1))))==0){
  tcgaID2<-tcgaID
  }else{
  nn2=which(tcgaID %in% names(which(table(as.character(pat_id))==1)))
  tcgaID2<-tcgaID[-nn2]
  }
  ### declare variables ###
  clones_dc=c()
  clones_dc2=c()
  proped=c()

  out_early=list()
  out_late=list()
  out_early_clonal=list()
  out_early_subclonal=list()
  out_late_clonal=list()
  out_late_subclonal=list()
  early_clonal_nr <- c()
  early_subclonal_nr <- c()
  late_clonal_nr <- c()
  late_subclonal_nr <- c()
  prop_cs <- c()
  mut_per_pat <- list()
  vaf_per_pat <- list()

  ##########################


  # loop over patients
  for(i in 1:length(tcgaID2)){

    dfout_early=list()
    dfout_late=list()

    pat_ind <- which(as.character(pat_id) %in% tcgaID2[i]) # choose i patient
    mut2 <- mut[pat_ind]
    # consider all muts after 0.5 as duplicated
    # i = 37 is long
    vaf2 <- as.numeric(vaf[pat_ind])

    # shift everything above VAFthr
    if(length(which(vaf2 > VAFthr)) > 0){
      vaf2[which(vaf2 > VAFthr)]<-vaf2[which(vaf2 > VAFthr)]/2
    }

    inp <- data.frame(mut2, vaf2)

    # find clones using gmm
    mod <- densityMclust(as.numeric(inp$vaf2))
    clstr <- mod$classification
    clone_num <- mod$G

    clones_dc2[i] <- clone_num/length(inp$vaf2)
    clones_dc[i] <- clone_num

    # if low number of clones was detected, then skip it!
    if(clones_dc2[i]>0.5){

      mut_per_pat[[i]] <- NaN
      vaf_per_pat[[i]] <- NaN

      out_early[[i]] <- NaN
      out_late[[i]] <- NaN

      out_early_subclonal[[i]] <- NaN
      out_early_clonal[[i]] <- NaN
      out_late_subclonal[[i]] <- NaN
      out_late_clonal[[i]] <- NaN

      early_clonal_nr[i] <- NaN
      early_subclonal_nr[i] <- NaN
      late_clonal_nr[i] <- NaN
      late_subclonal_nr[i] <- NaN

    # if it's ok do it:
    }else{

    mut_per_pat[[i]]<-mut2
    vaf_per_pat[[i]]<-vaf2

    vaf_clone_mean <- c()
    vaf_clone_median <- c()

    #loop over the clones and search for early and late events
    for(j in 1:clones_dc[i]){

      vaf_clone <- vaf2[which(clstr == j)]
      mut_clone <- mut2[which(clstr == j)]

      #in case there are the same mutations with the same VAF on different positions (it may happen as we use merged output from different software)
      if(length(unique(mut_clone))==1 & length(unique(vaf_clone))==1)
      {
        vaf_clone <- unique(vaf2[which(clstr == j)])
        mut_clone <- unique(mut2[which(clstr == j)])
      }

      # check if there are at least two mutations
      if(length(vaf_clone) >= 2){

        vaf_clone_mean[j] <- mean(vaf2[which(clstr == j)])
        vaf_clone_median[j] <- median(vaf2[which(clstr == j)])

        res_clone <- Ckmeans.1d.dp(vaf_clone, k=2, method = "linear")

        vaf_clusters <- res_clone$cluster
        vaf_clust_num <- length(unique(vaf_clusters))
        vaf_clone_clust <- data.frame(vaf_clusters, vaf_clone)

        # give names
        rownames(vaf_clone_clust) <- paste0(mut_clone, "_.",1:length(mut_clone))
        colnames(vaf_clone_clust) <- c("clust", "VAF")

        vaf_agg <- aggregate(vaf_clone_clust["VAF"], by = list(vaf_clusters), FUN = mean)
        vgr <- vaf_agg$Group.1[which.max(vaf_agg$VAF)]

        temp_early <- vaf_clone_clust[which(rownames(vaf_clone_clust) %in% rownames(vaf_clone_clust)[vaf_clone_clust$clust==vgr]),2]
        names(temp_early) <- as.character(gsub("\\_..*","",rownames(vaf_clone_clust)[vaf_clone_clust$clust==vgr]))
        dfout_early[[j]] <- temp_early

        temp_late <- vaf_clone_clust[which(rownames(vaf_clone_clust) %in% rownames(vaf_clone_clust)[vaf_clone_clust$clust!=vgr]),2]
        names(temp_late) <- as.character(gsub("\\_..*","",rownames(vaf_clone_clust)[vaf_clone_clust$clust!=vgr]))
        dfout_late[[j]]<-temp_late

      }else{
        dfout_early[[j]]=""
        dfout_late[[j]]=""
      }
    }

    out_early[[i]] <- unlist(dfout_early)
    out_late[[i]] <- unlist(dfout_late)

    #clonal - the closest to 0.5
    out_early_clonal[[i]] <- unlist(dfout_early[unique(c(which.min(abs(vaf_clone_mean-0.5)),which.max(vaf_clone_mean)))])
    out_late_clonal[[i]] <- unlist(dfout_late[unique(c(which.min(abs(vaf_clone_mean-0.5)),which.max(vaf_clone_mean)))])

    #subclonal - others, but no clonal

    ## one clone, means that there are no subclones
    if(clones_dc[i]==1){
      out_early_subclonal[[i]] <- NULL
      out_late_subclonal[[i]] <- NULL
    }else{
    out_early_subclonal[[i]] <- unlist(dfout_early[-unique(c(which.min(abs(vaf_clone_mean-0.5)),which.max(vaf_clone_mean)))])
    out_late_subclonal[[i]] <- unlist(dfout_late[-unique(c(which.min(abs(vaf_clone_mean-0.5)),which.max(vaf_clone_mean)))])
    }

    early_clonal_nr[i] <- length(unlist(out_early_clonal[i]))
    early_subclonal_nr[i] <- length(unlist(out_early_subclonal[i]))
    late_clonal_nr[i] <- length(unlist(out_late_clonal[i]))
    late_subclonal_nr[i] <- length(unlist(out_late_subclonal[i]))

    }

  }

  ## probability
  p_clonal = sum(sum(early_clonal_nr, na.rm = T), sum(late_clonal_nr, na.rm = T))/sum(sum(early_clonal_nr, na.rm = T),sum(early_subclonal_nr, na.rm = T), sum(late_clonal_nr, na.rm = T), sum(late_subclonal_nr, na.rm = T))
  p_subclonal = sum(sum(early_subclonal_nr, na.rm = T), sum(late_subclonal_nr, na.rm = T))/sum(sum(early_clonal_nr, na.rm = T),sum(early_subclonal_nr, na.rm = T), sum(late_clonal_nr, na.rm = T), sum(late_subclonal_nr, na.rm = T))

  smstat <- data.frame(early_clonal_nr, late_clonal_nr, early_subclonal_nr, late_subclonal_nr)

  mut_early <- table(as.character(gsub("\\_..*","",names(unlist(out_early)))))
  mut_late <- table(as.character(gsub("\\_..*","",names(unlist(out_late)))))

  mut_late_clonal <- table(as.character(gsub("\\_..*","",names(unlist(out_late_clonal)))))
  mut_early_clonal <- table(as.character(gsub("\\_..*","",names(unlist(out_early_clonal)))))
  mut_late_subclonal <- table(as.character(gsub("\\_..*","",names(unlist(out_late_subclonal)))))
  mut_early_subclonal <- table(as.character(gsub("\\_..*","",names(unlist(out_early_subclonal)))))

  #remove empty names
  if(names(mut_early[1])==""){
    mut_early <- mut_early[-1]
  }
  if(names(mut_late[1])==""){
    mut_late <- mut_late[-1]
  }
  if(names(mut_early_subclonal[1])==""){
    mut_early_subclonal <- mut_early_subclonal[-1]
  }
  if(names(mut_early_clonal[1])==""){
    mut_early_clonal <- mut_early_clonal[-1]
  }
  if(names(mut_late_subclonal[1])==""){
    mut_late_subclonal <- mut_late_subclonal[-1]
  }
  if(names(mut_late_clonal[1])==""){
    mut_late_clonal <- mut_late_clonal[-1]
  }

  #p1=1/mean(cl0)
  pvalbin1 <- c()
  pvalbin2 <- c()
  pvalbin3 <- c()
  pvalbin4 <- c()
  nam <- c()
  gecf <- c()
  glcf <- c()
  gesf <- c()
  glsf <- c()


  for(i in 1:length(mut_early_clonal)){
    #i = which(names(mut_early_clonal)=="KIT")
    gene_early_freq <- unname(mut_early)[which(names(mut_early) == names(mut_early_clonal[i]))] #all drivers early and late
    gene_late_freq <- unname(mut_late)[which(names(mut_late) == names(mut_early_clonal[i]))] #all drivers early and late

    if(length(gene_early_freq)==0){
      gene_early_freq<-0
    }

    if(length(gene_late_freq)==0){
      gene_late_freq<-0
    }

    gene_early_subclonal_freq <- unname(mut_early_subclonal)[which(names(mut_early_subclonal) == names(mut_early_clonal[i]))] #all drivers early and late
    gene_early_clonal_freq <- unname(mut_early_clonal[i]) #all early drivers
    gene_late_subclonal_freq <- unname(mut_late_subclonal)[which(names(mut_late_subclonal) == names(mut_early_clonal[i]))] #all drivers early and late
    gene_late_clonal_freq <- unname(mut_late_clonal)[which(names(mut_late_clonal) == names(mut_early_clonal[i]))] #all drivers early and late

    if(length(gene_early_subclonal_freq)==0){
      gene_early_subclonal_freq<-0
    }

    if(length(gene_early_clonal_freq)==0){
      gene_early_clonal_freq<-0
    }

    if(length(gene_late_subclonal_freq)==0){
      gene_late_subclonal_freq<-0
    }

    if(length(gene_late_clonal_freq)==0){
      gene_late_clonal_freq<-0
    }

    gene_all_freq <- gene_early_freq + gene_late_freq

    #type <- c(earlyFreq, notearlyFreq)
    #observed = c(earlyFreq, as.numeric(as.character(s3[i])))
    #pmut=(earlyFreq)/as.numeric(as.character(s3[i]))
    #expected = c(pmut, 1-pmut)


    ### STATISTIC TESTS ###

    csbin1<-binom.test(gene_early_clonal_freq, gene_all_freq, p = p_clonal, alternative = "greater")
    pvalbin1[i]=p.adjust(csbin1$p.value, "bonferroni")

    csbin2<-binom.test(gene_early_subclonal_freq, gene_all_freq, p = p_subclonal, alternative = "greater")
    pvalbin2[i]=p.adjust(csbin2$p.value, "bonferroni")

    csbin3<-binom.test(gene_late_clonal_freq, gene_all_freq, p = p_clonal, alternative = "greater")
    pvalbin3[i]=p.adjust(csbin3$p.value, "bonferroni")

    csbin4<-binom.test(gene_late_subclonal_freq, gene_all_freq, p = p_subclonal, alternative = "greater")
    pvalbin4[i]=p.adjust(csbin4$p.value, "bonferroni")



    #cs<-barnard.test(gene_early_clonal_freq, gene_early_subclonal_freq, gene_late_clonal_freq, gene_late_subclonal_freq, pooled=TRUE)
    #pval[i]=cs$p.value[2]


   # if(test=="fisher"){
   # csbin<-binom.test(gene_early_clonal_freq+gene_early_subclonal_freq, gene_late_clonal_freq+gene_late_subclonal_freq+gene_early_clonal_freq+gene_early_subclonal_freq)
   # pvalbin[i]=csbin$p.value

   # cs<-fisher.test(matrix(c(gene_early_clonal_freq, gene_late_clonal_freq, gene_early_subclonal_freq, gene_late_subclonal_freq), ncol = 2))
   # pval[i]=cs$p.value
   # }

    nam[i] <- names(mut_early_clonal[i])

    gecf[i] <- gene_early_clonal_freq
    glcf[i] <- gene_late_clonal_freq
    gesf[i] <- gene_early_subclonal_freq
    glsf[i] <- gene_late_subclonal_freq
  }

  dout <- c()
  # prepare output of early drivers
  dout <- data.frame(pvalbin1, pvalbin2, pvalbin3, pvalbin4, nam, gecf, glcf, gesf, glsf, gecf+glcf+gesf+glsf)
  colnames(dout) <- c("pval_EC","pval_ES","pval_LC","pval_LS","gene_name","early_clonal","late_clonal","early_subclonal","late_subclonal","all")

  dout2 <- dout[order(dout$pval_EC),]

  # early clonal VAF
  mut_early_clonal_vaf <- unlist(out_early_clonal)[!is.na(match(gsub("\\_..*","",names(unlist(out_early_clonal))),as.character(dout2$gene_name)))]
  nams_early_clonal <- gsub("\\_..*","",names(mut_early_clonal_vaf))
  vafs_early_clonal <- as.numeric(mut_early_clonal_vaf)
  tab_early_clonal <- data.frame(vafs_early_clonal, nams_early_clonal)
  tab_early_clonal_agg <-aggregate(tab_early_clonal["vafs_early_clonal"], by=list(nams_early_clonal), FUN=mean, na.rm=TRUE)

  # early subclonal VAF
  mut_early_subclonal <- unlist(out_early_subclonal)[!is.na(match(gsub("\\_..*","",names(unlist(out_early_subclonal))),as.character(dout2$gene_name)))]
  nams_early_subclonal <- gsub("\\_..*","",names(mut_early_subclonal))
  vafs_early_subclonal <- as.numeric(mut_early_subclonal)
  tab_early_subclonal <- data.frame(vafs_early_subclonal, nams_early_subclonal)
  tab_early_subclonal_agg <-aggregate(tab_early_subclonal["vafs_early_subclonal"], by=list(nams_early_subclonal), FUN=mean, na.rm=TRUE)

  # late clonal VAF
  mut_late_clonal_vaf <- unlist(out_late_clonal)[!is.na(match(gsub("\\_..*","",names(unlist(out_late_clonal))),as.character(dout2$gene_name)))]
  nams_late_clonal <- gsub("\\_..*","",names(mut_late_clonal_vaf))
  vafs_late_clonal <- as.numeric(mut_late_clonal_vaf)
  tab_late_clonal <- data.frame(vafs_late_clonal, nams_late_clonal)
  tab_late_clonal_agg <-aggregate(tab_late_clonal["vafs_late_clonal"], by=list(nams_late_clonal), FUN=mean, na.rm=TRUE)

  # early subclonal VAF
  mut_late_subclonal <- unlist(out_late_subclonal)[!is.na(match(gsub("\\_..*","",names(unlist(out_late_subclonal))),as.character(dout2$gene_name)))]
  nams_late_subclonal <- gsub("\\_..*","",names(mut_late_subclonal))
  vafs_late_subclonal <- as.numeric(mut_late_subclonal)
  tab_late_subclonal <- data.frame(vafs_late_subclonal, nams_late_subclonal)
  tab_late_subclonal_agg <-aggregate(tab_late_subclonal["vafs_late_subclonal"], by=list(nams_late_subclonal), FUN=mean, na.rm=TRUE)

  # write them DOOOWN
  col_early_clonal_agg <- tab_early_clonal_agg$vafs_early_clonal[match(dout2$gene_name, tab_early_clonal_agg$Group.1)]
  col_early_subclonal_agg <- tab_early_subclonal_agg$vafs_early_subclonal[match(dout2$gene_name, tab_early_subclonal_agg$Group.1)]
  col_late_clonal_agg <- tab_late_clonal_agg$vafs_late_clonal[match(dout2$gene_name, tab_late_clonal_agg$Group.1)]
  col_late_subclonal_agg <- tab_late_subclonal_agg$vafs_late_subclonal[match(dout2$gene_name, tab_late_subclonal_agg$Group.1)]

  # and put to data.frame
  dout3=data.frame(dout2, col_early_clonal_agg,  col_late_clonal_agg, col_early_subclonal_agg, col_late_subclonal_agg)

  newdout3=dout3[,1:4]
  newdout3[dout3[,1:4]>0.05]<-"ns"
  newdout3[dout3[,1:4]<=0.05]<-"*"
  newdout3[dout3[,1:4]<=0.01]<-"**"
  newdout3[dout3[,1:4]<=0.001]<-"***"




  #dout3[,1:4]<-newdout3

  dout4<-data.frame(dout3$gene_name, dout3[1:4], newdout3, dout3[c(6,8,7,9:14)])
  colnames(dout4)<-c("gene_name","pv_ec","pv_es","pv_lc","pv_ls",
                     "s_ec","s_es","s_lc","s_ls",
                     "n_ec","n_es","n_lc","n_ls","n_all",
                     "mvaf_ec","mvaf_es","mvaf_lc","mvaf_ls"
                     )

  return(list(mutList=dout4, clonSubclon=smstat))
}





