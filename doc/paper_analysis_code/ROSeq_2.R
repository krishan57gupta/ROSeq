ROSeq<-function(countData, condition, numCores, new_step)
{
  temp_data<-apply(countData, 2, function(x) as.numeric(x))
  colnames(temp_data)<-colnames(countData)
  rownames(temp_data)<-rownames(countData)
  countData<-temp_data
  labels<-unique(condition)
  cOne<-colnames(countData)[which(condition %in% labels[1])]
  cTwo<-colnames(countData)[which(condition %in% labels[2])]
  scgroups<-c(cOne,cTwo)
  geneIndex<-seq_len(nrow(countData))
  results<-list()
  plots_data_1=list()
  plots_data_2=list()
  for(gene in geneIndex)
  {
    # print(gene)
    # print(length(geneIndex))
    combine_result=initiateAnalysis(gene=gene, scdata=countData, scgroups=scgroups, classOne=cOne, classTwo=cTwo ,new_step)
    # print(combine_result)
    results[[gene]]<-combine_result$data
    plots_data_1[[gene]]<-combine_result$plots_data_1
    plots_data_2[[gene]]<-combine_result$plots_data_2
    # print(results[[gene]])
  }
  # results <- pbmcapply::pbmclapply(geneIndex, initiateAnalysis, scdata=countData, scgroups=scgroups, classOne=cOne, classTwo=cTwo, mc.cores=numCores ,new_step)
  pVals<-unlist(lapply(results,function(x) x[12]))
  log2FC<-unlist(lapply(results,function(x) x[13]))
  others_info=matrix(0,nrow=length(pVals),ncol=11)
  others_info[,1]<-unlist(lapply(results,function(x) x[1]))
  others_info[,2]<-unlist(lapply(results,function(x) x[2]))
  others_info[,3]<-unlist(lapply(results,function(x) x[3]))
  others_info[,4]<-unlist(lapply(results,function(x) x[4]))
  others_info[,5]<-unlist(lapply(results,function(x) x[5]))
  others_info[,6]<-unlist(lapply(results,function(x) x[6]))
  others_info[,7]<-unlist(lapply(results,function(x) x[7]))
  others_info[,8]<-unlist(lapply(results,function(x) x[8]))
  others_info[,9]<-unlist(lapply(results,function(x) x[9]))
  others_info[,10]<-unlist(lapply(results,function(x) x[10]))
  others_info[,11]<-unlist(lapply(results,function(x) x[11]))
  pAdj<-stats::p.adjust(pVals, method = "fdr")
  results<-cbind("pVals"=pVals,"pAdj"=pAdj,"log2FC"=log2FC,"Wald_T"=others_info[,11],
                 "F_a"=others_info[,1],"F_b"=others_info[,2],"F_A"=others_info[,3],"F_n_bins"=others_info[,4],"F_R2"=others_info[,5],
                 "S_a"=others_info[,6],"S_b"=others_info[,7],"S_A"=others_info[,8], "S_n_bins"=others_info[,9],"S_R2"=others_info[,10])
  rownames(results)<-rownames(countData)
  print("Done")
  return(list("stats"=results,"plots_1"=plots_data_1,"plots_2"=plots_data_2))
}
initiateAnalysis<-function(gene, scdata, scgroups, classOne, classTwo,new_step)
{
  sp<-scdata[gene, ]
  spOne<-scdata[gene, which(scgroups%in%classOne)]
  spTwo<-scdata[gene, which(scgroups%in%classTwo)]
  geneStats<-getDataStatistics(sp, spOne, spTwo)
  results_groupOne <-findParams(ds=spOne, geneStats,gene=gene,g=1,new_step)
  results_groupTwo <-findParams(ds=spTwo, geneStats,gene=gene,g=2,new_step)
  T<-tryCatch({computeDEG(results_groupOne$data, results_groupTwo$data)}, warning=function(w) {NA}, error=function(esp) {NA})
  pValues<-stats::pchisq(T, df=2, lower.tail=FALSE)
  combinedResults<-list("data"=c(results_groupOne$data[1], results_groupOne$data[2], results_groupOne$data[3], results_groupOne$data[4],
                                 results_groupOne$data[5], results_groupTwo$data[1], results_groupTwo$data[2], results_groupTwo$data[3], 
                                 results_groupTwo$data[4], results_groupTwo$data[5], T, pValues, geneStats[7]),"plots_data_1"=results_groupOne$plots_data,
                        "plots_data_2"=results_groupTwo$plots_data)
  return(combinedResults)
}
getDataStatistics<-function(sp, spOne, spTwo)
{
  maxds<-max(sp)
  minds<-min(sp)
  meands<-mean(sp)
  stdds<-stats::sd(sp)
  if(minds==maxds){
    maxds<-minds+.001
    meands<-minds+.0005
    stdds<-.0001
  }
  ceilds<-ceiling((maxds-meands)/stdds)
  floords<-floor((minds-meands)/stdds)
  f_m=mean(spOne)
  s_m=mean(spTwo)
  if(s_m==0)
  {
    s_m=0.0000001
  }
  if(f_m==0)
  {
    f_m=0.0000001
  }
  log2FC<-log2(f_m/s_m)
  geneStats<-c(maxds, minds, meands, stdds, ceilds, floords,log2FC)
  return (geneStats)
}

findParams<-function(ds, geneStats,gene=62,g=1,new_step)
{
  ############################################### most  imporetant factor if set very smaller then no error ###############################
  ################################################### as it set tradeoff between optimize value a,b and rank, as a,b in power of rank #####
  meands<-geneStats[3]
  stdds<-geneStats[4]
  ceilds<-geneStats[5]
  floords<-geneStats[6]
  step<-new_step
  if((ceilds-floords)/1000>step)
    step=(ceilds-floords)/1000
  binNumber<-length(seq(floords, ceilds-step, step))
  rs<-c(rep(NA), binNumber)
  count<-1
  for(i in seq(floords, ceilds-step, step)){
    LL<- meands+i*stdds
    UL <-meands+(i+step)*stdds
    rs[count]<-length(intersect(which(ds<UL), which(ds>=LL)))
    count<-count+1
  }
  if(sum(is.na(rs))>0){
    print("nul values found")
    print(rs)
    rs<-rs[!is.na(rs)]
  }
  fds<-rs
  number_of_bins<-length(fds)
  rank<-seq_len(number_of_bins)
  read_count_sorted<-sort(fds, decreasing=TRUE)
  normalized_read_count_sorted<-read_count_sorted/sum(read_count_sorted)
  model<-stats::optim(par = c(0.25, 3), minimizeNLL, r=rank, readCount=normalized_read_count_sorted, method = "L-BFGS-B", lower=c(-50,-50), upper=c(50,50))
  a<-model$par[1]
  b<-model$par[2]
  A<-1/sum((number_of_bins+1-rank)^b/(rank^a))
  f<-A*((number_of_bins+1-rank)^b)/(rank^a)
  if(sum(is.na(f))>0){
    print("finding very high or low value")
    print(f[seq_len(10)])}
  SS_res<-sum((normalized_read_count_sorted-f)^2)
  SS_tot<-sum((normalized_read_count_sorted-mean(normalized_read_count_sorted))^2)
  R2<-as.numeric(1-SS_res/SS_tot)
  A<-as.numeric(A)
  results<-list("data"=c(a, b, A, number_of_bins, R2),"plots_data"=list("normalized_read_count_sorted"=normalized_read_count_sorted,
                                                                        "f"=f,"rank"=rank,"ds"=ds))
  return(results)
}


minimizeNLL<-function(coefficients, r, readCount)
{
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  sumReadCount<-sum(readCount)
  A <-1/sum(((N+1-r)^b)/(r^a))
  NLL<-a*sum(readCount*log(r)) - b*sum(readCount*log(N+1-r)) - sumReadCount*log(A)
  NLL<-as.numeric(NLL)
  # print(NLL)
  return (NLL)
}


computeDEG<-function(results_1, results_2)
{
  I_1<-getI(results_1)
  I_2<-getI(results_2)
  I_1<-as.numeric(I_1)
  I_2<-as.numeric(I_2)
  I1<-matrix(c(I_1[1], I_1[2], I_1[3], I_1[4]), nrow = 2, ncol=2)
  I2<-matrix(c(I_2[1], I_2[2], I_2[3], I_2[4]), nrow = 2, ncol=2)
  V1<-solve(I1, tol = 1e-20)
  V2<-solve(I2, tol = 1e-20)
  m<-results_1[4]
  n<-results_2[4]
  a1<-results_1[1]
  b1<-results_1[2]
  a2<-results_2[1]
  b2<-results_2[2]
  w<-n/(m + n)
  T<-(m*n/(m+n))* t(matrix(c(a1-a2,b1-b2), nrow=2, ncol=1)) %*% solve(w *V1 + (1-w)* V2, tol = 1e-20) %*% matrix(c(a1-a2,b1-b2), nrow=2, ncol=1)
  return(T)
}


getI<-function(results)
{
  rank<-seq_len(results[4])
  coefficients<-c(results[1], results[2])
  u1<-getu1(coefficients, rank)
  v<-getv(coefficients, rank)
  u2<-getu2(coefficients, rank)
  du1da<-getdu1da(coefficients, rank)
  du1db<-getdu1db(coefficients, rank)
  du2da<-getdu2da(coefficients, rank)
  du2db<-getdu2db(coefficients, rank)
  dvda<-getdvda(coefficients, rank)
  dvdb<-getdvdb(coefficients, rank)
  d2logAda2<-getd2logAda2( u1, v, du1da, dvda)
  d2logAdb2<-getd2logAdb2( u2, v, du2db, dvdb)
  d2logAdbda<-getd2logAdbda( u1, v, du1db, dvdb)
  d2logAdadb<-getd2logAdadb( u2, v, du2da, dvda)
  I<-c(-d2logAda2, -d2logAdadb, -d2logAdbda, -d2logAdb2)
  return(I)
}


getu1<-function(coefficients, r)
{
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-log(r)
  den1<-r^a
  u1<-sum(num1*num2/den1)
  return(u1)
}


getv<-function( coefficients, r)
{
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  den1<-r^a
  v<-sum(num1/den1)
  return(v)
}


getu2<-function(coefficients, r)
{
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-log(N+1-r)
  den1<-r^a
  u2<-(-sum(num1*num2/den1))
  return(u2)
}


getdu1da<-function(coefficients, r)
{
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-(log(r))^2
  den1<-r^a
  du1da<- -sum(num1*num2/den1)
  return (du1da)
}


getdu1db<-function(coefficients, r)
{
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-log(r)
  num3<-log(N+1-r)
  den1<-r^a
  du1db<-sum(num1*num2*num3/den1)
  return(du1db)
}

getdu2da<-function(coefficients, r)
{
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-log(r)
  num3<-log(N+1-r)
  den1<-r^a
  du2da<-sum(num1*num2*num3/den1)
  return(du2da)
}


getdu2db<-function( coefficients, r)
{
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-(log(N+1-r))^2
  den1<-r^a
  du2db<-(-sum(num1*num2/den1))
  return(du2db)
}


getdvda<-function(coefficients, r)
{
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-log(r)
  den1<-r^a
  dvda<-(-sum(num1*num2/den1))
  return(dvda)
}


getdvdb<-function( coefficients, r)
{
  a<-coefficients[1]
  b<-coefficients[2]
  N<-length(r)
  num1<-(N+1-r)^b
  num2<-log(N+1-r)
  den1<-r^a
  dvdb<-sum(num1*num2/den1)
  return(dvdb)
}


getd2logAda2<-function(u1, v, du1da, dvda)
{
  num1<-v*du1da
  num2<-u1*dvda
  den1<-v^2
  d2logAda2<-(num1-num2)/den1
  return(d2logAda2)
}


getd2logAdadb<-function( u2, v, du2da, dvda)
{
  num1<-v*du2da
  num2<-u2*dvda
  den1<-v^2
  d2logAdadb<-(num1-num2)/den1
  return(d2logAdadb)
}


getd2logAdb2<-function( u2, v, du2db, dvdb)
{
  num1<-v*du2db
  num2<-u2*dvdb
  den1<-v^2
  d2logAdb2<-(num1-num2)/den1
  return(d2logAdb2)
}

getd2logAdbda<-function( u1, v, du1db, dvdb)
{
  num1<-v*du1db
  num2<-u1*dvdb
  den1<-v^2
  d2logAdbda<-(num1-num2)/den1
  return(d2logAdbda)
}

