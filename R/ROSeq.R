##' @title ROSeq - A rank based approach to modeling gene expression with filtered and normalized read count matrix
##' @description Takes in the complete filtered and normalized read count matrix, the location of the two sub-populations and the number of cores to be used
##' @param countData The normalised and filtered, read count matrix, with row names as genes name/ID and column names as sample id/name
##' @param condition Labels for the two sub-populations
##' @param nbits optional, 10 as more accuracy with little slow (due to calculation on big number, Rmpfr library used), default is 0, means no Rmpfr library used for calculation
##' @param numCores The number of cores to be used
##' @return pValues A vector containing FDR adjusted p significance values
##' @examples countData<-matrix(sample(c(seq_len(100)),1000,replace = TRUE),nrow=100,ncol=10)
##' rownames(countData)<-paste("G",seq_len(100),sep="_")
##' colnames(countData)<-paste("S",seq_len(10),sep="_")
##' condition<-c(1,1,1,1,1,2,2,2,2,2)
##' countData<-apply(countData,2,function(x) as.numeric(x))
##' g_keep <- apply(countData,1,function(x) sum(x>0)>5)
##' countData<-edgeR::cpm(countData) # optioanl, can be used other normalization
##' library(ROSeq)
##' output<-ROSeq(countData=countData, condition = condition, nbits=0, numCores=1)
##' output
##'
##' @export ROSeq
ROSeq<-function(countData, condition, nbits=0, numCores)
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
    # results<-list()
    # for(gene in geneIndex)
    # {
    #     print(gene)
    #     print(length(geneIndex))
    #     results[[gene]]<-initiateAnalysis(gene=gene, scdata=countData, scgroups=scgroups, classOne=cOne, classTwo=cTwo, nbits)
    #     print(results[[gene]])
    # }
    results <- pbmcapply::pbmclapply(geneIndex, initiateAnalysis, scdata=countData, scgroups=scgroups, classOne=cOne, classTwo=cTwo, mc.cores=numCores, nbits)
    pVals<-unlist(lapply(results,function(x) x[[12]]))
    pAdj<-stats::p.adjust(pVals, method = "fdr")
    results<-cbind(pVals,pAdj)
    rownames(results)<-rownames(countData)
    return(results)
}

##' @title Computes differential analysis for a given gene
##'
##' @description Takes in the row index which corresponds to a gene and evaluates for differential expression across two cell types.
##' @param gene The row index of the normalised and filtered, read count matrix
##' @param scdata The normalised and filtered, read count matrix
##' @param scgroups The location of the two sub-populations
##' @param classOne The location of the first sub-population, for example, sample names as given as columns names
##' @param classTwo The location of thesecond sub-population, for example, sample names as given as columns names
##' @param nbits number of bits for mpfr function
##' @return combinedResults A vector containing 12 values (gr1: a, g1: b, gr1: A, gr1: number of bins, gr1: R2, gr2: a, gr2: b, gr2: A, gr2: number of bins, gr2: R2, T, p)
initiateAnalysis<-function(gene, scdata, scgroups, classOne, classTwo, nbits)
{
    sp<-scdata[gene, ]
    spOne<-scdata[gene, which(scgroups%in%classOne)]
    spTwo<-scdata[gene, which(scgroups%in%classTwo)]
    geneStats<-getDataStatistics(sp, spOne, spTwo, nbits)
    results_groupOne <-findParams(spOne, geneStats, nbits)
    results_groupTwo <-findParams(spTwo, geneStats, nbits)
    T<-tryCatch({computeDEG(results_groupOne, results_groupTwo, nbits)}, warning=function(w) {NA}, error=function(esp) {NA})
    pValues<-stats::pchisq(T, df=2, lower.tail=FALSE)
    combinedResults<-c(results_groupOne[1], results_groupOne[2], results_groupOne[3], results_groupOne[4], results_groupOne[5], results_groupTwo[1], results_groupTwo[2], results_groupTwo[3], results_groupTwo[4], results_groupTwo[5], T, pValues)
    return(combinedResults)
}

##' @title Evaluates statistics of the read counts corresponding to the gene
##'
##' @description Takes in the complete read count vector corresponding to the gene (sp) and also the data corresponding to the two sub-populations (sp1 and sp2)
##'
##' @param sp The complete (normalized and filtered) read count data corresponding to the gene in question
##' @param spOne The (normalized and filtered) read count data corresponding to the first sub-population
##' @param spTwo The (normalized and filtered) read count data corresponding to the second sub-population
##' @param nbits number of bits for mpfr function
##' @return geneStats A vector containing 7 values corresponding to the gene data (maximum, minimum, mean, standard deviation, upper multiple of standard deviation, lower multiple of standard deviation and log2(fold change))
getDataStatistics<-function(sp, spOne, spTwo, nbits)
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
    geneStats<-c(maxds, minds, meands, stdds, ceilds, floords)
    return (geneStats)
}

##' @title Finds the optimal values of parameters a and b that model the probability distribution of ranks, by Maximising the Log-Likelihood
##'
##' @description Takes in as input the read count data corresponding to one sub-population and the typical gene statistics. 
##' Then it splits the entire range into equally sized bins of size \eqn{k * \sigma}, where k is a scalar with a default value of 0.05, and \eqn{\sigma} is the standard deviation of the pulled expression estimates across the cell-groups. 
##' Each of these bins corresponds to a rank. Therefore, for each group, cell frequency for each bin maps to a rank.  These frequencies are normalized group-wise by dividing by the total cell count within a concerned group.
##' @param ds The (normalized and filtered) read count data corresponding to a sub-population
##' @param nbits number of bits for mpfr function
##' @param geneStats A vector containing 7 values corresponding to the gene data (maximum, minimum, mean, standard deviation, upper multiple of standard deviation, lower multiple of standard deviation and log_{2}(fold change))
##' @return results A vector containing 5 values (a, b, A, number of bins, R2)
findParams<-function(ds, geneStats, nbits)
{
    ############################################### most  imporetant factor if set very smaller then no error ###############################
    ################################################### as it set tradeoff between optimize value a,b and rank, as a,b in power of rank #####
    step<-.05
    meands<-geneStats[3]
    stdds<-geneStats[4]
    ceilds<-geneStats[5]
    floords<-geneStats[6]
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
    if(nbits>0){
        model<-stats::optim(par = c(0.25, 3), minimizeNLL, r=rank, readCount=normalized_read_count_sorted, nbits=nbits, method = "Nelder-Mead")}
    if(nbits<=0){
        model<-stats::optim(par = c(0.25, 3), minimizeNLL, r=rank, readCount=normalized_read_count_sorted, nbits=nbits, method = "L-BFGS-B", lower=c(-50,-50), upper=c(50,50))}
    a<-model$par[1]
    b<-model$par[2]
    if(nbits>0){
        A<-1/sum((number_of_bins+1-rank)^Rmpfr::mpfr(b,nbits)/(rank^Rmpfr::mpfr(a,nbits)))
        f<-A*((number_of_bins+1-rank)^Rmpfr::mpfr(b,nbits))/(rank^Rmpfr::mpfr(a,nbits))}
    if(nbits<=0){
        A<-1/sum((number_of_bins+1-rank)^b/(rank^a))
        f<-A*((number_of_bins+1-rank)^b)/(rank^a)}
    if(sum(is.na(f))>0){
        print("finding very high or low value")
        print(f[seq_len(10)])}
    SS_res<-sum((normalized_read_count_sorted-f)^2)
    SS_tot<-sum((normalized_read_count_sorted-mean(normalized_read_count_sorted))^2)
    R2<-as.numeric(1-SS_res/SS_tot)
    A<-as.numeric(A)
    results<-c(a, b, A, number_of_bins, R2)
    return(results)
}

##' @title Minimizes the Negative Log-Likelihood by iterating across values of parameters a and b
##'
##' @description Takes in as input a vector of values (coefficients), the number of bins and the normalized probability dsitribution of ranks
##' @param coefficients A vector containing two values for a and b
##' @param r The number of bins
##' @param readCount A vector of (normalized) frequency of read counts that fall within each bin
##' @return NLL Negative-Log Likelihood for the input coefficients
##' @param nbits number of bits for mpfr function
##' @seealso \code{\link{findParams}}
minimizeNLL<-function(coefficients, r, readCount, nbits)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    sumReadCount<-sum(readCount)
    if(nbits>0)
    {
        A <-1/sum(((N+1-r)^Rmpfr::mpfr(b,nbits))/(r^Rmpfr::mpfr(a,nbits)))
        NLL<-a*sum(readCount*log(r)) - b*sum(readCount*log(N+1-r)) - sumReadCount*log(A)
    }
    if(nbits<=0)
    {
        A <-1/sum(((N+1-r)^b)/(r^a))
        NLL<-a*sum(readCount*log(r)) - b*sum(readCount*log(N+1-r)) - sumReadCount*log(A)
    }
    NLL<-as.numeric(NLL)
    # print(NLL)
    return (NLL)
}

##' @title Computes differential expression for the gene in question, by comparing the optimal parameters for sub-populations one and two
##' @description  Uses the (asymptotically) optimum two-sample Wald test  based on the MLE of the parameters and its asymptotic variances given by the inverse of the Fisher information matrix
##' @param results_1 A vector corresponding to sub-population one and containing 5 values (a, b, A, number of bins, R2)
##' @param results_2 A vector corresponding to sub-population two and containing 5 values (a, b, A, number of bins, R2)
##' @param nbits number of bits for mpfr function
##' @return T  The Wald test statistic for testing the null hypothesis
##' @seealso \code{\link{getI}}, \code{\link{findParams}}
computeDEG<-function(results_1, results_2, nbits)
{
    I_1<-getI(results_1, nbits)
    I_2<-getI(results_2, nbits)
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

##' @title Computes the Fisher Information Matrix
##' @description The Fisher Information Matrix and its derivatives are essential to evulate the minima of log likelihood
##' @param results A vector containing 5 values (a, b, A, number of bins, R2)
##' @param nbits number of bits for mpfr function
##' @return I  The Fisher Information Matrix
getI<-function(results, nbits)
{
    rank<-seq_len(results[4])
    coefficients<-c(results[1], results[2])
    u1<-getu1(coefficients, rank, nbits)
    v<-getv(coefficients, rank, nbits)
    u2<-getu2(coefficients, rank, nbits)
    du1da<-getdu1da(coefficients, rank, nbits)
    du1db<-getdu1db(coefficients, rank, nbits)
    du2da<-getdu2da(coefficients, rank, nbits)
    du2db<-getdu2db(coefficients, rank, nbits)
    dvda<-getdvda(coefficients, rank, nbits)
    dvdb<-getdvdb(coefficients, rank, nbits)
    d2logAda2<-getd2logAda2( u1, v, du1da, dvda)
    d2logAdb2<-getd2logAdb2( u2, v, du2db, dvdb)
    d2logAdbda<-getd2logAdbda( u1, v, du1db, dvdb)
    d2logAdadb<-getd2logAdadb( u2, v, du2da, dvda)
    I<-c(-d2logAda2, -d2logAdadb, -d2logAdbda, -d2logAdb2)
    return(I)
}

##' @title Computes u1
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param nbits number of bits for mpfr function
##' @param r the rank vector
##' @return u1
getu1<-function(coefficients, r, nbits)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    if(nbits>0)
    {
        num1<-Rmpfr::mpfr(N+1-r,nbits)^Rmpfr::mpfr(b,nbits)
        num2<-log(r)
        den1<-Rmpfr::mpfr(r,nbits)^Rmpfr::mpfr(a,nbits)
    }
    if(nbits<=0)
    {
        num1<-(N+1-r)^b
        num2<-log(r)
        den1<-r^a
    }
    u1<-sum(num1*num2/den1)
    return(u1)
}

##' @title Computes v
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param nbits number of bits for mpfr function
##' @param r the rank vector
##' @return v
getv<-function( coefficients, r, nbits)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    if(nbits>0)
    {
        num1<-Rmpfr::mpfr(N+1-r,nbits)^Rmpfr::mpfr(b,nbits)
        den1<-Rmpfr::mpfr(r,nbits)^Rmpfr::mpfr(a,nbits)
    }
    if(nbits<=0)
    {
        num1<-(N+1-r)^b
        den1<-r^a
    }
    v<-sum(num1/den1)
    return(v)
}

##' @title Computes u2
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param nbits number of bits for mpfr function
##' @param r the rank vector
##' @return u2
getu2<-function(coefficients, r, nbits)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    if(nbits>0)
    {
        num1<-Rmpfr::mpfr(N+1-r,nbits)^Rmpfr::mpfr(b,nbits)
        num2<-log(N+1-r)
        den1<-Rmpfr::mpfr(r,nbits)^Rmpfr::mpfr(a,nbits)
    }
    if(nbits<=0)
    {
        num1<-(N+1-r)^b
        num2<-log(N+1-r)
        den1<-r^a
    }
    u2<-(-sum(num1*num2/den1))
    return(u2)
}

##' @title Finds the first derivative of u1 with respect to a. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param nbits number of bits for mpfr function
##' @param r the rank vector
##' @return du1da
getdu1da<-function(coefficients, r, nbits)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    if(nbits>0)
    {
        num1<-Rmpfr::mpfr(N+1-r,nbits)^Rmpfr::mpfr(b,nbits)
        num2<-(log(r))^2
        den1<-Rmpfr::mpfr(r,nbits)^Rmpfr::mpfr(a,nbits)
    }
    if(nbits<=0)
    {
        num1<-(N+1-r)^b
        num2<-(log(r))^2
        den1<-r^a
    }
    du1da<- -sum(num1*num2/den1)
    return (du1da)
}

##' @title Finds the first derivative of u1 with respect to b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param nbits number of bits for mpfr function
##' @param r the rank vector
##' @return du1db
getdu1db<-function(coefficients, r, nbits)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    if(nbits>0)
    {
        num1<-Rmpfr::mpfr(N+1-r,nbits)^Rmpfr::mpfr(b,nbits)
        num2<-log(r)
        num3<-log(N+1-r)
        den1<-Rmpfr::mpfr(r,nbits)^Rmpfr::mpfr(a,nbits)
    }
    if(nbits<=0)
    {
        num1<-(N+1-r)^b
        num2<-log(r)
        num3<-log(N+1-r)
        den1<-r^a
    }
    du1db<-sum(num1*num2*num3/den1)
    return(du1db)
}

##' @title Finds the first derivative of u2 with respect to a. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param nbits number of bits for mpfr function
##' @param r the rank vector
##' @return du2da
getdu2da<-function(coefficients, r, nbits)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    if(nbits>0)
    {
        num1<-Rmpfr::mpfr(N+1-r,nbits)^Rmpfr::mpfr(b,nbits)
        num2<-log(r)
        num3<-log(N+1-r)
        den1<-Rmpfr::mpfr(r,nbits)^Rmpfr::mpfr(a,nbits)
    }
    if(nbits<=0)
    {
        num1<-(N+1-r)^b
        num2<-log(r)
        num3<-log(N+1-r)
        den1<-r^a
    }
    du2da<-sum(num1*num2*num3/den1)
    return(du2da)
}

##' @title Finds the first derivative of u2 with respect to b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param nbits number of bits for mpfr function
##' @param r the rank vector
##' @return du2db
getdu2db<-function( coefficients, r, nbits)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    if(nbits>0)
    {
        num1<-Rmpfr::mpfr(N+1-r,nbits)^Rmpfr::mpfr(b,nbits)
        num2<-(log(N+1-r))^2
        den1<-Rmpfr::mpfr(r,nbits)^Rmpfr::mpfr(a,nbits)
    }
    if(nbits<=0)
    {
        num1<-(N+1-r)^b
        num2<-(log(N+1-r))^2
        den1<-r^a
    }
    du2db<-(-sum(num1*num2/den1))
    return(du2db)
}

##' @title Finds the first derivative of v with respect to a. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param nbits number of bits for mpfr function
##' @param r the rank vector
##' @return dvda
getdvda<-function(coefficients, r, nbits)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    if(nbits>0)
    {
        num1<-Rmpfr::mpfr(N+1-r,nbits)^Rmpfr::mpfr(b,nbits)
        num2<-log(r)
        den1<-Rmpfr::mpfr(r,nbits)^Rmpfr::mpfr(a,nbits)
    }
    if(nbits<=0)
    {
        num1<-(N+1-r)^b
        num2<-log(r)
        den1<-r^a
    }
    dvda<-(-sum(num1*num2/den1))
    return(dvda)
}

##' @title Finds the first derivative of v with respect to b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param coefficients the optimal values of a and b
##' @param nbits number of bits for mpfr function
##' @param r the rank vector
##' @return dvdb
getdvdb<-function( coefficients, r, nbits)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    if(nbits>0)
    {
        num1<-Rmpfr::mpfr(N+1-r,nbits)^Rmpfr::mpfr(b,nbits)
        num2<-log(N+1-r)
        den1<-Rmpfr::mpfr(r,nbits)^Rmpfr::mpfr(a,nbits)
    }
    if(nbits<=0)
    {
        num1<-(N+1-r)^b
        num2<-log(N+1-r)
        den1<-r^a
    }
    dvdb<-sum(num1*num2/den1)
    return(dvdb)
}

##' @title Finds the double derivative of A with with respect to a. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param u1 u1
##' @param v v
##' @param du1da First derivative of u1 with respect to parameter a
##' @param dvda First derivative of v with respect to parameter a
##' @return d2logAda2
getd2logAda2<-function(u1, v, du1da, dvda)
{
    num1<-v*du1da
    num2<-u1*dvda
    den1<-v^2
    d2logAda2<-(num1-num2)/den1
    return(d2logAda2)
}

##' @title Finds the double derivative of A with with respect to a and b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param u2 u2
##' @param v v
##' @param du2da First derivative of u2 with respect to a
##' @param dvda First derivative of v with respect to a
##' @return d2logAdadb
getd2logAdadb<-function( u2, v, du2da, dvda)
{
    num1<-v*du2da
    num2<-u2*dvda
    den1<-v^2
    d2logAdadb<-(num1-num2)/den1
    return(d2logAdadb)
}

##' @title Finds the double derivative of A with with respect to b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param u2 u2
##' @param v v
##' @param du2db First derivative of u2 with respect to b
##' @param dvdb First derivative of v with respect to
##' @return d2logAdb2
getd2logAdb2<-function( u2, v, du2db, dvdb)
{
    num1<-v*du2db
    num2<-u2*dvdb
    den1<-v^2
    d2logAdb2<-(num1-num2)/den1
    return(d2logAdb2)
}

##' @title Finds the double derivative of A with with respect to a and b. This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating the first and second order derivatives of A with respect to parameters a and b
##' @param u1 u1
##' @param v v
##' @param du1db First derivative of u1 with respect to b
##' @param dvdb First derivative of v with respect to b
##' @return d2logAdbda
getd2logAdbda<-function( u1, v, du1db, dvdb)
{
    num1<-v*du1db
    num2<-u1*dvdb
    den1<-v^2
    d2logAdbda<-(num1-num2)/den1
    return(d2logAdbda)
}
