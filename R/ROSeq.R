##' @title Modeling expression ranks for noise-tolerant differential 
##' expression analysis of scRNA-Seq data
##' @description Takes in the complete filtered and normalized read count
##' matrix, the location of the two sub-populations and the number of cores
##' to be used
##' @param countData The normalised and filtered, read count matrix, with
##' row names as genes name/ID and column names as sample id/name
##' @param condition Labels for the two sub-populations
##' @param numCores The number of cores to be used
##' @return pValues and FDR adjusted p significance values
##' @examples 
##' countData<-list()
##' countData$count<-ROSeq::L_Tung_single$NA19098_NA19101_count
##' countData$group<-ROSeq::L_Tung_single$NA19098_NA19101_group
##' head(countData$count)
##' gene_names<-rownames(countData$count)
##' countData$count<-apply(countData$count,2,function(x) as.numeric(x))
##' rownames(countData$count)<-gene_names
##' countData$count<-countData$count[,colSums(countData$count> 0) > 2000]
##' g_keep <- apply(countData$count,1,function(x) sum(x>2)>=3)
##' countData$count<-countData$count[g_keep,]
##' countData$count<-limma::voom(ROSeq::TMMnormalization(countData$count))
##' output<-ROSeq(countData=countData$count$E, condition = countData$group)
##' output
##' @export ROSeq
ROSeq<-function(countData, condition, numCores = 1)
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
    results <- pbmcapply::pbmclapply(geneIndex, initiateAnalysis, 
        scdata=countData,scgroups=scgroups, classOne=cOne, classTwo=cTwo, 
            mc.cores=numCores)
    pVals<-unlist(lapply(results,function(x) x[[12]]))
    pAdj<-stats::p.adjust(pVals, method = "fdr")
    results<-cbind("pVals"=pVals,"pAdj"=pAdj)
    rownames(results)<-rownames(countData)
    return(results)
}

##' @title Computes differential analysis for a given gene
##'
##' @description Takes in the row index which corresponds to a gene and 
##' evaluates for differential expression across two cell types.
##' @param gene The row index of the normalised and filtered, 
##' read count matrix
##' @param scdata The normalised and filtered, read count matrix
##' @param scgroups The location of the two sub-populations
##' @param classOne The location of the first sub-population, for example, 
##' sample names as given as columns names
##' @param classTwo The location of thesecond sub-population, for example, 
##' sample names as given as columns names
##' @return combinedResults A vector containing 12 values (gr1: a, g1: b, 
##' gr1: A, gr1: number of bins, gr1: R2, gr2: a, gr2: b, gr2: A, 
##' gr2: number of bins, gr2: R2, T, p)
initiateAnalysis<-function(gene, scdata, scgroups, classOne, classTwo)
{
    sp<-scdata[gene, ]
    spOne<-scdata[gene, which(scgroups%in%classOne)]
    spTwo<-scdata[gene, which(scgroups%in%classTwo)]
    geneStats<-getDataStatistics(sp, spOne, spTwo)
    results_groupOne <-findParams(spOne, geneStats)
    results_groupTwo <-findParams(spTwo, geneStats)
    T<-tryCatch(
    {
        computeDEG(results_groupOne, results_groupTwo)
    }, 
    warning=function(w) {NA}, error=function(esp) {NA})
    pValues<-stats::pchisq(T, df=2, lower.tail=FALSE)
    combinedResults<-c(results_groupOne[1], results_groupOne[2], 
    results_groupOne[3], results_groupOne[4], results_groupOne[5], 
    results_groupTwo[1], results_groupTwo[2], results_groupTwo[3], 
    results_groupTwo[4], results_groupTwo[5], T, pValues)
    return(combinedResults)
}

##' @title Evaluates statistics of the read counts corresponding to the gene
##'
##' @description Takes in the complete read count vector corresponding to the
##' gene (sp) and also the data corresponding to the two sub-populations 
##' (sp1 and sp2)
##'
##' @param sp The complete (normalized and filtered) read count data 
##' corresponding to the gene in question
##' @param spOne The (normalized and filtered) read count data corresponding 
##' to the first sub-population
##' @param spTwo The (normalized and filtered) read count data corresponding 
##' to the second sub-population
##' @return geneStats A vector containing 6 values corresponding to the gene 
##' data(maximum, minimum, mean, standard deviation, upper multiple of standard
##' deviation and lower multiple of standard deviation)
getDataStatistics<-function(sp, spOne, spTwo)
{
    maxds<-max(sp)
    minds<-min(sp)
    meands<-mean(sp)
    stdds<-stats::sd(sp)
    if(minds==maxds)
    {
        stop("Please remove genes with constant expression across the cells.")
    }
    ceilds<-ceiling((maxds-meands)/stdds)
    floords<-floor((minds-meands)/stdds)
    geneStats<-c(maxds, minds, meands, stdds, ceilds, floords)
    return (geneStats)
}

##' @title Finds the optimal values of parameters a and b that model 
##' the probability distribution of ranks, by Maximising the Log-Likelihood
##'
##' @description Takes in as input the read count data corresponding 
##' to one sub-population and the typical gene statistics. 
##' Then it splits the entire range into equally sized bins of 
##' size \eqn{k * \sigma}, where k is a scalar with a default value 
##' of 0.05, and \eqn{\sigma} is the standard deviation of the pulled 
##' expression estimates across the cell-groups. 
##' Each of these bins corresponds to a rank. Therefore, for each group, cell
##' frequency for each bin maps to a rank. These frequencies are normalized 
##' group-wise by dividing by the total cell count within a concerned group.
##' @param ds The (normalized and filtered) read count data corresponding to 
##' a sub-population
##' @param geneStats A vector containing 7 values corresponding to the gene 
##' data (maximum, minimum, mean, standard deviation, upper multiple of the
##' standard deviation, lower multiple of standard deviation and 
##' log_{2}(fold change))
##' @return results A vector containing 5 values (a, b, A, number of bins, R2)
findParams<-function(ds, geneStats)
{
    meands<-geneStats[3]
    stdds<-geneStats[4]
    ceilds<-geneStats[5]
    floords<-geneStats[6]
    step<-.05
    binNumber<-length(seq(floords, ceilds-step, step))
    opt_limit<-as.integer(log(1.7e300,binNumber)/2)
    rs<-vapply(seq(floords, ceilds-step, step),function(i)
    {
        LL<- meands+i*stdds
        UL <-meands+(i+step)*stdds
        length(intersect(which(ds<UL), which(ds>=LL)))
    }, numeric(1))
    fds<-rs
    number_of_bins<-length(fds)
    rank<-seq_len(number_of_bins)
    read_count_sorted<-sort(fds, decreasing=TRUE)
    normalized_read_count_sorted<-read_count_sorted/sum(read_count_sorted)
    model<-stats::optim(par = c(0.25, 3), minimizeNLL, r=rank, 
        readCount=normalized_read_count_sorted, method = "L-BFGS-B", 
            lower=c(-opt_limit,-opt_limit), upper=c(opt_limit,opt_limit))
    a<-model$par[1]
    b<-model$par[2]
    A<-1/sum((number_of_bins+1-rank)^b/(rank^a))
    f<-A*((number_of_bins+1-rank)^b)/(rank^a)
    SS_res<-sum((normalized_read_count_sorted-f)^2)
    SS_tot<-sum((normalized_read_count_sorted-
    mean(normalized_read_count_sorted))^2)
    R2<-as.numeric(1-SS_res/SS_tot)
    A<-as.numeric(A)
    results<-c(a, b, A, number_of_bins, R2)
    return(results)
}

##' @title Minimizes the Negative Log-Likelihood by iterating across 
##' values of parameters a and b
##'
##' @description Takes in as input a vector of values (coefficients), 
##' the number of bins and the normalized probability dsitribution of ranks
##' @param coefficients A vector containing two values for a and b
##' @param r The number of bins
##' @param readCount A vector of (normalized) frequency of read counts that 
##' fall within each bin
##' @return NLL Negative-Log Likelihood for the input coefficients
##' @seealso \code{\link{findParams}}
minimizeNLL<-function(coefficients, r, readCount)
{
    a<-coefficients[1]
    b<-coefficients[2]
    N<-length(r)
    sumReadCount<-sum(readCount)
    A <-1/sum(((N+1-r)^b)/(r^a))
    NLL<-a*sum(readCount*log(r)) - b*sum(readCount*log(N+1-r)) - 
        sumReadCount*log(A)
    NLL<-as.numeric(NLL)
    return (NLL)
}

##' @title Computes differential expression for the gene in question, 
##' by comparing the optimal parameters for sub-populations one and two
##' @description  Uses the (asymptotically) optimum two-sample Wald test  
##' based on the MLE of the parameters and its asymptotic variances given 
##' by the inverse of the Fisher information matrix
##' @param results_1 A vector corresponding to sub-population one and 
##' containing 5 values (a, b, A, number of bins, R2)
##' @param results_2 A vector corresponding to sub-population two and 
##' containing 5 values (a, b, A, number of bins, R2)
##' @return T  The Wald test statistic for testing the null hypothesis
##' @seealso \code{\link{getI}}, \code{\link{findParams}}
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
    T<-(m*n/(m+n))* t(matrix(c(a1-a2,b1-b2), nrow=2, ncol=1)) %*% 
        solve(w *V1 + (1-w)* V2, tol = 1e-20) %*% matrix(c(a1-a2,b1-b2),
            nrow=2, ncol=1)
    return(T)
}

##' @title Computes the Fisher Information Matrix
##' @description The Fisher Information Matrix and its derivatives are 
##' essential to evulate the minima of log likelihood
##' @param results A vector containing 5 values (a, b, A, number of bins, R2)
##' @return I  The Fisher Information Matrix
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
    d2logAda2<-getd2logAda2( u1, v=v, du1da=du1da, dvda=dvda)
    d2logAdb2<-getd2logAdb2( u1=u2, v=v, du1da=du2db, dvda=dvdb)
    d2logAdbda<-getd2logAdbda( u1=u1, v=v, du1da=du1db, dvda=dvdb)
    d2logAdadb<-getd2logAdadb( u1=u2, v=v, du1da=du2da, dvda=dvda)
    I<-c(-d2logAda2, -d2logAdadb, -d2logAdbda, -d2logAdb2)
    return(I)
}

##' @title Computes u1
##' @description u1, v and u2 constitute the equations required for evaluating 
##' the first and second order derivatives of A with respect to parameters 
##' a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return u1
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

##' @title Computes v
##' @description u1, v and u2 constitute the equations required for evaluating 
##' the first and second order derivatives of A with respect to parameters 
##' a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return v
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

##' @title Computes u2
##' @description u1, v and u2 constitute the equations required for evaluating 
##' the first and second order derivatives of A with respect to parameters 
##' a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return u2
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

##' @title Finds the first derivative of u1 with respect to a. 
##' This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating 
##' the first and second order derivatives of A with respect to parameters 
##' a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return du1da
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

##' @title Finds the first derivative of u1 with respect to b. 
##' This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating 
##' the first and second order derivatives of A with respect to parameters 
##' a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return du1db
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

##' @title Finds the first derivative of u2 with respect to a. 
##' This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating 
##' the first and second order derivatives of A with respect to parameters 
##' a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return du2da
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

##' @title Finds the first derivative of u2 with respect to b. 
##' This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating 
##' the first and second order derivatives of A with respect to parameters 
##' a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return du2db
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

##' @title Finds the first derivative of v with respect to a. 
##' This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating 
##' the first and second order derivatives of A with respect to parameters 
##' a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return dvda
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

##' @title Finds the first derivative of v with respect to b. 
##' This first derivative is evaluated at the optimal (a_hat, b_hat).
##' @description u1, v and u2 constitute the equations required for evaluating
##' the first and second order derivatives of A with respect to parameters 
##' a and b
##' @param coefficients the optimal values of a and b
##' @param r the rank vector
##' @return dvdb
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

##' @title Finds the double derivative of A
##' @description Finds the double derivative of A with with respect to 
##' a, (a, b), b , (a, b) in respective templates from right to left. 
##' This first derivative is evaluated at the optimal (a_hat, b_hat).
##' u1, v and u2 constitute the equations required for evaluating 
##' the first and second order derivatives of A with respect to parameters 
##' a and b
##' @param u1 u1
##' @param v v
##' @param du1da First derivative of u1 with respect to parameter a
##' @param dvda First derivative of v with respect to parameter a
##' @return d2logAda2
getd<-function(u1, v, du1da, dvda)
{
    num1<-v*du1da
    num2<-u1*dvda
    den1<-v^2
    d2logAda2<-(num1-num2)/den1
    return(d2logAda2)
}
getd2logAdbda<-getd2logAdb2<-getd2logAdadb<-getd2logAda2<-getd2logAda2<-getd

##' @title TMM Normalization.
##' @description Trimmed Means of M values (TMM) normalization 
##' (on the basis of edgeR package)
##' @param countTable The filtered, read count matrix, with row names
##' as genes name/ID and column names as sample id/name
##' @return countTableTMM
##' @examples 
##' countData<-list()
##' countData$count<-ROSeq::L_Tung_single$NA19098_NA19101_count
##' countData$group<-ROSeq::L_Tung_single$NA19098_NA19101_group
##' head(countData$count)
##' gene_names<-rownames(countData$count)
##' countData$count<-apply(countData$count,2,function(x) as.numeric(x))
##' rownames(countData$count)<-gene_names
##' countData$count<-countData$count[,colSums(countData$count> 0) > 2000]
##' g_keep <- apply(countData$count,1,function(x) sum(x>2)>=3)
##' countData$count<-countData$count[g_keep,]
##' countTableTMM<-ROSeq::TMMnormalization(countData$count)
##' countTableTMM
##' @export TMMnormalization
TMMnormalization <- function(countTable)
{
    cname<-colnames(countTable)
    rname<-rownames(countTable)
    nf<-edgeR::calcNormFactors(countTable ,method= "TMM")
    nf<- colSums(countTable)*nf
    scalingFactors <- nf/mean(nf)
    countTableTMM <- t(t(countTable)/scalingFactors)
    colnames(countTableTMM)<-cname
    rownames(countTableTMM)<-rname
    return(countTableTMM)
}
