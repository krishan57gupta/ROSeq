\name{L_Tung_single}
\alias{L_Tung_single}
\docType{data}
\title{
    Single cell samples for DE genes analysis
}
\description{
    Three replicates from three human induced pluripotent stem cell (iPSC)
    lines were considered.  96  single  cells  were considered  in  each  of
    the three replicates corresponding to one of the three individuals
    (these shall be referred to by their labels NA19098,NA19101 and NA19239)
}
\usage{data("L_Tung_single")}
\format{
    The format is: list of cells corresponding NA19098 versus NA19101
    and groups labels.
}
\details{
    filtered and normalized  data
}
\source{
    Tung, P.-Y.et al.Batch effects and the effective design of 
    single-cell geneexpression studies.Scientific reports7, 39921 (2017).
}
\references{
    Tung, P.-Y.et al.Batch effects and the effective design of 
    single-cell geneexpression studies.Scientific reports7, 39921 (2017).
}
\examples{
    data(L_Tung_single)
    ## summary(ROSeq::L_Tung_single)
}
\keyword{datasets}
