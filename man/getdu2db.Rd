% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROSeq.R
\name{getdu2db}
\alias{getdu2db}
\title{Finds the first derivative of u2 with respect to b. 
This first derivative is evaluated at the optimal (a_hat, b_hat).}
\usage{
getdu2db(coefficients, r)
}
\arguments{
\item{coefficients}{the optimal values of a and b}

\item{r}{the rank vector}
}
\value{
du2db
}
\description{
u1, v and u2 constitute the equations required for evaluating 
the first and second order derivatives of A with respect to parameters 
a and b
}
