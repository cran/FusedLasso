\name{lasso}
\Rdversion{1.1}
\alias{lasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for fitting the lasso model
}
\description{
This function fits the standard version of the lasso. It can take a general matrix X and
    provides for possible weights on the \code{lambda1} penalty. 
}
\usage{
lasso(X, y, lambda1, wObs = NULL, wLambda1 = NULL, betaStart = NULL,
maxIterInner = 10000, maxActivateVars = 10, accuracy = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
      The predictor matrix. Can be a regular matrix or of type \code{dgCMatrix}.
}
  \item{y}{
The respone vector.
}
  \item{lambda1}{
A vector of values for \code{lambda1}. Has to be non-negative and same length as \code{lambda2}.
}
  \item{wObs}{
Is a vector with weights for each observation. Defaults to 1 if \code{NULL}.
}
  \item{wLambda1}{
A vector with weights for \code{lambda1}. Currently not implemented (to be added soon).
}
  \item{betaStart}{
Starting value for the coefficients in the algorithm. 
}
  \item{maxIterInner}{
Maximum number of iterations in the inner loop to run.
}
\item{maxActivateVars}{
Maximum number of previously inactive variables to activate at the same time
}
\item{accuracy}{
Accuracy at which the algorithm will stop.
}
}
\details{
Will have to edit some details here later.
}
\value{
A list with the elements
\item{beta}{A sparse matrix of type \code{dgCMatrix} that returns the solutions for each value of
\code{lambda1}. Each column is a separate solution.}
\item{lambda1}{The vector of the lambda1 values that were used.}
\item{converged}{A logical vector inidicating if the algorithm converged.}
}
\author{
Holger Hoefling
}

    
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models }
\keyword{regression}% __ONLY ONE__ keyword per line
\keyword{multivariate}

