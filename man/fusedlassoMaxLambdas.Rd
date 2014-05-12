\name{fusedlassoMaxLambdas}
\Rdversion{1.1}
\alias{fusedlassoMaxLambdas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for fitting the fused lasso model
}
\description{
This function fits the standard version of the fused lasso. It can take a general matrix X and
    provides for possible weights on the \code{lambda1} and \code{lambda2} penalties. 
}
\usage{
fusedlassoMaxLambdas(X, y, family=c("gaussian", "binomial", "cox"),
wObs=NULL, wLambda1=NULL, graph=NULL, maxIterInner=10000,
maxIterOuter=100, maxActivateVars=10, maxNonZero = 2*length(y),
accuracy=1e-6, addIntercept=TRUE, Xunpenalized=NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
      The predictor matrix. Can be a regular matrix or of type \code{dgCMatrix}.
}
  \item{y}{
The respone vector.
}
\item{family}{
    Type of response, gaussian or binomial
}
  \item{wObs}{
Is a vector with weights for each observation. Defaults to 1 if \code{NULL}.
}
  \item{wLambda1}{
A vector with weights for \code{lambda1}. Currently not implemented (to be added soon).
}
  \item{graph}{
Description of the graph that corresponds to the \code{lambda2} penalty structure. If \code{NULL} is given a straight line is assumed. If a vector of length 2 is given, a 2-dimensional graph with equal weights is assumed, with dimensions as given in the vector. Otherwise a connection object has to be provided (add manual for this later).
}
  \item{maxIterInner}{
Maximum number of iterations in the inner loop to run.
}
  \item{maxIterOuter}{
Maximum number of iterations in the outer loop to run.
}
  \item{maxActivateVars}{
Maximum number of previously inactive variables to activate at the same time
}
\item{maxNonZero}{
Maximum number of beta that are ever non-zero
}
  \item{accuracy}{
Accuracy at which the algorithm will stop.
}
\item{addIntercept}{
    Should the model use an intercept
}
\item{Xunpenalized}{
Matrix of predictors that should not be penalized
}
}
\details{
Find the value for lambda1 and lambda2 at which the first penalized variable enters the models or
    any fused group breaks apart.
}
\value{
A list with the elements
\item{maxLambda1}{A sparse matrix of type \code{dgCMatrix} that returns the solutions for each value of
\code{maxLambda2} and \code{lambda2}. Each column is a separate solution.}
}
\author{
Holger Hoefling
}

\seealso{
\code{\link{testData}}
}

    
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models }
\keyword{regression}% __ONLY ONE__ keyword per line
\keyword{nonlinear}
\keyword{survival}
\keyword{multivariate}

