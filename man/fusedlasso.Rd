\name{fusedlasso}
\Rdversion{1.1}
\alias{fusedlasso}
\alias{graphlassoL2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for fitting the fused lasso model
}
\description{
This function fits the standard version of the fused lasso. It can take a general matrix X and
    provides for possible weights on the \code{lambda1} and \code{lambda2} penalties. 
}
\usage{
fusedlasso(X, y, lambda2, numLambda1=100, minFracOfMaxLambda1=1e-3,
lambda1=NULL, family = c("gaussian", "binomial"), wObs = NULL,
wLambda1 = NULL, graph = NULL, betaStart = NULL, interceptStart=NULL,
unpenalizedStart=NULL, maxIterInner = 10000, maxIterOuter = 100,
maxActivateVars = 10, accuracy = 1e-06, maxNonZero=2*nrow(X),
addIntercept=TRUE, Xunpenalized=NULL,
fusionCheck=c("all", "active", "none","naive","huber"), verbose=FALSE)
graphlassoL2(X, y, lambda1, lambda2, family = c("gaussian", "binomial"),
wObs = NULL, wLambda1 = NULL, graph = NULL, betaStart = NULL,
interceptStart=NULL, unpenalizedStart=NULL, maxIterInner = 10000,
maxIterOuter = 100, maxActivateVars = 10, accuracy = 1e-06,
maxNonZero=2*nrow(X), addIntercept=TRUE, Xunpenalized=NULL, verbose=FALSE)
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
A vector of values for \code{lambda1}. Has to be non-negative. IF NULL, an automatic
    sequence of length \code{numLambda1} is being generated from a sufficiently large value expentionally decreasing to
    \code{minFracOfMaxLambda1} of that value.
}
\item{numLambda1}{length of the automatic, exponentional sequence}
\item{minFracOfMaxLambda1}{Down to which fraction of the maximum value of lambda1 should the
    automatic sequence go}
\item{lambda2}{
A single value for \code{lambda2} or a vector of the same length as \code{lambda1}. Has to be non-negative.
}
\item{family}{
    Should a gaussian or a binomial distribution be used for the regression.
}
  \item{wObs}{
Is a vector with weights for each observation. Defaults to 1 if \code{NULL}.
}
  \item{wLambda1}{
A vector with weights for \code{lambda1}.
}
  \item{graph}{
Description of the graph that corresponds to the \code{lambda2} penalty structure. If \code{NULL} is given a straight line is assumed. If a vector of length 2 is given, a 2-dimensional graph with equal weights is assumed, with dimensions as given in the vector. Otherwise a connection object has to be provided (add manual for this later).
}
  \item{betaStart}{
Starting value for the coefficients in the algorithm. 
}
\item{interceptStart}{
Starting value for the intercept.
}
\item{unpenalizedStart}{
Starting value for the unpenalized part of X that is given
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
  \item{accuracy}{
Accuracy at which the algorithm will stop.
}
\item{maxNonZero}{
How many coordinates of beta should at most be non-zero. Breaks off if there are too many.
}
\item{addIntercept}{
        Should the model use an intercept
    }
\item{Xunpenalized}{
Matrix of predictors that are unpenalized and added to the model
}
\item{fusionCheck}{
    Should the fused sets be checked for breaking up? 
}
\item{verbose}{Should the function give some output what it is doing?}
}
\details{
Will have to edit some details here later.
}
\value{
A list with the elements
\item{beta}{A sparse matrix of type \code{dgCMatrix} that returns the solutions for each value of
\code{lambda1} and \code{lambda2}. Each column is a separate solution.}
\item{lambda1}{The vector of the lambda1 values that were used.}
\item{lambda2}{The vector of the lambda2 values that were used.}
\item{success}{A logical vector inidicating if the algorithm converged.}
\item{Intercept}{The intercept of the model; 0 if no intercept was included}
}
\author{
Holger Hoefling
}

\seealso{
\code{\link{testData}}
}
\examples{
\dontrun{
data(testData)
testData$Exn10p20Res <- fusedlasso(testData$Exn10p20X,
testData$Exn10p20y,lambda1=testData$Exn10p20Lambda1,
lambda2=testData$Exn10p20Lambda2, addIntercept=FALSE)
testData$Exn10p20Correct <- (testData$Exn10p20Res$beta -testData$Exn10p20Sol) < 2e-4
testData$Exn10p20Correct[as.logical(is.na(testData$Exn10p20Correct))] <- TRUE
stopifnot(as.vector(testData$Exn10p20Correct))

testData$Exn20p1000Res <- fusedlasso(testData$Exn20p1000X, testData$Exn20p1000y,
lambda1=testData$Exn20p1000Lambda1,lambda2=testData$Exn20p1000Lambda2, addIntercept=FALSE)
testData$Exn20p1000Correct <- (testData$Exn20p1000Res$beta -testData$Exn20p1000Sol) < 2e-4
testData$Exn20p1000Correct[as.logical(is.na(testData$Exn20p1000Correct))] <- TRUE
stopifnot(as.vector(testData$Exn20p1000Correct))

testData$Exn20p100Res <- fusedlasso(testData$Exn20p100X, testData$Exn20p100y,
lambda1=testData$Exn20p100Lambda1, lambda2=testData$Exn20p100Lambda2, addIntercept=FALSE)
testData$Exn20p100Correct <- (testData$Exn20p100Res$beta -testData$Exn20p100Sol) < 2e-4
testData$Exn20p100Correct[as.logical(is.na(testData$Exn20p100Correct))] <- TRUE
stopifnot(as.vector(testData$Exn20p100Correct))

testData$Exn20p100Res2Dim <- fusedlasso(testData$Exn20p100X, testData$Exn20p100y,
lambda1=testData$Exn20p100Lambda1, lambda2=testData$Exn20p100Lambda2,
graph=c(10,10), addIntercept=FALSE)
testData$Exn20p100Correct2Dim <- (testData$Exn20p100Res2Dim$beta -testData$Exn20p100Sol2Dim) < 2e-4
testData$Exn20p100Correct2Dim[as.logical(is.na(testData$Exn20p100Correct2Dim))] <- TRUE
stopifnot(as.vector(testData$Exn20p100Correct2Dim))

testData$Exn50p20Res <- fusedlasso(testData$Exn50p20X, testData$Exn50p20y,
lambda1=testData$Exn50p20Lambda1, lambda2=testData$Exn50p20Lambda2, addIntercept=FALSE)
testData$Exn50p20Correct <- (testData$Exn50p20Res$beta -testData$Exn50p20Sol) < 2e-4
testData$Exn50p20Correct[as.logical(is.na(testData$Exn50p20Correct))] <- TRUE
stopifnot(as.vector(testData$Exn50p20Correct))

testData$Exn50p100LogRes <- fusedlasso(testData$Exn50p100LogX, testData$Exn50p100Logy,
lambda1=testData$Exn50p100LogLambda1, lambda2=testData$Exn50p100LogLambda2,
family="binomial", addIntercept=TRUE)
testData$Exn50p100LogCorrect <- (testData$Exn50p100LogRes$beta -testData$Exn50p100LogSol) < 2e-4
testData$Exn50p100LogCorrect[as.logical(is.na(testData$Exn50p100LogCorrect))] <- TRUE
stopifnot(as.vector(testData$Exn50p100LogCorrect))
stopifnot((testData$Exn50p100LogSolInt - testData$Exn50p100LogRes$Intercept) < 1e-3)

}
}

    
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models }
\keyword{regression}% __ONLY ONE__ keyword per line
\keyword{nonlinear}
\keyword{survival}
\keyword{multivariate}

