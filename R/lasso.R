### X is a matrix; can be given as a matrix or as a matrix of type dgCMatrix from the Matrix package
### y is a vector; can not be given in sparse format
### wObs is a vector for the weights of the observations
### wLambda1 are weights associated with lambda1
### betaStart - starting position for the beta vector
### graph; following options are possible
###     NULL; then a straight line with weights 1 is assumed
###     vector of length 2 - giving the 2 dimensions of the graph; assumes equal weights
###     list with components - connections and weights; used for arbitrary graphs
### maxIterInner - maximum number of iterations until stop in inner iteration
### maxIterOuter - maximum number of iterations in outer loop
### accuracy - up to what error should the algorihtm run
### maxActivateVars - how many variables should be activated at the same time
### lambda1, lambda2 - vectors of the same length giving the values of lambda1 and lambda2; both have to be positive

lasso <- function(X, y, lambda1, wObs=NULL, wLambda1=NULL, betaStart=NULL, maxIterInner=10000, maxActivateVars=10, accuracy=1e-6) {
    ### check the variables that they have the correct format
    ### matrix should be sparse of class dgCMatrix
    if(is.matrix(X)) {
        X = as(X, "dgCMatrix")
    }
    
    ### check if y has the correct size
    if(length(y) != dim(X)[1]) {
        stop("y has incorrect length")
    }

    ### check that lambda1 and lambda2 have the right size and are all nonzero
    if(sum(lambda1 <= 0) > 0) {
        warning("Lambda1 has to have only positive elements; increasing to 1e-3")
        lambda1[lambda1 <= 0] = 1e-3
    }

    if(is.null(wObs)) {
        wObs = rep(1, dim(X)[1])
    }

    if(is.null(wLambda1)) {
        wLambda1 = rep(1, dim(X)[2])
    }
    else {
        stop("weights for lambda1 are currently not implemented")
    }

    if(is.null(betaStart)) {
        betaStart = rep(0, dim(X)[2])
    }

    ### adjust lambda for the number of observations
    n <- dim(X)[1]
    lambda1 <- n * lambda1

    maxIterInner = as.integer(maxIterInner)
    maxActivateVars = as.integer(maxActivateVars)

    res <- .Call("LassoWrapper", X, y, wObs, betaStart, wLambda1, maxIterInner, accuracy, maxActivateVars, lambda1)

    return(res)
}




