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

fusedlasso <- function(X, y, lambda2, numLambda1=100, minFracOfMaxLambda1=1e-3, lambda1=NULL, family=c("gaussian", "binomial"), wObs=NULL, wLambda1=NULL, graph=NULL, betaStart=NULL, interceptStart=NULL, unpenalizedStart=NULL, maxIterInner=10000, maxIterOuter=100, maxActivateVars=10, accuracy=1e-6, maxNonZero=2*nrow(X), addIntercept=TRUE, Xunpenalized=NULL, fusionCheck=c("all","active","none", "naive", "huber"), verbose=FALSE) {

    addIntercept <- eval(addIntercept)
    interceptStart <- eval(interceptStart)
    Xunpenalized <- eval(Xunpenalized)
    unpenalizedStart <- eval(unpenalizedStart)

    addNodes <- 0
    ### check the variables that they have the correct format
    checkInputVariables()
    verbose <- as.logical(verbose)   

    fusionCheck <- match.arg(fusionCheck, c("all", "active", "none", "naive","huber"))
    if(fusionCheck == "all") {
        maxFusionLevel <- as.integer(2)
    }
    else if(fusionCheck == "active") {
        maxFusionLevel <- as.integer(1)
    }
    else if(fusionCheck == "none") {
        maxFusionLevel <- as.integer(0)
    }
    else if(fusionCheck == "naive") {
        maxFusionLevel <- as.integer(-1)
    }
    else {
        maxFusionLevel <- as.integer(-2)
    }

    ### check the lambdas for consistency 
    if(is.null(lambda1)) {
        adjustWithMax <- TRUE
        lambda1 <- expGrid(numLambda1, minFracOfMaxLambda1)
    }
    else {
        adjustWithMax <- FALSE
        lambda1 <- sort(lambda1, decreasing=TRUE)
    }

    if(length(lambda2) != 1 && length(lambda1) != length(lambda2)) {
        stop("Give either 1 value for lambda2 or same number as lambda2")
    }

    if(lambda2[1] <= 0) {
        stop("Lambda2 must be positive")
    }


    if(sum(lambda2 <= 0) > 0) {
        stop("Lambda2 has to have only positive elements")
    }

    if(length(lambda2)==1) {
        lambda2 <- rep(lambda2[1], length(lambda1))
    }
    lambda1 <- nrow(X) * lambda1
    lambda2 <- nrow(X) * lambda2

    res <- .Call("FusedLassoWrapperL1", X, y, wObs, betaStart, wLambda1, graph, maxIterInner, maxIterOuter, accuracy, maxActivateVars, lambda1, lambda2, maxNonZero, addNodes, family, adjustWithMax, maxFusionLevel, verbose)

    ### if an intercept was added, present it seperately, or set it to 0 otherwise
    addedBetas <- matrix(numeric(0), ncol=1)
    if(addNodes > 0) {
        addedBetas <- res$beta[(nrow(res$beta)-addNodes+1):nrow(res$beta),,drop=FALSE]
        res$beta <- res$beta[1:(nrow(res$beta)-addNodes),,drop=FALSE]
    }
    if(addIntercept) {
        res$Intercept <- addedBetas[1,]
        addedBetas <- addedBetas[-1,,drop=FALSE]

    }
    else {
        res$Intercept <- rep(0, length(lambda1))
    }
    if(nrow(addedBetas) > 0) {
        res$unpenalized <- addedBetas
    }
    else {
        res$unpenalized <- NULL
    }

    res$lambda1 <- res$lambda1 / nrow(X)
    res$lambda2 <- res$lambda2 / nrow(X)

    res$maxLambda1 <- res$maxLambda1 / nrow(X)
    res$maxLambda2 <- res$maxLambda2 / nrow(X)

    class(res) <- "fusedlasso"

    return(res)
}

graphlassoL2 <- function(X, y, lambda1, lambda2, family=c("gaussian", "binomial"), wObs=NULL, wLambda1=NULL, graph=NULL, betaStart=NULL, interceptStart=NULL, unpenalizedStart=NULL, maxIterInner=10000, maxIterOuter=100, maxActivateVars=10, accuracy=1e-6, maxNonZero=2*nrow(X), addIntercept=TRUE, Xunpenalized=NULL, verbose=FALSE) {

    addIntercept <- eval(addIntercept)
    interceptStart <- eval(interceptStart)
    Xunpenalized <- eval(Xunpenalized)
    unpenalizedStart <- eval(unpenalizedStart)

    addNodes <- 0

    checkInputVariables()
    verbose <- as.logical(verbose)   

    if(length(lambda2) != 1 && length(lambda1) != 1 && length(lambda1) != length(lambda2)) {
        stop("Vectors lambda1 and lambda2 have to be of same length or length 1")
    }

    if(sum(lambda1 <= 0) > 0) {
        stop("Lambda1 has to have only positive elements")
    }

    if(sum(lambda2 <= 0) > 0) {
        stop("Lambda2 has to have only positive elements")
    }

    if(length(lambda1)==1) {
        lambda1 <- rep(lambda1[1], length(lambda2))
    }
    if(length(lambda2)==1) {
        lambda2 <- rep(lambda2[1], length(lambda1))
    }

    lambda1 <- nrow(X) * lambda1
    lambda2 <- nrow(X) * lambda2

    res <- .Call("FusedLassoWrapperL2", X, y, wObs, betaStart, wLambda1, graph, maxIterInner, maxIterOuter, accuracy, maxActivateVars, lambda1, lambda2, maxNonZero, addNodes, family, verbose)

    ### if an intercept was added, present it seperately, or set it to 0 otherwise
    addedBetas <- matrix(numeric(0), ncol=1)
    if(addNodes > 0) {
        addedBetas <- res$beta[(nrow(res$beta)-addNodes+1):nrow(res$beta),,drop=FALSE]
        res$beta <- res$beta[1:(nrow(res$beta)-addNodes),,drop=FALSE]
    }
    if(addIntercept) {
        res$Intercept <- addedBetas[1,]
        addedBetas <- addedBetas[-1,,drop=FALSE]

    }
    else {
        res$Intercept <- rep(0, length(lambda1))
    }
    if(nrow(addedBetas) > 0) {
        res$unpenalized <- addedBetas
    }
    else {
        res$unpenalized <- NULL
    }

    res$lambda1 <- res$lambda1 / nrow(X)
    res$lambda2 <- res$lambda2 / nrow(X)

    res$maxLambda1 <- res$maxLambda1 / nrow(X)
    res$maxLambda2 <- res$maxLambda2 / nrow(X)

    class(res) <- "fusedlasso"

    return(res)
}


fusedlassoMaxLambdas <- function(X, y, family = c("gaussian", "binomial", "cox"), wObs=NULL, wLambda1=NULL, graph=NULL, maxIterInner=10000, maxIterOuter=100, maxActivateVars=10, maxNonZero=2*length(y), accuracy=1e-6, addIntercept=TRUE, Xunpenalized=NULL) {
    ### check the variables that they have the correct format
    addIntercept <- eval(addIntercept)
    Xunpenalized <- eval(Xunpenalized)

    addNodes <- 0
    
    checkInputVariables()

    ### which variables are to be excluded
    exemptVars <- as.integer(which(wLambda1 <= 1e-4))

    res <- .Call("FusedLassoWrapperMaxLambdas", X, y, wObs, wLambda1, graph, maxIterInner, accuracy, maxActivateVars, exemptVars, addNodes, family)
    
    res$maxLambda1 <- res$maxLambda1 / nrow(X)
    res$maxLambda2 <- res$maxLambda2 / nrow(X)

    return(res)
}


expGrid <- function(numGrid, low=0.001) {
    mult <- log(low)/numGrid
    return(exp((1:numGrid) * mult))
}


checkInputVariables <- function() {
    env <- parent.frame()
    ### first, set default values for all variables where necessary
    if(evalq(is.null(wObs),env)) {
        evalq(wObs <- rep(1, nrow(X)),env)
    }
   
    if(evalq(is.null(wLambda1),env)) {
        evalq(wLambda1 <- rep(1, ncol(X)),env)
    }

    if(evalq(is.null(graph),env)) {
        evalq(graph <- as.integer(ncol(X)),env)
    }
    else if(evalq(is.vector(graph) && !is.list(graph),env)) {
        if(evalq(length(graph) == 1, env)) {
            evalq(graph <- as.integer(graph), env)
            if(evalq(graph[1] != ncol(X), env)) {
                stop("One-dimensional graph has to have same size as vector beta")
            }
        }
        else if(evalq(length(graph) == 2 || length(graph) == 3,env)) {
            evalq(graph <- as.integer(graph),env)
            if(evalq(prod(graph) != ncol(X),env)) {
                stop("2-dimensional graph has to have same size as vector beta")
            } 
        }
        else {
            stop("If a vector is given for graph, it has to have length 2")
        }
    }
    else if(evalq(is.list(graph),env)) {
        cat("Custom graph given\n")
    }

    if(evalq(!exists("betaStart") || is.null(betaStart),env)) {
        evalq(betaStart <- rep(0, ncol(X)),env)
    }
    if(evalq(!exists("interceptStart") || is.null(interceptStart), env)) {
        evalq(interceptStart <- 0.001, env)
    }

    ### matrix should be sparse of class dgCMatrix
    if(evalq(is.matrix(X),env)) {
        evalq(X <- as(X, "dgCMatrix"), env)
    }


    evalq(family <- match.arg(family, c("gaussian", "binomial")), env)
    ### now check that the input variables are ok
    if(evalq(family=="gaussian", env)) {
        ### input y has to be a numeric vector
        evalq(y <- as.numeric(y), env)
    }
    else if(evalq(family=="binomial", env)) {
        y <- get("y", envir=env)
        ### if it is a vector, it is assumed to be either 0 or 1
        if(is.vector(y)) {
            evalq(y <- as.numeric(y), env)
            if(any(y<0) || any(y>1)) {
                stop("y has to be between 0 and 1")
            }
        }
        else if(is.matrix(y)) {
            if(dim(y)[2] != 2) {
                stop("y has to be a 2-column matrix")
            }
            fracs <- as.numeric(y[,2] / (y[,1] + y[,2]))
            numTrial <- as.numeric(y[,1] + y[,2])
            assign("y", fracs, envir=env)
            wObs <- get("wObs", envir=env)
            wObs <- wObs * numTrial
            assign("wObs", wObs, envir=env)
        }
        else {
            stop("y has to be a vector or a 2-column matrix")
        }
    }
    else {
        stop("Not a valid family")
    }

    ### adjust lambda for the number of observations
    evalq(maxIterInner <- as.integer(maxIterInner), env)
    evalq(maxIterOuter <- as.integer(maxIterOuter), env)
    evalq(maxActivateVars <- as.integer(maxActivateVars), env)
    evalq(maxNonZero <- as.integer(maxNonZero), env)

    addIntercept <- get("addIntercept", env)
    if(addIntercept) {
        ### add an intercept as the last variable
        evalq(X <- cBind(X, rep(1,dim(X)[1])), env)
        evalq(betaStart <- c(betaStart, interceptStart[1]), env) # so it is in the model from the beginning
        evalq(wLambda1 <- c(wLambda1, 1e-6), env)
        evalq(addNodes <- 1, env)
    }
    Xunpenalized <- get("Xunpenalized", env)
    if(evalq(!is.null(Xunpenalized), env)) {
        if(evalq(is.null(unpenalizedStart), env)) {
            evalq(unpenalizedStart <- rep(0, ncol(Xunpenalized)), env)
        }
        evalq(X <- cBind(X, Xunpenalized), env)
        evalq(X <- as(X, "dgCMatrix"), env)
        evalq(betaStart <- c(betaStart, unpenalizedStart), env)
        evalq(wLambda1 <- c(wLambda1, rep(1e-6, ncol(Xunpenalized))), env)
        evalq(addNodes <- addNodes + ncol(Xunpenalized), env)
    }
    evalq(addNodes <- as.integer(addNodes), env)
}

