### write a predict function
### assumes that lambda2 is constant (can be chosen exactly in the fusedlasso function
### anyway)
### also assumes that lambda1 is decreasing
### oldX and Surv only needed for baseline prediction
### could move baseline prediction into fusedlasso function as standard to eliminate this
predict.fusedlasso <- function(fl, newX, lambda1Vec, Surv=NULL, oldX=NULL) {
    ### first calculate the intermediate position of the betas
    maxLambda1 <- max(fl$lambda1)
    minLambda1 <- min(fl$lambda1)

    ### now check that the lambda1 are in the correct interval
    outOfBound <- lambda1Vec > maxLambda1 || lambda1Vec < minLambda1
    if(any(outOfBound)) {
        cat("Some lambda1 are out of bounds")
        print(lambda1Vec[outOfBound])
    }

    lambda1Vec <- lambda1Vec[!outOfBound]

    ### now calculate the interpolation 
    flLambda1 <- fl$lambda1
  
    counter <- 1 
    res <- matrix(numeric(length(lambda1Vec) * nrow(newX)), nrow=nrow(newX))
    counter <- 1
    for(lambda1 in lambda1Vec) {
        ### for each lambda, find a an interpolation
        inter <- interpolatePoints(lambda1, fl$lambda1)
        newBeta <- fl$beta[,inter$smallerPos]*inter$smallerWeight + fl$beta[,inter$largerPos]*inter$largerWeight
        newIntercept <- fl$Intercept[inter$smallerPos]*inter$smallerWeight + fl$Intercept[inter$largerPos]*inter$largerWeight
        if(fl$family=="gaussian") {
            res[,counter] <- predictGaussian(newBeta, newIntercept, newX)
        }
        else if(fl$family=="binomial") {
            res[,counter] <- predictBinomial(newBeta, newIntercept, newX)
        }
        else if(fl$family=="Cox") {
            res[,counter] <- predictCox(newBeta, newX, Surv, oldX)
        }
        else {
            stop("Unknown distribution!")
        }
        res[,counter] <-  
        counter <- counter + 1
    }     
    
    return(list(prediction=res, lambda1=lambda1Vec))
}


predictBinomial <- function(beta, intercept, newX) {
    newEff <- intercept + newX %*% beta
    newp <- 1/(1+exp(-newEff))
    return(newp)
}

predictGaussian <- function(beta, intercept, newX) {
    return(newX %*% beta + intercept)
}

predictCox <- function(beta, newX, Surv, oldX) {
    stop("This feature has not yet been implemented")
}

interpolatePoints <-function(x, xVec) {
    foo <- nextLarger(x, xVec)
    largerPos <- foo$pos
    largerVal <- foo$val
    foo <- nextSmaller(x, xVec)
    smallerPos <- foo$pos
    smallerVal <- foo$val
    largerWeight <- 1
    smallerWeight <- 0

    if(smallerPos != largerPos) {
        largerWeight <- (x-smallerVal)/(largerVal-smallerVal)
        smallerWeight <- 1-largerWeight
    }

    return(list(smallerPos=smallerPos, largerPos=largerPos, smallerWeight=smallerWeight, largerWeight=largerWeight))
}

nextLarger <- function(x, xVec) {
    foo <- xVec -x
    foo[foo < 0] <- 1e10
    pos <- which.min(foo)
    val <- xVec[pos]
    return(list(pos=pos, val=val))
}

nextSmaller <- function(x, xVec) {
    foo <- xVec - x
    foo[foo > 0] <- -1e10
    pos <- which.max(foo)
    val <- xVec[pos]
    return(list(pos=pos, val=val))
}

