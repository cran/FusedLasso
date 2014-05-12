#include <iostream>
#include <algorithm>
#include "quadraticDerivativeLogistic.h"
#include "GeneralFunctions.h"
#include <math.h>
#include <stdlib.h>

QuadraticDerivativeLogistic::QuadraticDerivativeLogistic(const SparseMatrix& X, const vector<double>& y, const vector<double>&  w, const vector<double>& beta) {
    // set number of rows and cols and store X, y and w and beta
    n = y.size();
    p = beta.size();
  
    this->X = X;
    this->y = y;
    this->w = w;
    this->beta = beta;

    constructorWorker();
}

QuadraticDerivativeLogistic::QuadraticDerivativeLogistic(const vector<double>& X, const vector<double>& y, const vector<double>& w, const vector<double>& beta) {
    // set number of rows and cols and store X, y and w and beta
    n = y.size();
    p = beta.size();
  
    this->X = SparseMatrix(X, n, p);
    this->y = y;
    this->w = w;
    this->beta = beta;

    constructorWorker();
}


void QuadraticDerivativeLogistic::constructorWorker() {
    // precalculate the current WXbeta
    probIsNaN = false;
    vector<int>::iterator vecIt;
    Xbeta.clear();
    Xbeta.resize(n, 0);
    double curBeta;
    int counter;
    for(int pos = 0; pos < beta.size(); ++pos) {
        X.addMultOfColumn(Xbeta, pos, beta[pos]);
    }
    
    z.resize(n, 0);
    // adjust the weights and the response
    for(int i = 0; i < n; ++i) {
        double prob = getProb(i);
        w[i] *= prob * (1 - prob);
        z[i] = Xbeta[i] + ((y[i] - prob) / (prob * (1 - prob))); 
        if(isnan(z[i])) {
            probIsNaN = true;
        }
    }

    // create WX
    WX = X;
    WX.multByDiag(w);

    // it is useful to pre-calculate the diagonal of XTX and XTWy
    diagXTWX.clear();
    diagXTWX.resize(p, 0);
    XTWz.clear();
    XTWz.resize(p, 0);
    for(int i = 0; i < p;  ++i) {
        diagXTWX[i] = X.innerProd(i, w);
        XTWz[i] = WX.innerProd(z,i);
    }

    cutoff = 0.01;
}


double QuadraticDerivativeLogistic::getDerivative(int pos) {
    // calculate the rest of the expression XTWXbeta
    double res = WX.innerProd(Xbeta, pos);
    if(isnan(res)) {
      error("NAN1\n");
    }
    // now subtract XTWy
    res -= XTWz[pos];
    if(isnan(res)) {
      error("NAN2\n");
    }

    return(res);
}

inline double QuadraticDerivativeLogistic::getProb(const int pos) const {
    return(1/(1+exp(-Xbeta[pos])));
}

bool QuadraticDerivativeLogistic::isExtreme() {
    if(probIsNaN) {
        return(true);
    }
    for(int i = 0; i < n; ++i) {
        if(getProb(i) > cutoff || getProb(i) < 1 - cutoff) {
            return(false);
        }
    }
    return(true);
}


void QuadraticDerivativeLogistic::updateBeta(int pos, double newBeta) {
    // update beta and save the difference
    double betaDiff = newBeta - beta[pos];
    if(betaDiff == 0) return;

    beta[pos] = newBeta;

    // now update Xbeta
    X.addMultOfColumn(Xbeta, pos, betaDiff);
}

vector<double> QuadraticDerivativeLogistic::getDerivativeVec() {
    vector<double> res(p,0);
    for(int i = 0; i < p; ++i) {
        res[i] = getDerivative(i);
    }
    return(res);
}

