#include <iostream>
#include <algorithm>
#include "quadraticDerivativeDiagonal.h"
#include "GeneralFunctions.h"

QuadraticDerivativeDiagonal::QuadraticDerivativeDiagonal(const SparseMatrix& X, const vector<double>& y, const vector<double>&  w, const vector<double>& beta) {
    // set number of rows and cols and store X, y and w and beta
    n = y.size();
    p = beta.size();
  
    this->X = X;
    this->y = y;
    this->w = w;
    this->beta = beta;

    constructorWorker();
}

QuadraticDerivativeDiagonal::QuadraticDerivativeDiagonal(const vector<double>& X, const vector<double>& y, const vector<double>& w, const vector<double>& beta) {
    // set number of rows and cols and store X, y and w and beta
    n = y.size();
    p = beta.size();
  
    this->X = SparseMatrix(X, n, p);
    this->y = y;
    this->w = w;
    this->beta = beta;

    constructorWorker();
}


void QuadraticDerivativeDiagonal::constructorWorker() {
    // precalculate the current WXbeta
    vector<int>::iterator vecIt;
    // create WX
    WX = X;
    WX.multByDiag(w);
    WXbeta.clear();
    WXbeta.resize(n, 0);
    double curBeta;
    int counter;
    for(int pos = 0; pos < beta.size(); ++pos) {
        WX.addMultOfColumn(WXbeta, pos, beta[pos]);
    }
    
    // it is useful to pre-calculate the diagonal of XTX and XTWy
    diagXTWX.clear();
    diagXTWX.resize(p, 0);
    XTWy.clear();
    XTWy.resize(p, 0);
    for(int i = 0; i < p;  ++i) {
        diagXTWX[i] = X.innerProd(i, w);
        XTWy[i] = WX.innerProd(y,i);
    }
}


double QuadraticDerivativeDiagonal::getDerivative(int pos) {
    // calculate the rest of the expression XTWXbeta
    double res = X.innerProd(WXbeta, pos);

    // now subtract XTWy
    res -= XTWy[pos];
    return(res);
}


void QuadraticDerivativeDiagonal::updateBeta(int pos, double newBeta) {
    // update beta and save the difference
    double betaDiff = newBeta - beta[pos];
    if(betaDiff == 0) return;

    beta[pos] = newBeta;

    // now update Xbeta
    WX.addMultOfColumn(WXbeta, pos, betaDiff);
}

vector<double> QuadraticDerivativeDiagonal::getDerivativeVec() {
    vector<double> res(p,0);
    for(int i = 0; i < p; ++i) {
        res[i] = getDerivative(i);
    }
    return(res);
}

