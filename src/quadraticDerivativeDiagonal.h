#ifndef _QUADRATICDERIVATIVEDIAGONAL_H_
#define _QUADRATICDERIVATIVEDIAGONAL_H_

#include <iostream>
#include <vector>
#include "SparseMatrix.h"
#include "quadraticDerivative.h"

using namespace std;

// the class that keeps the data for calculating the derivative
class QuadraticDerivativeDiagonal : public QuadraticDerivative {
    public:
    int n; // number of rows in X
    int p; // number of cols in X

    SparseMatrix X; // copy of vector X
    vector<double> y; // copy of vector y
    vector<double> w; // copy of vector beta
    vector<double> beta; // copy of R object that stores the current value of beta
    SparseMatrix WX; // stores a sparse matrix of WX
    vector<double> WXbeta; // stores the current value of Xbeta
    vector<double> diagXTWX;
    vector<double> XTWy;

public:
	/*
	    Returns the value of the quadratic derivative
	    Encapsulates all necessary pre-computation in the Constructor
        Assumes that the dimensions of the objects are correct
	*/
    QuadraticDerivativeDiagonal() {};

    // Constructor
    QuadraticDerivativeDiagonal(const vector<double>& X, const vector<double>& y, const vector<double>& w, const vector<double>& beta);
    QuadraticDerivativeDiagonal(const SparseMatrix& X, const vector<double>& y, const vector<double>& w, const vector<double>& beta);
   
    // make the appropriate calculations in the constructor after the data has been set
    void constructorWorker();

    ~QuadraticDerivativeDiagonal() {};

    /*
        calculate the derivative for the position pos
    */
    double getDerivative(int pos);	

    void updateBeta(int pos, double newBeta);	

    inline double getBeta(int pos) { return beta[pos]; };

    inline vector<double> getBetaVec() { return beta; };

    inline int getBetaSize() { return beta.size(); };

    inline double getHessian(int pos) { return diagXTWX[pos]; }

    inline int getn() const { return X.n; }

    inline int getp() const { return X.p; }

    inline void activate(int pos) {};

    vector<double> getDerivativeVec();

    bool isExtreme() { return false; };
};

#endif // _QUADRATICDERIVATIVEDIAGONAL_H_

