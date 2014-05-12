#ifndef _QUADRATICDERIVATIVELOGISTIC_H_
#define _QUADRATICDERIVATIVELOGISTIC_H_

#include <iostream>
#include <vector>
#include "SparseMatrix.h"
#include "quadraticDerivative.h"

using namespace std;

// the class that keeps the data for calculating the derivative
// intended for logistic regression
// only differs from the standard class 
class QuadraticDerivativeLogistic : public QuadraticDerivative {
    public:
    int n; // number of rows in X
    int p; // number of cols in X

    SparseMatrix X; // copy of vector X
    vector<double> y; // copy of vector y
    vector<double> z; // adjusted response
    vector<double> w; // copy of vector beta
    vector<double> beta; // copy of R object that stores the current value of beta
    SparseMatrix WX; // stores a sparse matrix of WX
    vector<double> Xbeta; // stores the current value of Xbeta
    vector<double> diagXTWX;
    vector<double> XTWz;

    double cutoff;
    bool probIsNaN;
public:
	/*
	    Returns the value of the quadratic derivative
	    Encapsulates all necessary pre-computation in the Constructor
        Assumes that the dimensions of the objects are correct
	*/
    QuadraticDerivativeLogistic() {};

    // Constructor
    QuadraticDerivativeLogistic(const vector<double>& X, const vector<double>& y, const vector<double>& w, const vector<double>& beta);
    QuadraticDerivativeLogistic(const SparseMatrix& X, const vector<double>& y, const vector<double>& w, const vector<double>& beta);
   
    // make the appropriate calculations in the constructor after the data has been set
    void constructorWorker();

    ~QuadraticDerivativeLogistic() {};

    /*
        calculate the derivative for the position pos
    */
    double getDerivative(int pos);	
    
    inline double getProb(const int pos) const;

    void updateBeta(int pos, double newBeta);	

    inline double getBeta(int pos) { return beta[pos]; };

    inline vector<double> getBetaVec() { return beta; };

    inline int getBetaSize() { return beta.size(); };

    inline double getHessian(int pos) { return diagXTWX[pos]; }

    inline int getn() const { return X.n; }

    inline int getp() const { return X.p; }

    inline void activate(int pos) {};

    vector<double> getDerivativeVec();

    bool isExtreme();

};

#endif // _QUADRATICDERIVATIVELOGISTIC_H_

