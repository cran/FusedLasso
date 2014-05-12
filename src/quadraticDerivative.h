#ifndef _QUADRATICDERIVATIVE_H_
#define _QUADRATICDERIVATIVE_H_

#include <iostream>
#include <vector>
#include "SparseMatrix.h"

using namespace std;

// the class that keeps the data for calculating the derivative
class QuadraticDerivative {

public:
    /*
        calculate the derivative for the position pos
    */
    virtual ~QuadraticDerivative() {};

    virtual double getDerivative(int pos) = 0;	

    virtual void updateBeta(int pos, double newBeta) = 0;	

    virtual double getBeta(int pos) = 0;

    virtual vector<double> getBetaVec() = 0;

    virtual int getBetaSize() = 0;

    virtual double getHessian(int pos) = 0;

    virtual int getn() const = 0;

    virtual int getp() const = 0;

    virtual void activate(int pos) = 0;

    virtual vector<double> getDerivativeVec() = 0;

    virtual bool isExtreme() {return(false);};

    virtual void print(vector<int> active) {};

    virtual double lossFuncChange(const vector<vector<int> > &connections, const vector<double> & wLambda1, const vector<vector<double> > &wLambda2) { return(-1); };
};

#endif // _QUADRATICDERIVATIVE_H_
