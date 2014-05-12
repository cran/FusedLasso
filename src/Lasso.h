#ifndef _LASSO_H_
#define _LASSO_H_

#include <iostream>
#include <utility>
#include <vector>
#include "quadraticDerivativeDiagonal.h"
#include "Steps.h"
#include "Graph.h"
#include "SparseMatrix.h"
#include "FusedLassoCoordinate.h"
#include "counted_ptr.h"

using namespace std;

// main class of the program; stores important data and administrate other important classes
class Lasso {

public:
    SparseMatrix X;
    vector<double> y; // for normal usage
    vector<double> wObs;
    vector<double> beta;
    
    vector<double> wLambda1;
    // object for calculating the quadratic derivatives under the current fusions
    counted_ptr<QuadraticDerivative> quadratic;

    int n;
    int p;
    double lambda1;
    double lambda2;

    // these are used to initialize flc
    double accuracy;
    int maxIterInner;
    int maxActivateVars;
public:

    // calculate the pull on every variable (not taking into account 
    // variables that are on the same level i.e. could be fused)

    // Helper function for the constructor
    void createObject(int maxIterInner, double accuracy, int maxActivateVars, double lambda1); 

public:

    // Constructor
    // Currently the penalty graph is the default
    Lasso(const double *X, const double* y, const double* wObs, const double* beta, const double* wLambda1, int n, int p, int maxIterInner, double accuracy, int maxActivateVars, double lambda1);

    // construct with vectors as input
    Lasso(const vector<double> &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, int maxIterInner, double accuracy, int maxActivateVars, double lambda1);

    // constuct with sparse matrix as input
    Lasso(const SparseMatrix &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, int maxIterInner, double accuracy, int maxActivateVars, double lambda1);
  
    ~Lasso();

    // After initialization, run the whole algorithm
    bool runAlgorithm();

    SparseMatrix runAlgorithm(const vector<double>& lambda1Vec, vector<bool>& success);

    // set new values for lambda1 and lambda2
    void setNewLambdas(double lambda1);

    vector<double> getBeta() { return beta; };

    // adds a column to betaMat with the new solution
    void getBeta(SparseMatrix& betaMat);
};

#endif // _LASSO_H_

