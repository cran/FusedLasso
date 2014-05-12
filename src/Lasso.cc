#include "Lasso.h"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <map>
#include <utility>
#include "GeneralFunctions.h"
#include "Exceptions.h"
#include "quadraticDerivativeDiagonal.h"

using namespace std;

Lasso::Lasso(const double* X, const double* y, const double*  wObs, const double* beta, const double* wLambda1, int n, int p, int maxIterInner, double accuracy, int maxActivateVars, double lambda1) {
    // set number of rows and cols and store X, y and w
    (this->y).resize(n, 0);
    (this->wObs).resize(n, 0);
    (this->wLambda1).resize(p, 0);
    (this->beta).resize(p, 0);

    // now copy the elements of the vectors over and create the Sparse Matrix
    this->X = SparseMatrix(X, n, p);

    for(int i = 0; i < n; ++i) {
        this->y[i] = y[i];
        this->wObs[i] = wObs[i];
    }

    for(int i = 0; i < p; ++i) {
        this->beta[i] = beta[i];
        this->wLambda1[i] = wLambda1[i];
    }

    this->n = n;
    this->p = p;

    createObject(maxIterInner, accuracy, maxActivateVars, lambda1);
}



Lasso::Lasso(const vector<double> &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, int maxIterInner, double accuracy, int maxActivateVars, double lambda1) {
    this->n = y.size();
    this->p = beta.size();
    this->X = SparseMatrix(X, n, p);
    this->y = y;
    this->wObs = wObs;
    this->beta = beta;
    this->wLambda1 = wLambda1;

    createObject(maxIterInner, accuracy, maxActivateVars, lambda1);
}

Lasso::Lasso(const SparseMatrix &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, int maxIterInner, double accuracy, int maxActivateVars, double lambda1) {
    
    this->X = X;
    this->y = y;
    this->wObs = wObs;
    this->beta = beta;
    this->n = X.n;
    this->p = X.p;
    this->wLambda1 = wLambda1;

    createObject(maxIterInner, accuracy, maxActivateVars, lambda1);
}

// does all the rest of the work after X, y, w and beta have been set
void Lasso::createObject(int maxIterInner, double accuracy, int maxActivateVars, double lambda1) {
    this->lambda1 = lambda1;

    this->maxIterInner = maxIterInner;
    this->accuracy = accuracy;
    this->maxActivateVars = maxActivateVars;

    // now set the quadratic object
    quadratic = counted_ptr<QuadraticDerivative>(new QuadraticDerivativeDiagonal(this->X, this->y, this->wObs, this->beta));
}

Lasso::~Lasso() {
}

bool Lasso::runAlgorithm() {
    vector<vector<int> > connections(p, vector<int>());
    vector<vector<double> > wLambda2(p, vector<double>());

    FusedLassoCoordinate flc(quadratic, wLambda1, connections, wLambda2, maxIterInner, accuracy, maxActivateVars, lambda1, 0);
   
    bool result = flc.runAlgorithm(); 
    beta = flc.getBeta();

    return(result);
}


SparseMatrix Lasso::runAlgorithm(const vector<double>& lambda1Vec, vector<bool>& success) {
    SparseMatrix betaSols(beta.size());
    success.clear(); success.resize(lambda1Vec.size());

    for(int i = 0; i < lambda1Vec.size(); ++i) {
//        cout << "Working on Lambda1: " << lambda1Vec[i] << endl;
        
        setNewLambdas(lambda1Vec[i]);
        success[i] = runAlgorithm();
        betaSols.addColumn(getBeta());
    }
    return betaSols;
}



void Lasso::setNewLambdas(double lambda1) {
    // when setting the new lambdas, only the lambdas
    // have to be reset; no other changes are necessary
    this->lambda1 = lambda1;
}

void Lasso::getBeta(SparseMatrix& betaMat) {
    betaMat.addColumn(beta);
}


