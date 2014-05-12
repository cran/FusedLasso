#ifndef _FUSEDLASSOCOORDINATE_H_
#define _FUSEDLASSOCOORDINATE_H_

#include <iostream>
#include <utility>
#include <vector>
#include "quadraticDerivative.h"
#include "Steps.h"
#include "Graph.h"
#include "SparseMatrix.h"
#include <set>
#include "counted_ptr.h"

using namespace std;


enum penEnum {L1, Huber, L2};

// main class of the program; stores important data and administrate other important classes
class FusedLassoCoordinate {

public:
    vector<double> wObs;
    vector<double> wLambda1Mult; // multiplied by lambda1
    vector<vector<int> > connections; // currently will be initialized to a straight line
    vector<vector<double> > wLambda2Mult; // currently will be initialized to 1, multiplied by lambda2

    // no option yet for weights related to lambda2
    vector<double> beta;
    
    vector<int> active; // which fused groups are active

    // object for calculating the quadratic derivatives under the current fusions
    counted_ptr<QuadraticDerivative> quadratic; 

    int n;
    int p;
    double lambda1;
    double lambda2;
    double accuracy;
    int maxIterInner;
    int maxActivateVars;
    int iterNum;
    bool checkLoss;
    double huberParam;
    penEnum penType;

public:

    /* 
     * the following functions are helper functions for running a single iteration
     */

    // make a single step for variable pos
    void singleStep(int pos);

    // make one iteration for the variable pos
    // returns the L2 distance to the previous iteration
    double singleIteration(const int iterNum);

    // activate inactive variables where appropriate
    // but only up to a certain number of variables
    // returns the number of newly activated variables
    int activateVariables();

    void derivAdjustment(double& adjDeriv, double& zeroPenalty, int pos);

    // deactivate variables (as many as necessary and only if exactly 0)
    // returns the number of deactivated variables
    int deactivateVariables();

    /*
     * The following functions are helper functions for determining
     * the next sets of fused variables
     */

    // for the fused variables adds a single beta to the list of active variables
    bool activateBeta(int pos);

    // for the fused variables, removes a single beta from the list of active variables
    bool deactivateBeta(int pos);

public:

    // Constructor
    // construct with sparse matrix as input
    FusedLassoCoordinate(counted_ptr<QuadraticDerivative> quadDer, const vector<double>& wLambda1, const vector<vector<int> > &connections, const vector<vector<double> > &wLambda2, int maxIterInner, double accuracy, int maxActivateVars, double lambda1, double lambda2, penEnum penType=L1, double huberParam=100);

    // After initialization, run the whole algorithm
    bool runAlgorithm();

    int getIterNum() const;

    vector<double> getBeta() { return beta; };

    vector<double> getBetaOriginal(vector<int>& fusions);

    void printBeta(ostream& outStream);

    void printBetaActive(ostream& outStream); 

    void printDerivActive(ostream& outStream);

    void printDerivs(ostream& outStream);

    void printConnectionsWeights(ostream& outStream);

    void printPosSingleStepInfo(int pos, ostream& outStream);

};

#endif // _FUSEDLASSOCOORDINATE_H_
