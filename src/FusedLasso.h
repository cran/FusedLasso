#ifndef _FUSEDLASSO_H_
#define _FUSEDLASSO_H_

#include <iostream>
#include <utility>
#include <vector>
#include "quadraticDerivativeDiagonal.h"
#include "Steps.h"
#include "Graph.h"
#include "SparseMatrix.h"
#include "FusedLassoCoordinate.h"
#include "counted_ptr.h"

enum regEnum {GAUSSIAN, BINOMIAL, COX};

using namespace std;

// main class of the program; stores important data and administrate other important classes
class FusedLasso {

public:
    SparseMatrix X;
    vector<double> y; // for normal usage
    vector<bool> censored; // for survival use
    vector<double> wObs;
    vector<double> beta;
    
    vector<int> fusions; // stores how the beta are currently fused
    vector<int> fusedGroupSize; // how large is each of the groups

    vector<int> lastDetailedFusions; //stores the result of the last time the detailed fusion was run
    vector<double> wLambda1;
    // object for calculating the quadratic derivatives under the current fusions
    counted_ptr<QuadraticDerivative> quadratic;

    Graph pg;

    int n;
    int p;
    double lambda1;
    double lambda2;

    // these are used to initialize flc
    double accuracy;
    int maxIterInner;
    int maxIterOuter;
    int innerIterNum;
    int outerIterNum;
    int maxActivateVars;
    int fusionLevel;

    regEnum regType;

public:

    // calculate the pull on every variable (not taking into account 
    // variables that are on the same level i.e. could be fused)
    vector<double> getPulls();

    // find the fusions from the univariate tensions
    void getEqualFusions(vector<int>& newFusions, vector<int>& newFusedGroupSize, bool zeroSingle=true, double accFactor=1);
    void getSplitFusionsActive(vector<int>& newFusions, vector<int>& newFusedGroupSize);
    void getSplitFusionsInactive(vector<int>& newFusions, vector<int>& newFusedGroupSize);

    // check if two fusions are equal 
    bool areFusionsEqual(vector<int> &fusion1, vector<int> &fusion2);

    // adjusts the nodePull so that the maximum-flowgraph can be run on it
    void makePullAdjustment(const vector<double>& beta, vector<double>& nodePull, double lambda2);
    
    // Helper function for the constructor
    void createObject(const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2); 

    // given nodes and groups, adds the complementary groups as defined by the given graph
    void addComplementaryGroups(vector<list<int> >& groups, list<int>& nodes);
    // calculat the mean of the current group
    double calcGroupAverage(int pos, vector<double>& nodePull, vector<double>& beta); 

    void sortAllLists(vector<list<int> >& x); 
 
public:

    // Constructor
    // Currently the penalty graph is the default
    // Here regType can be GAUSSIAN or LOGISTIC
    FusedLasso(const double *X, const double* y, const double* wObs, const double* beta, const double* wLambda1, int n, int p, const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2, regEnum regType=GAUSSIAN);

    // construct with vectors as input
    FusedLasso(const vector<double> &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2, regEnum regType=GAUSSIAN);

    // constuct with sparse matrix as input
    FusedLasso(const SparseMatrix &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2, regEnum regType=GAUSSIAN);
  
    ~FusedLasso();

    // translate the beta to a fused beta using the fusions saved in the class
    vector<double> translateOriginalToFused();

    // find the new fusions and set them
    // returns false if they are the same as the old ones
    // i.e. the algorithm has converged
    bool identifyNewFusions(double lastIterChange, int maxFusionLevel);
    bool identifyNewFusionsHuber();

    // initialize the FusedLassoCoordinate object
    // run it to convergence and update all the objects in the current class
    // to the new values of beta
    bool runFused(penEnum penType);

    // run the fused lasso for survival data
    bool runFusedGeneral(penEnum penType);

    // After initialization, run the whole algorithm
    bool runAlgorithm(int maxFusionLevel);
    bool runAlgorithmHuber();
    bool runAlgorithmL2();

    SparseMatrix runAlgorithm(penEnum penType, int maxFusionLevel, const vector<double>& lambda1Vec, const vector<double>& lambda2Vec, const int maxNonZero, vector<bool>& success, vector<int>& outerIterNumVec, vector<int>& innerIterNumVec, bool verbose=false);

    // set new values for lambda1 and lambda2
    void setNewLambdas(double lambda1, double lambda2);

    // getOuterIterNum
    int getOuterIterNum();
    // getInnerIterNum
    int getInnerIterNum();

    // finds the maximum lambda1 at which the first non-exempt variable
    // enters the model
    double findMaxLambda1(const list<int>& exemptVars);

    // finds the value of lambda2 at which the first group of variables breaks apart
    double findMaxLambda2(double lambda1);


    vector<double> getBeta() { return beta; };

    // adds a column to betaMat with the new solution
    void getBeta(SparseMatrix& betaMat);

    void printBetaAndGroup(ostream& outStream);

    void printDerivs(ostream& outStream);

    void printAdjustedDerivs(ostream& outStream);

    void printGroups(vector<list<int> > x, ostream& outStream);

    void checkSolution(ostream& outStream);
};

#endif // _FUSEDLASSO_H_
