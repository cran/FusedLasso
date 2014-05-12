#ifndef _STEPS_H_
#define _STEPS_H_

#include <utility>
#include <vector>

using namespace std;

/*
 * Given the current betas, it finds the root of the function
 */

class Steps {
inline static bool updateBeta(vector<double>& beta, const int pos, double& deriv, const double& slope, const double& zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepSize);
 
public:
    // most important application -> find the root

    static void findRoot(vector<double> &beta, double derivQuad, const double slope, const int pos, const double zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepSize);
    
    static void findRootL1(vector<double> &beta, double derivQuad, const double hessian, const int pos, const double zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepSize);
    
    static void findRootL2(vector<double> &beta, double derivQuad, const double hessian, const int pos, const double zeroStepSize, const vector<int>& l2Idx, const vector<double>& l2Mult);

    static void findRootHuber(vector<double> &beta, double derivQuad, const double hessian, const int pos, const double zeroStepSize, const vector<int>& huberIdx, const vector<double>& huberMult, double huberParam);

};

#endif // _STEPS_H_
