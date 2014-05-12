
#include "FusedLassoCoordinate.h"
#include <R_ext/Utils.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <utility>
#include "GeneralFunctions.h"
#include "Exceptions.h"
#include <map>

#define DEBUG false
#define DEBUG_DETAILED false //(lambda1 <= 3.92 && lambda2 <= 1.06)

using namespace std;

FusedLassoCoordinate::FusedLassoCoordinate(counted_ptr<QuadraticDerivative> quadDer, const vector<double>& wLambda1, const vector<vector<int> > &connections, const vector<vector<double> > &wLambda2, int maxIterInner, double accuracy, int maxActivateVars, double lambda1, double lambda2, penEnum penType, double huberParam) : quadratic(quadDer) {
 
    // set all parameters that are easily accessible  
    this->wObs = wObs;
    this->wLambda1Mult.resize(wLambda1.size());
    for(int i = 0; i < wLambda1.size(); ++i) {
        wLambda1Mult[i] = wLambda1[i] * lambda1;
    }
    this->p = quadDer->getp();
    this->n = quadDer->getn();
    this->lambda1 = lambda1;
    this->lambda2 = lambda2;
    this->accuracy = accuracy;
    this->huberParam = huberParam;
    this->penType = penType;
    this->maxIterInner = maxIterInner;
    this->maxActivateVars = maxActivateVars;
    this->beta = quadDer->getBetaVec();
//    this->checkLoss = checkLoss;

    // now initialize the active variables; all that are non-zero;
    active.clear();
    active.reserve(min(n,p));
    for(int i = 0; i < beta.size(); ++i) {
        if(beta[i] != 0) {
            active.push_back(i);
            quadDer->activate(i);
        }
    }

    this->connections = connections;
    this->wLambda2Mult = wLambda2;
    
    // now multiply everything by lambda2
    for(int i = 0; i < wLambda2Mult.size(); ++i) {
        for(int j = 0; j < wLambda2Mult[i].size(); ++j) {
            wLambda2Mult[i][j] *= lambda2;
        }
    }

    // now create the derivative
    quadratic = quadDer;
}



bool FusedLassoCoordinate::activateBeta(int pos) {
    // first check if the element is already included
    vector<int>::iterator vecIt = find(active.begin(), active.end(), pos);
    if(vecIt == active.end()) {// could not find the element
        active.push_back(pos);
        quadratic->activate(pos);
        return true;
    }
    return false;
}


bool FusedLassoCoordinate::deactivateBeta(int pos) {
    vector<int>::iterator vecIt = find(active.begin(), active.end(), pos);
    if(vecIt != active.end()) {
        active.erase(vecIt); // delete the element
        return true;
    }
    return false;
}

void FusedLassoCoordinate::singleStep(int pos) {
    // here, view the problem as only having variable pos
    // find the optimal position and move beta there
    double loss;
    
//    if(DEBUG_DETAILED) {
//        quadratic->print(active);
//    }

    if(DEBUG_DETAILED) {
        cout << "Getting derivative" << endl;
    }
    double deriv = quadratic->getDerivative(pos);
    if(DEBUG_DETAILED) {
        cout << "Getting hessian" << endl;
    }
    double slope = quadratic->getHessian(pos); 

    double oldBeta = beta[pos]; // only needed in conjunction with diagnostics

    if(DEBUG_DETAILED) { // removed 
        cout << "Finding root" << endl;
        cout << "Pos: " << pos << endl;
        cout << "deriv: " << deriv << endl;
        cout << "Slope: " << slope << endl;
        cout << "wlambda1: " << wLambda1Mult[pos] << endl;
        cout << "connections:";
        for(int i = 0; i < connections[pos].size(); ++i) {
            cout << connections[pos][i] << ", ";
        }
        cout << endl;
        cout << "wlambda2: ";
        for(int i = 0; i < wLambda2Mult[pos].size(); ++i) {
            cout << wLambda2Mult[pos][i] << ", ";
        }
        cout << endl;
        cout << "Beta: " << beta[pos] << endl;
//        for(int i = 0; i < beta.size(); ++i) {
//            cout << beta[i] << ", ";
//        }
        cout << endl;
    }
    switch(penType) {
        case L1: Steps::findRootL1(beta, deriv, slope, pos, wLambda1Mult[pos], connections[pos], wLambda2Mult[pos]); break;
        case Huber: Steps::findRootHuber(beta, deriv, slope, pos, wLambda1Mult[pos], connections[pos], wLambda2Mult[pos], huberParam); break;
        case L2: Steps::findRootL2(beta, deriv, slope, pos, wLambda1Mult[pos], connections[pos], wLambda2Mult[pos]); break;
    }
    if(DEBUG_DETAILED) {
        cout << "Updating beta:" << beta[pos] << endl;

    }
    quadratic->updateBeta(pos, beta[pos]);
    if(DEBUG_DETAILED) {
        cout << "After update" << endl;

        cout << "End" << endl;
    }
}

double FusedLassoCoordinate::singleIteration(const int iterNum) {
    // check if there is anything to do first
    if(active.size() == 0) {
        return 0;
    }

    // go through all active variables in single steps
    // calculate the change in L2 norm and return it
    vector<int>::iterator vecIt;
    double maxNorm = 0;
    double oldBeta;

    if(DEBUG_DETAILED) {
        cout << "New single iteration" << endl;
    }

    for(vecIt = active.begin(); vecIt != active.end(); ++vecIt) {
        oldBeta = beta[*vecIt];
        if(DEBUG_DETAILED) { // removed 
            cout << "Step: " << *vecIt << endl;
        }
        singleStep(*vecIt);
        double foo = (beta[*vecIt] - oldBeta);
        if(fabs(foo) > maxNorm) { 
            maxNorm = fabs(foo); 
        }
    }
    if(DEBUG_DETAILED) {
        cout << "Iteration: " << iterNum << " Error: " << maxNorm << endl;
    }
   
    // now normalize the l2norm and return it

    return maxNorm;

}


int FusedLassoCoordinate::activateVariables() {
    // go through all variables that are equal to 0
    // calculate their penalty adjusted derivative and sort them
    // by absolute value of derivative
    // then active the top variables, but at most maxVars
    multimap<double, int> activationCandidates;
    vector<double> derivVec = quadratic->getDerivativeVec(); 
    vector<bool> isActive(beta.size(), false);
    for(int i = 0; i < active.size(); ++i) {
        isActive[active[i]] = true;
    }

    double adjDeriv;
    for(int pos = 0; pos < beta.size(); ++pos) {
        if(!isActive[pos]) {
            double zeroPenalty = wLambda1Mult[pos];

            adjDeriv = derivVec[pos];
            derivAdjustment(adjDeriv, zeroPenalty, pos);
            // now if in absolute value > lambda1, add to list of
            // possible activations
            if(fabs(adjDeriv) > zeroPenalty) {
                activationCandidates.insert(make_pair(fabs(adjDeriv), pos));
            }
        }
    }

    // now go through the activation candidates
    multimap<double, int>::reverse_iterator mapIt;
    int activateCount;
    for(mapIt = activationCandidates.rbegin(), activateCount = 0; mapIt !=activationCandidates.rend() && activateCount < maxActivateVars; ++mapIt, ++ activateCount) {
        active.push_back(mapIt->second);       
        quadratic->activate(mapIt->second); 
        if(DEBUG) {
            cout << "Activate: " << mapIt->second << " Deriv: " << mapIt->first << endl;
        }
    }

    // and sort the vector of active variables
    sort(active.begin(), active.end());

    return activateCount;
}

void FusedLassoCoordinate::derivAdjustment(double& adjDeriv, double& zeroPenalty, int pos) {
    if(penType == L1) {
        for(int i = 0; i < connections[pos].size(); ++i) {
            if(beta[connections[pos][i]] > 0) {
                adjDeriv -= wLambda2Mult[pos][i];
            }
            else if(beta[connections[pos][i]] < 0) {
                adjDeriv += wLambda2Mult[pos][i];
            }
            else {
                zeroPenalty += wLambda2Mult[pos][i];
            }
        }
        return;
    }
    
    if(penType == Huber) {
        for(int i = 0; i < connections[pos].size(); ++i) {
            if(beta[connections[pos][i]] > 1/huberParam) {
                adjDeriv -= wLambda2Mult[pos][i];
            }
            else if(beta[connections[pos][i]] < -1/huberParam) {
                adjDeriv += wLambda2Mult[pos][i];
            }
            else{
                adjDeriv -= huberParam * beta[connections[pos][i]] * wLambda2Mult[pos][i];
            }
        }
        return;
    }

    if(penType == L2) {
        for(int i = 0; i < connections[pos].size(); ++i) {
            adjDeriv -= beta[connections[pos][i]] * 2 * wLambda2Mult[pos][i];
        }
        return;
    }
    error("Should by one of the previous lTypes\n");

}


int FusedLassoCoordinate::deactivateVariables() {
    // go through all active variables and remove the
    // ones from being active that are equal to 0

    int deleteCount = 0;
    vector<int>::iterator vecIt;
    for(vecIt = active.begin(); vecIt != active.end(); ) {
       if(beta[*vecIt] == 0) { // yes, deactivate
            if(DEBUG) {
                cout << "Deactivate: " << *vecIt << " ";
            }
          vecIt = active.erase(vecIt); // due to deletion don't need to move one element ahead
          deleteCount++;
       }
       else {
           vecIt++; 
       }
    }
    if(DEBUG) {
        cout << endl;
    }
    return deleteCount;
}

bool FusedLassoCoordinate::runAlgorithm() {
    // in here at the end the beta and unfused object should be updated
    // then all other public functions can always consider this
    // as being up to date

    bool finished = false;
    int newActivations;
    int deactivations;
    iterNum = 0;
    double curError;

    if(DEBUG) {
        cout << "New runAlgorithm" << endl;
    }

    while(!finished && iterNum < maxIterInner) {
        curError = 1;
        
        if(DEBUG_DETAILED) {
            printBetaActive(cout);
            printDerivActive(cout);
        }

        while(curError > accuracy && iterNum < maxIterInner) {
//            if(checkLoss & iterNum % 100 == 0) {
//                double change = quadratic->lossFuncChange(connections, wLambda1Mult, wLambda2Mult);
//                if(change > -accuracy) {
//                    break;
//                }
//            }

            if(DEBUG && (iterNum % 100 == 0)) {
                cout << "Iteration started" << endl;
            }
            curError = singleIteration(iterNum);
            R_CheckUserInterrupt(); // checks if the user has cancelled the computation
            if(DEBUG && (iterNum % 100 == 0)) {
                cout << "Iteration: " << iterNum << " Error: " << curError << endl;
                if(curError > 1e5) {
                    printBeta(cout);            
                }
            }
            iterNum++;
            if(curError > 1e5) {
                return(false);
            }
        }

//        deactivations = deactivateVariables();
        deactivations = 0;
        newActivations = activateVariables();

        if(DEBUG_DETAILED) {       
            cout << "Deactivate: " << deactivations << " Activations: " << newActivations << endl;
        }
//        checkSolution(cout);

        if(newActivations == 0) { // has converged
            finished = true;
        }
        if( iterNum >= maxIterInner) {
            finished = false;
        }
    }

    if(iterNum>=maxIterInner && DEBUG) {
        printBetaActive(cout);
        printDerivActive(cout);
    }
    if(DEBUG) {
        cout << "End run Algorithm with iterations " << iterNum << endl;
    }
//    cout << "Number of inner iterations: " << curInnerIter << endl;
//    cout << "Accuracy: " << accuracy << " Error: " << curError << endl;

    return(finished);
}

int FusedLassoCoordinate::getIterNum() const {
    return iterNum;
}


vector<double> FusedLassoCoordinate::getBetaOriginal(vector<int>& fusions) {
    vector<double> betaOrig(fusions.size());

    for(int i = 0; i < fusions.size(); ++i) {
        if(fusions[i] > beta.size() - 1) {
            throw OutOfBoundsException(i);
        }
        betaOrig[i] = beta[fusions[i]];
    }
    return betaOrig;
}

void FusedLassoCoordinate::printBetaActive(ostream& outStream) {
    outStream << "===================== Beta Active =====================" << endl;
    for(int i = 0; i < beta.size(); ++i) {
        if(beta[i] != 0) {
            outStream << i << ":" << beta[i] << " | ";
        }
    }
    outStream << endl << "======================== End Beta Active ================" << endl;
}

void FusedLassoCoordinate::printDerivActive(ostream& outStream) {
    outStream << "===================== Deriv Active =====================" << endl;
    for(int i = 0; i < beta.size(); ++i) {
        if(beta[i] != 0) {
            outStream << i << ":" << quadratic->getDerivative(i) << " | ";
        }
    }
    outStream << endl << "======================== End Deriv Active ================" << endl;
}




void FusedLassoCoordinate::printBeta(ostream& outStream) {
    outStream << "======================= Beta ===============" << endl;
    for(int i = 0; i < beta.size(); ++i) {
        outStream << beta[i] << ", ";
    }
    outStream << endl << "================ End Beta ================" << endl;
}

void FusedLassoCoordinate::printDerivs(ostream& outStream) {
    outStream << "============= Derivs =============" << endl;
    for(int i = 0; i < beta.size(); ++i) {
        outStream << quadratic->getDerivative(i) << " " ;
    }
    outStream << endl << "============= End Derivs =========" << endl;
}


void FusedLassoCoordinate::printConnectionsWeights(ostream& outStream) {
    outStream << "============ Conn Weights ============" << endl;
    for(int i = 0; i < connections.size(); ++i) {
        outStream << "Node:" << i << endl; 
        for(int j = 0; j < connections[i].size(); ++j) {
            outStream << connections[i][j] << " " << wLambda2Mult[i][j] << ";";
        }
        outStream << endl;
    }

}



void FusedLassoCoordinate::printPosSingleStepInfo(int pos, ostream& outStream) {
    outStream << "Info at Step: " << pos << endl;
    outStream << "Deriv: " << quadratic->getDerivative(pos) << endl;    
    outStream << "Hessian: " << quadratic->getHessian(pos) << endl;
    outStream << "OldBeta: " << beta[pos] << endl;
    outStream << "wLambda1: " << wLambda1Mult[pos] << endl;
    outStream << "Connections: " << endl;
    printVector(connections[pos], outStream);
    outStream << "wLambda2Mult: " << endl;
    printVector(wLambda2Mult[pos], outStream);

}

    

