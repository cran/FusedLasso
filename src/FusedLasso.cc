#include "FusedLasso.h"
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <map>
#include <utility>
#include "GeneralFunctions.h"
#include "Exceptions.h"
#include "quadraticDerivativeDiagonal.h"
#include "quadraticDerivativeLogistic.h"

#define DEBUG false
#define DEBUG2 true
#define DEBUG_DETAILED false
#define DEBUG_VERY_DETAILED false


using namespace std;

FusedLasso::FusedLasso(const double* X, const double* y, const double*  wObs, const double* beta, const double* wLambda1, int n, int p, const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2, regEnum regType) {
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
    
    this->regType = regType;

    createObject(graph, maxIterInner, maxIterOuter, accuracy, maxActivateVars, lambda1, lambda2);
}



FusedLasso::FusedLasso(const vector<double> &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2, regEnum regType) {
    this->n = y.size();
    this->p = beta.size();
    this->X = SparseMatrix(X, n, p);
    this->y = y;
    this->wObs = wObs;
    this->beta = beta;
    this->wLambda1 = wLambda1;

    this->regType = regType;

    createObject(graph, maxIterInner, maxIterOuter, accuracy, maxActivateVars, lambda1, lambda2);
}

FusedLasso::FusedLasso(const SparseMatrix &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2, regEnum regType) {
    
    this->X = X;
    this->y = y;
    this->wObs = wObs;
    this->beta = beta;
    this->n = X.n;
    this->p = X.p;
    this->wLambda1 = wLambda1;

    this->regType = regType;

    createObject(graph, maxIterInner, maxIterOuter, accuracy, maxActivateVars, lambda1, lambda2);
}


// does all the rest of the work after X, y, w and beta have been set
void FusedLasso::createObject(const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2) {
    this->pg = graph;
    this->lambda1 = lambda1;
    this->lambda2 = lambda2;

    this->maxIterInner = maxIterInner;
    this->maxIterOuter = maxIterOuter;
    this->accuracy = accuracy;
    this->maxActivateVars = maxActivateVars;

    // set the beta and the default for the fusions
    this->fusedGroupSize.resize(p, 1);
    this->fusions.clear(); this->fusions.resize(p);
    for(int i = 0; i < p; ++i) {
            this->fusions[i] = i;
    }

    // now set the quadratic object
    switch(regType) {
        case GAUSSIAN: 
            quadratic = counted_ptr<QuadraticDerivative>(new QuadraticDerivativeDiagonal(this->X, this->y, this->wObs, this->beta));
            break;
        case BINOMIAL:
            quadratic = counted_ptr<QuadraticDerivative>(new QuadraticDerivativeLogistic(this->X, this->y, this->wObs, this->beta));
    }

    lastDetailedFusions  = fusions;
    fusionLevel = 0;

//    graph.printGraph(cout);
}


FusedLasso::~FusedLasso() {
}

// assumes that beta has been updated from the fusedBeta before -> translateFusedToOriginal()
vector<double> FusedLasso::getPulls() {
    vector<double> pulls(beta.size());
    for(int i=0; i < pulls.size(); ++i) {
        pulls[i] = quadratic->getDerivative(i);
    }
    
    return pulls;
}

// gives the current fusions based on the current betas (no splits will be performed)
void FusedLasso::getEqualFusions(vector<int>& newFusions, vector<int>& newFusedGroupSize, bool zeroSingle, double accFactor) {
    newFusions.clear(); newFusedGroupSize.clear();
    newFusions.resize(beta.size(), -1);
    list<int> nodes;
    list<int>::iterator listIt;
    int groupNum = 0;

    double accHere = accuracy * accFactor;

    for(int i = 0; i < beta.size(); ++i) {
        if(newFusions[i] == -1) {
            nodes = pg.connectedWithSameValue(i, beta, accHere);
            if(zeroSingle && beta[i]==0) {
                for(listIt = nodes.begin(); listIt != nodes.end(); ++listIt) {
                    newFusions[*listIt] = groupNum;
                    groupNum++;
                    newFusedGroupSize.push_back(1);
                }
            }
            else {
                newFusedGroupSize.push_back(nodes.size());
                for(listIt = nodes.begin(); listIt != nodes.end(); ++listIt) {
                    newFusions[*listIt] = groupNum;
                }
                ++groupNum;
            }
        }
    }
}

bool FusedLasso::areFusionsEqual(vector<int> &fusion1, vector<int> &fusion2) {
    if(fusion1.size() != fusion2.size()) {
        return false;
    }

    for(int i = 0; i < fusion1.size(); ++i) {
        if(fusion1[i] != fusion2[i]) {
            if(DEBUG_VERY_DETAILED) {
	      REprintf("============ ARE FUSIONS EQUAL ====================\n %d\n ========== END ARE FUSIONS EQUAL ==========\n", i);
            }
            return(false);
        }
    }

    return true;
}

void FusedLasso::getSplitFusionsActive(vector<int>& newFusions, vector<int>& newFusedGroupSize) {

    newFusions.clear(); newFusedGroupSize.clear();
    newFusions.resize(p, -1);
    
    // first, we need to get the pulls correctly
    vector<double> nodePull = getPulls();
    // adjust the nodePulls appropriately
    makePullAdjustment(beta, nodePull, lambda2);
    vector<list<int> > groups;
    list<int> nodes;
    list<int>::iterator listIt;
    int newNumFused = 0;

    for(int i = 0; i < p; ++i) {
        if(newFusions[i] == -1) { // not yet grouped
            nodes = pg.connectedWithSameValue(i, beta, accuracy);
            if(beta[i] == 0) { // all variables are singles
                for(listIt = nodes.begin(); listIt != nodes.end(); ++listIt) {
                    beta[*listIt] = 0; // sets some betas not more than accuracy different
                    // from 0 back to zero
                    newFusedGroupSize.push_back(1);
                    newFusions[*listIt] = newNumFused;
                    newNumFused++;
                }
                continue;
            }
            else if(beta[i] > 0) {
                groups = pg.splitGroup(nodes, nodePull, lambda1, lambda2, false);
            }
            else {
                groups = pg.splitGroup(nodes, nodePull, -lambda1, lambda2, false);
            }

            try {
            addComplementaryGroups(groups, nodes);
            }
            catch(DoubleGroupedException e) {
	      //                cout << "Beta: " << beta[i] << endl;
              //  cout << "Lambda1: " << lambda1 << endl;
              //  cout << "Lambda2: " << lambda2 << endl;
              //  pg.splitGroup(nodes, nodePull, -lambda1, lambda2, false);
              //  cout << "====================== First Graph ====================" << endl;
              //  pg.printGraph(cout);
              //  pg.splitGroup(nodes, nodePull, lambda1, lambda2, true);
              //  cout << "====================== Second Graph ====================" << endl;
              //  pg.printGraph(cout);
                
              //  nodePull = getPulls();
              //  cout << "NodePull (no Adj:)" << endl;
              //  for(listIt = nodes.begin(); listIt != nodes.end(); ++listIt) {
              //      cout << *listIt << ": " << nodePull[*listIt] << endl;
              //  }
              //  makePullAdjustment(beta, nodePull, lambda2);
              //  nodes = pg.connectedWithSameValue(i, beta, accuracy);
              //  cout << "NodePull: " << endl;
              //  for(listIt = nodes.begin(); listIt != nodes.end(); ++listIt) {
              //      cout << *listIt << ": " << nodePull[*listIt] << endl;
              //  }
              //  cout << endl;
              //  throw e;
            }

            if(DEBUG_VERY_DETAILED && groups.size() > 1) {
              //  cout << "======== SPLIT =========" << endl;
            }

            // now mark all the sets in here
            for(int j = 0; j < groups.size(); ++j) {
                if(DEBUG_VERY_DETAILED && groups.size() > 1) {
                  //  cout << "GROUP " << j << " " << groups[j].front() << " " << groups[j].back() << endl;
                }

                // setting this here is not completely correct 
                newFusedGroupSize.push_back(groups[j].size());
                for(listIt = groups[j].begin(); listIt != groups[j].end(); ++listIt) {
                    newFusions[*listIt] = newNumFused;
                }
                newNumFused++;
            }
            if(DEBUG_VERY_DETAILED && groups.size() > 1) {
              //  cout << "============= END SPLIT =========" << endl;
            }
    
        }
    }
}

void FusedLasso::getSplitFusionsInactive(vector<int>& newFusions, vector<int>& newFusedGroupSize) {

    newFusions.clear(); newFusedGroupSize.clear();
    newFusions.resize(p, -1);
    
    // first, we need to get the pulls correctly
    vector<double> nodePull = getPulls();
    // adjust the nodePulls appropriately
    makePullAdjustment(beta, nodePull, lambda2);
    vector<list<int> > groups;
    list<int> nodes;
    list<int>::iterator listIt;
    int newNumFused = 0;

    for(int i = 0; i < p; ++i) {
        if(newFusions[i] == -1) { // not yet grouped
            nodes = pg.connectedWithSameValue(i, beta, accuracy);
            if(beta[i] == 0) {
                vector<list<int> > foo;
                groups = pg.splitGroup(nodes, nodePull, -lambda1, lambda2, false);
                foo = pg.splitGroup(nodes, nodePull, lambda1, lambda2, true);
                groups.insert(groups.end(), foo.begin(), foo.end());
            }
            else { // leave everything the same as before
                groups.resize(1);
                groups[0] = nodes;
            }

            try {
            addComplementaryGroups(groups, nodes);
            }
            catch(DoubleGroupedException e) {
              /*  cout << "Beta: " << beta[i] << endl;
                cout << "Lambda1: " << lambda1 << endl;
                cout << "Lambda2: " << lambda2 << endl;
                pg.splitGroup(nodes, nodePull, -lambda1, lambda2, false);
                cout << "====================== First Graph ====================" << endl;
                pg.printGraph(cout);
                pg.splitGroup(nodes, nodePull, lambda1, lambda2, true);
                cout << "====================== Second Graph ====================" << endl;
                pg.printGraph(cout);
                
                nodePull = getPulls();
                cout << "NodePull (no Adj:)" << endl;
                for(listIt = nodes.begin(); listIt != nodes.end(); ++listIt) {
                    cout << *listIt << ": " << nodePull[*listIt] << endl;
                }
                makePullAdjustment(beta, nodePull, lambda2);
                nodes = pg.connectedWithSameValue(i, beta, accuracy);
                cout << "NodePull: " << endl;
                for(listIt = nodes.begin(); listIt != nodes.end(); ++listIt) {
                    cout << *listIt << ": " << nodePull[*listIt] << endl;
                }
                cout << endl; */
                throw e;
            }

            if(DEBUG_VERY_DETAILED && groups.size() > 1) {
	      // cout << "======== SPLIT =========" << endl;
            }

            // now mark all the sets in here
            for(int j = 0; j < groups.size(); ++j) {
                if(DEBUG_VERY_DETAILED && groups.size() > 1) {
		  // cout << "GROUP " << j <<  " " << groups[j].front() << " " << groups[j].back() << endl;
                }
                // setting this here is not completely correct 
                newFusedGroupSize.push_back(groups[j].size());
                for(listIt = groups[j].begin(); listIt != groups[j].end(); ++listIt) {
                    newFusions[*listIt] = newNumFused;
                }
                newNumFused++;
            }
            if(DEBUG_VERY_DETAILED && groups.size() > 1) {
	      // cout << "============= END SPLIT =========" << endl;
            }
    
        }
    }
}


void FusedLasso::sortAllLists(vector<list<int> >& x) {
    for(int i = 0; i < x.size(); ++i) {
        x[i].sort();
    }
}


bool listComp(const list<int>& list1, const list<int>& list2) {
    int val1, val2;
    val1 = list1.front();
    val2 = list2.front();
    return val1 < val2;
}

void FusedLasso::addComplementaryGroups(vector<list<int> >& groups, list<int>& nodes) {
    list<int> complNodes = nodes;
    vector<list<int> > foo;
    vector<list<int> > groupsCopy = groups;
    vector<list<int> > groupsSave = groups;
    list<int> groupsMerged;
    list<int>::iterator listIt;

    // first find the complementary nodes
    for(int i = 0; i < groups.size(); ++i) {
            groupsMerged.splice(groupsMerged.end(), groupsCopy[i]);
    }

    groupsMerged.sort(); 
    int groupsMergedSize = groupsMerged.size();
    groupsMerged.unique();
    if(groupsMergedSize > groupsMerged.size()) {
      /* cout << "Nodes: "; 
        printList(nodes, cout);
        for(int i = 0; i < groupsSave.size(); ++i) {
            cout << endl;
            printList(groupsSave[i], cout);
        }
        cout << endl; */
        throw DoubleGroupedException(0); 
    }

    complNodes = pg.getComplement(groupsMerged, nodes);

    // now make the rest into groups and add it to the groups
    foo = pg.identifyConnectedGroups(complNodes);
    groups.insert(groups.end(), foo.begin(), foo.end()); 
    sortAllLists(groups);
    sort(groups.begin(), groups.end(), listComp); 
}



vector<double> FusedLasso::translateOriginalToFused() {
    vector<double> fusedBeta(fusedGroupSize.size());
    for(int i = 0; i < fusions.size(); ++i) {
        fusedBeta[fusions[i]] = beta[i];
    }
    return fusedBeta;
}

bool FusedLasso::identifyNewFusionsHuber() {
    vector<int> newFusions;
    vector<int> newFusedGroupSize;

    double huberParam = 1000;

    // create the object to run the coordinate descent algorithm
    vector<int> singleFusions(p);
    for(int i = 0; i < p; ++i) {
        singleFusions[i] = i;
    }

    vector<vector<int> > connectionsSingle;
    vector<vector<double> > wLambda2Single;

    pg.getFusedConnectionsWeights(singleFusions, p, connectionsSingle, wLambda2Single);

    // now first run Huber, then normal
    FusedLassoCoordinate flcHuber(quadratic, wLambda1, connectionsSingle, wLambda2Single, 100, accuracy, 100000, lambda1, lambda2, Huber, huberParam);
    bool resultHuber = flcHuber.runAlgorithm();
//    cout << "Huber" << endl;
//    printVector(flcHuber.getBetaOriginal(singleFusions), cout);
//    printVector(quadratic->getDerivativeVec(), cout);

    FusedLassoCoordinate flc(quadratic, wLambda1, connectionsSingle, wLambda2Single, maxIterInner, accuracy, maxActivateVars, lambda1, lambda2, L1);
    bool result = flc.runAlgorithm();
    beta = flc.getBetaOriginal(singleFusions);
//    cout <<" Single " << endl;
//    printVector(beta, cout);

    getEqualFusions(newFusions, newFusedGroupSize);
//    cout << "Fusions" << endl;
//    printVector(newFusions, cout);
    if(areFusionsEqual(fusions, newFusions)) {
        return(false);
    }
    else {
        fusions = newFusions;
        fusedGroupSize = newFusedGroupSize;
        return(true);
    }

}


bool FusedLasso::identifyNewFusions(double lastIterChange, int maxFusionLevel) {
    // first, find the new fusions, compare them to the old ones
    // and set them if they are different
    vector<int> newFusions;
    vector<int> newFusedGroupSize;
  
    if(maxFusionLevel==-1) {
        return(false); // nothing should be done
    }

    if(lastIterChange > accuracy) {
        fusionLevel = 0;
    }
    else {
        ++fusionLevel;
    }

    if(DEBUG_DETAILED) {
      // cout << "Change last iteration: " << lastIterChange << " Level: " << fusionLevel << endl;
    }

    if(fusionLevel == 0 && fusionLevel <= maxFusionLevel) {
//        cout << "Level 0" << endl;
        getEqualFusions(newFusions, newFusedGroupSize);
        if(areFusionsEqual(fusions, newFusions)) {
            ++fusionLevel;
        }
    }
    
    if(fusionLevel == 1 && fusionLevel <= maxFusionLevel) {
//        cout << "Level 1" << endl;
        getSplitFusionsInactive(newFusions, newFusedGroupSize);
        if(areFusionsEqual(fusions, newFusions)) {
            ++fusionLevel;
        }
    }

    if(fusionLevel == 2 && fusionLevel <= maxFusionLevel) {
//        cout << "Level 2" << endl;
        // getSplitFusionsActive(newFusions, newFusedGroupSize);
        getSplitFusionsActive(newFusions, newFusedGroupSize);
        if(areFusionsEqual(fusions, newFusions)) {
            ++fusionLevel;
        }
    } 
    
    if(fusionLevel > maxFusionLevel || fusionLevel == 3) {
        return(false);
    }
    else {
        fusions = newFusions;
        fusedGroupSize = newFusedGroupSize;
        return(true);
    }
}


bool FusedLasso::runFused(penEnum penType) {
    // create the object to run the coordinate descent algorithm
    vector<double> betaFused = translateOriginalToFused();
    vector<vector<int> > fusedGroups(fusedGroupSize.size());
    for(int i = 0; i < fusions.size(); ++i) {
        fusedGroups[fusions[i]].push_back(i);
    }

    SparseMatrix fusedX  = X.createFusedX(fusedGroups);
    vector<double> wLambda1Coordinate(fusedGroupSize.size());
    for(int i = 0; i < wLambda1Coordinate.size(); ++i) {
        wLambda1Coordinate[i] = 0;
        for(int j = 0; j < fusedGroups[i].size(); ++j) {
            wLambda1Coordinate[i] += wLambda1[fusedGroups[i][j]];
        }
    }

    vector<vector<int> > connectionsFused;
    vector<vector<double> > wLambda2Fused;

    pg.getFusedConnectionsWeights(fusions, fusedGroupSize.size(), connectionsFused, wLambda2Fused);

    counted_ptr<QuadraticDerivative> quadDer(new QuadraticDerivativeDiagonal(fusedX, y, wObs, betaFused));
    FusedLassoCoordinate flc(quadDer, wLambda1Coordinate, connectionsFused, wLambda2Fused, maxIterInner, accuracy/10, maxActivateVars, lambda1, lambda2, penType);

    bool result = flc.runAlgorithm();
    innerIterNum += flc.getIterNum();

    // now translate the result back
    beta = flc.getBetaOriginal(fusions);
    for(int i = 0; i< beta.size(); ++i) {
        quadratic->updateBeta(i, beta[i]);
    }
    outerIterNum++;
    return result;

}


bool FusedLasso::runFusedGeneral(penEnum penType) {
  // leave if there are already too many interations
  // cout << "Start runFusedGeneral" << endl;
  if(outerIterNum > maxIterOuter) {
    return(false);
  }
  //  printVector(beta, cout);

    vector<double> betaFused = translateOriginalToFused();
    vector<vector<int> > fusedGroups(fusedGroupSize.size());
    for(int i = 0; i < fusions.size(); ++i) {
        fusedGroups[fusions[i]].push_back(i);
    }

    if(DEBUG) {
      // cout << "Before createFusedX" << endl;
    }
    SparseMatrix fusedX = X.createFusedX(fusedGroups);
    vector<double> wLambda1Fused(fusedGroupSize.size());
    for(int i = 0; i < fusedGroups.size(); ++i) {
        wLambda1Fused[i] = 0;
        for(int j = 0; j < fusedGroups[i].size(); ++j) {
            wLambda1Fused[i] += wLambda1[fusedGroups[i][j]];
        }
    }

    vector<vector<int> > connectionsFused;
    vector<vector<double> > wLambda2Fused;

    if(DEBUG) {
      // cout << "Before getConnectionWeights" << endl;
    }
    pg.getFusedConnectionsWeights(fusions, fusedGroupSize.size(), connectionsFused, wLambda2Fused);
    if(DEBUG) {
      // cout << "After getConnectionWeights" << endl;
    }
    if(DEBUG) {
      // cout << "Before define variables" << endl;
    }
    vector<double> oldBeta;
    double maxBetaChange = accuracy + 1;
    counted_ptr<QuadraticDerivative> quadDer;
    counted_ptr<FusedLassoCoordinate> flc;
    bool result = true;
    
    if(DEBUG) {
      /* cout << "After define variables" << endl;

      cout << "maxBetaChange " << maxBetaChange << endl;
      cout << "Acccuracy " << accuracy << endl;
      cout << "Result " << result << endl;
      cout << "outerIterNum " << outerIterNum << endl;
      cout << "maxIterOuter " << maxIterOuter << endl; */
    }

    while(maxBetaChange > accuracy && result && outerIterNum <= maxIterOuter) {
        oldBeta = betaFused;
	if(DEBUG) {
	  //cout << "Entered while loop in runFusedGeneral" << endl;
	}

        if(regType==BINOMIAL) {
	  if(DEBUG) {
	    //cout << "Before new QuadraticDerivative" << endl;
	  }
            quadDer = counted_ptr<QuadraticDerivative>(new QuadraticDerivativeLogistic(fusedX, y, wObs, betaFused));
	    if(DEBUG) {
	      //cout << "After new QuadraticDerivative" << endl;
	    }
        }
        else {
	  error("Invalid regression type\n");exit(1);
        }

        if(quadDer -> isExtreme()) {
            return(false);
        }

        flc = counted_ptr<FusedLassoCoordinate> (new FusedLassoCoordinate(quadDer, wLambda1Fused, connectionsFused, wLambda2Fused, maxIterInner, accuracy/10 , maxActivateVars, lambda1, lambda2, penType));
        if(DEBUG) { 
	  REprintf("Before runAlgorithm()\n");
        }
        result =  flc->runAlgorithm();
        if(DEBUG) {
	  REprintf("After runAlgorithm()\n");
        }
        betaFused = flc->getBeta();
        innerIterNum += flc->getIterNum();

        maxBetaChange = maxDiffDoubleVec(oldBeta, betaFused);
        outerIterNum++;
        if(DEBUG_DETAILED) {
	  /* cout << "runFusedGeneral: " << outerIterNum << " Change " << maxBetaChange << " Result: " << result  << "Acc: " << accuracy << endl;
            cout << "NumNonZero: " << numNonZero(flc->getBetaOriginal(fusions)) << endl;
            cout << "Beta";
	            printVector(betaFused, cout); */
        }
    }

    // cout << "Before get beta original" << endl;
    beta = flc->getBetaOriginal(fusions);
   
    // cout << "After get beta original" << endl;
    if(regType==BINOMIAL) {
      // cout << "Before get quadratic derivative" << endl;
        quadratic = counted_ptr<QuadraticDerivative>(new QuadraticDerivativeLogistic(X, y, wObs, beta)); 
	// cout << "after get quadratic derivative" << endl;
    }
    else {
      error("INVALID TYPE IN runFusedGeneral\n");
    }
   if(quadDer -> isExtreme()) {
       return(false);
   }

   // cout << "Return value:" << (result && outerIterNum < maxIterOuter) << endl;
   return (result && outerIterNum < maxIterOuter);
}

bool FusedLasso::runAlgorithmL2() {
    outerIterNum = 0;
    innerIterNum = 0;
    bool lastRunOK = true;
    // set single fusions
    fusedGroupSize.resize(p, 1);
    for(int i = 0; i < p; ++i) {
        fusions[i] = i;
    }

    if(regType != GAUSSIAN) {
        lastRunOK = runFusedGeneral(L2);
    }
    else {
        lastRunOK = runFused(L2);
    }
    
    if(DEBUG) {
      REprintf("Inner iter num: %d\n", innerIterNum);
    }

    if(lastRunOK) {
        return true;
    }
    else {
        return false;
    }
}


bool FusedLasso::runAlgorithmHuber() {
    
    outerIterNum = 0;
    innerIterNum = 0;
    bool lastRunOK = true;
    double lastIterChange;

    vector<double> oldBeta;

    if(regType != GAUSSIAN) {
        lastRunOK = runFusedGeneral(L1);
    }
    else {
        lastRunOK = runFused(L1);
    }
    // check if we have an extreme result now
    if(quadratic->isExtreme()) {
        return false;
    }

//    while(outerIterNum < maxIterOuter && lastRunOK) {
        oldBeta = beta;

        //newFusionDifferent = identifyNewFusions(lastIterChange);
        if(DEBUG) {
	  REprintf("Before identifyNewFusionsHuber\n");
        }
        bool newFusionDifferent = identifyNewFusionsHuber();
        if(DEBUG) {
	  REprintf("After identifyNewFusionsHuber\n");
        }

        vector<int> newFusions = fusions;
        vector<int> newFusedGroupSize = fusedGroupSize;

        do { 
            oldBeta = beta;
            if(DEBUG) {
	      REprintf("Iteration number %d\n", outerIterNum);
            }
            fusions = newFusions;
            fusedGroupSize = newFusedGroupSize;
            if(regType != GAUSSIAN) {
	      if(DEBUG) {
                REprintf("Before runFusedGeneral(L1)\n");
	      }

                lastRunOK = runFusedGeneral(L1);
		if(DEBUG) {
		  REprintf("After runFusedGeneral(L1)\n");
		}
            }
            else {
                lastRunOK = runFused(L1);

            }
	    //            cout << "Beta" << endl;
	    //            printVector(beta, cout);
	    //            cout << "Fusions" << endl;
	    //            printVector(fusions, cout);

            // check if we have an extreme result now
            if(quadratic->isExtreme()) {
                return false;
            }
            getEqualFusions(newFusions, newFusedGroupSize);

            fusions = newFusions;
            fusedGroupSize = newFusedGroupSize;
	    //            cout << "Fusions " << endl;
	    //            printVector(fusions, cout);
	    //            printVector(beta, cout);
            if(DEBUG && outerIterNum > 100) {
	      /* cout << newFusedGroupSize.size() << endl;
                for(int i = 0; i < fusions.size(); ++i) {
                    if(fusions[i] != newFusions[i]) {
		       cout << i << ":" << fusions[i] << " " << newFusions[i] << "; ";
			 cout << "Beta: " << beta[i] << "Oldbeta: " << oldBeta[i] << "; ";  
                    }
                }
                cout << endl; */
            }
            lastIterChange = maxDiffDoubleVec(oldBeta, beta);
        } while(!areFusionsEqual(fusions, newFusions) && lastIterChange > accuracy);

//        if(lastIterChange < accuracy) {
//            break;
//        }
//    }
    
    if(DEBUG) {
      REprintf("Number of outer iterations: %d\n", outerIterNum);
    }
    if(lastRunOK && outerIterNum < maxIterOuter) {
        return true;
    }
    else {
        return false;
    }
}

bool FusedLasso::runAlgorithm(int maxFusionLevel) {
    outerIterNum = 0;
    innerIterNum = 0;
    bool lastRunOK = true;
    bool newFusionDifferent;
    vector<double> oldBeta; 
    fusionLevel = 0;

    while(outerIterNum < maxIterOuter && lastRunOK) {
        if(DEBUG) {
	  REprintf("Outer iterNum %d\n", outerIterNum);
        }

        oldBeta = beta;
        if(regType != GAUSSIAN) {
            lastRunOK = runFusedGeneral(L1);
        }
        else {
            lastRunOK = runFused(L1);
        }

        // check if we have an extreme result now
        if(quadratic->isExtreme()) {
            lastRunOK = false;
            break;
        }

        double lastIterChange = maxDiffDoubleVec(oldBeta, beta);
        //newFusionDifferent = identifyNewFusions(lastIterChange);
        newFusionDifferent = identifyNewFusions(lastIterChange, maxFusionLevel);
       
        //cout << "In runAlgorithm" << endl;
        // printVector(beta, cout);
        //REprintf("Number of outer iterations: %d\n", outerIterNum);

        if(!newFusionDifferent) { // something has changed
            break; // algorithm has converged
        }
    }
    
    if(DEBUG) {
      REprintf("Number of outer iterations: %d\n", outerIterNum);
    }


    if(lastRunOK && outerIterNum < maxIterOuter) {
        return true;
    }
    else {
        return false;
    }

    return false;
}


SparseMatrix FusedLasso::runAlgorithm(penEnum penType, int maxFusionLevel, const vector<double>& lambda1Vec, const vector<double>& lambda2Vec, const int maxNonZero, vector<bool>& success, vector<int>& outerIterNumVec, vector<int>& innerIterNumVec, bool verbose) {
    SparseMatrix betaSols(beta.size());
    success.clear(); success.resize(lambda1Vec.size());
    outerIterNumVec.clear(); outerIterNumVec.resize(lambda1Vec.size());
    innerIterNumVec.clear(); innerIterNumVec.resize(lambda1Vec.size());

    for(int i = 0; i < lambda1Vec.size(); ++i) {
        if(verbose || DEBUG) {
	  REprintf("Working on Lambda1:i: %d Lambda1: %f Lambda2: %f\n", lambda2Vec[i], i, lambda1Vec[i], lambda2Vec[i]);
        }    
        setNewLambdas(lambda1Vec[i], lambda2Vec[i]);
        switch(penType) {
            case L1: if(maxFusionLevel >= -1) {
                        success[i] = runAlgorithm(maxFusionLevel);        
                     }
                     else {
                        success[i] = runAlgorithmHuber();
                     }
                     break;
            case L2: success[i] = runAlgorithmL2(); break;
        }
        outerIterNumVec[i] = getOuterIterNum();
        innerIterNumVec[i] = getInnerIterNum();
        betaSols.addColumn(getBeta());
        if(verbose || DEBUG) {
	  REprintf("Number active beta: %d\n", numNonZero(getBeta()));
        }
        if(numNonZero(getBeta()) > maxNonZero) {
            break; // number really saved can be read from the returned matrix
        }
        if(!success[i]) {
            break; // algorithm did not finish successfully
        }
    }
    return betaSols;
}

int FusedLasso::getOuterIterNum() {
    return outerIterNum;
}

int FusedLasso::getInnerIterNum() {
    return innerIterNum;
}

void FusedLasso::setNewLambdas(double lambda1, double lambda2) {
    // when setting the new lambdas, only the lambdas
    // have to be reset; no other changes are necessary
    this->lambda1 = lambda1;
    this->lambda2 = lambda2;
}


double FusedLasso::findMaxLambda1(const list<int>& exemptVars) {
    if(DEBUG) {
      REprintf("Beginning of findMaxLambda1");
    }
    vector<vector<int> > connectionsEmpty(p,vector<int>(0));
    vector<vector<double> > wLambda2Empty(p, vector<double>(0));
    vector<double> wLambda1Extreme(p, 1e6);
    list<int>::const_iterator listIt;
    for(listIt = exemptVars.begin(); listIt != exemptVars.end(); ++listIt) {
        wLambda1Extreme[*listIt] = 1e-4;
    }

    FusedLassoCoordinate flc(quadratic, wLambda1Extreme, connectionsEmpty, wLambda2Empty, maxIterInner, accuracy, maxActivateVars, 1, 1, L1);

    bool result = flc.runAlgorithm();

    this->beta = flc.getBeta();

    switch(regType) {
        case BINOMIAL:
            quadratic = counted_ptr<QuadraticDerivative>(new QuadraticDerivativeLogistic(this->X, this->y, this->wObs, this->beta));
    }

    double highLambda1 = 0;
    double foo;
    // now search through the 0 variables which will start next
    for(int i = 0; i < p; ++i) {
        if(quadratic->getBeta(i) == 0) {
            foo = fabs(quadratic->getDerivative(i) / wLambda1[i]);
            if(highLambda1 < foo) {
                highLambda1 = foo;
            }
        }
    }
    if(DEBUG) {
      REprintf("End of findMaxLambda1");
    }
    return highLambda1;
}


// identify lambda2
double FusedLasso::findMaxLambda2(double lambda1_here) {
    if(p > 100000) {
      REprintf("Too big to estimate lambda2; returning 1");
        return n;
    }

    if(DEBUG) {
      REprintf("Beginning of findMaxLambda2");
    }
    outerIterNum=0; // reset as this is run several times
    // first, set beta equal to 0.1
    for(int i = 0; i < p; i++) {
        beta[i] = 0.1;
        quadratic->updateBeta(i, 0.1);
    }

    vector<int> newFusions;
    vector<int> newFusedGroupSize;
    
    // fuse all variables and run them to equilibrium
    lambda1 = lambda1_here;
    lambda2 = n * 2e6;
    setNewLambdas(lambda1, lambda2);
    getEqualFusions(newFusions, newFusedGroupSize, false, 0.001/accuracy);
    //cout << "===============After get equalFusions" << endl;
    //cout << "===============before first runAlgorithmHuber" << endl;
    runAlgorithmHuber();
    //cout << "===============After first runAlgorithmHuber" << endl;

    vector<double> oldBeta(p,0);
    for(int i=0; i < p; ++i) {
      oldBeta[i] = beta[i];
    }
    bool hasChanged;
    hasChanged=false;
    //    printVector(oldBeta,cout);
    //    cout <<  "Now start" << endl;

//    cout << "numFusedGroups " << numFusedGroups << endl;
    lambda2/=2;
//    cout << "Beta" << endl;
//    printVector(beta, cout);

    for(int i = 0; i < 30; ++i) {
        if(DEBUG) {
	  //cout << "Lambda1: " << lambda1 << " Lambda2 " << lambda2 << "Accuracy" << accuracy << " i: " << i << endl;
        }
        setNewLambdas(lambda1, lambda2);
        runAlgorithmHuber();
	//        getEqualFusions(newFusions, newFusedGroupSize, false, 0.001/accuracy);
	for(int i = 0; i < p; ++i) {
	  if(fabs(oldBeta[i] - beta[i]) > accuracy * 10) {
	    hasChanged = true;
	  }
	}
	//        printVector(beta, cout);
//        cout << "FusedGroupSize " << newFusedGroupSize.size() << endl;
        // check that there is only group, otherwise, last step was maximal
	//        if(newFusedGroupSize.size() > numFusedGroups) {
	if(hasChanged) {
            return(lambda2*2);
        } 
        else {
            lambda2 /= 2;
        }
    } 
    return(lambda2);
}


void FusedLasso::getBeta(SparseMatrix& betaMat) {
    betaMat.addColumn(beta);
}


void FusedLasso::printBetaAndGroup(ostream& outStream) {
     outStream << "========== Fusions =================" << endl;

    int curGroup = -1;
    for(int i = 0; i < p; ++i) {
        if(fusions[i] > curGroup) { // a new group has started
            curGroup++;
            outStream << " Group: " << i << " Beta: " << beta[i] ;
        }
    }

    outStream << endl << "========== End Fusions =============" << endl;
}


void FusedLasso::makePullAdjustment(const vector<double>& beta, vector<double>& nodePull, double lambda2) {
    pg.makePullAdjustment(beta, nodePull, lambda2, accuracy);
}

double FusedLasso::calcGroupAverage(int pos, vector<double>& nodePull, vector<double>& beta) {
    double curBeta = beta[pos];
    double average = 0;
    double groupSize = 0;
    while(fabs(beta[pos] - curBeta) < accuracy) {
        groupSize++;
        average += nodePull[pos];
        ++pos;
    }
    return average/groupSize;
} 


void FusedLasso::printGroups(vector<list<int> > x, ostream& outStream) {
    list<int>::iterator listIt;
    for(int i = 0; i < x.size(); ++i) {
        for(listIt = x[i].begin(); listIt != x[i].end(); ++listIt) {
            outStream << "i: " << i << " : " << *listIt << endl;
        }
    }
}





void FusedLasso::checkSolution(ostream& outStream) {
    vector<double> nodePull = getPulls();
    vector<double> nodePullOriginal = nodePull;
    vector<double> nodePullAdj = nodePull;
    makePullAdjustment( beta, nodePullAdj, lambda2);
    vector<int> newFusions;
    vector<int> newFusedGroupSize;

    // make all groups that are equal
    getEqualFusions(newFusions, newFusedGroupSize);

    // first adjust those for lambda1
    int curPos = 0;
    int curGroup = newFusions[curPos];
    double adjustment;

    outStream << "Lambda1: " << lambda1 << endl;
    outStream << "Lambda2: " << lambda2 << endl;

    while(curPos < nodePull.size()) {
    
        if(beta[curPos] > 0) {
            nodePull[curPos] += lambda1;
            curPos++;
        }
        else if(beta[curPos] < 0) {
            nodePull[curPos] -= lambda1;
            curPos++;
        }
        else { // is 0, so just subtract the average of the group
            curGroup = newFusions[curPos];
            adjustment = calcGroupAverage(curPos, nodePullAdj, beta);
            if(fabs(adjustment) > lambda1) {
                outStream << "=============== Problem =================" << endl;
                outStream << "Group " << curGroup << " at position " << curPos <<
                    " has adjustment " << adjustment << endl;
            }
            while(newFusions[curPos] == curGroup && curPos < nodePull.size()) {
                nodePull[curPos] -= adjustment;
                curPos++;
            }
        }
    }

    // now run the maximum flow graph and check that 
    list<int> allNodes = pg.allNodes();
    // now adjust nodePull by dividing by lambda2
    for(int i = 0; i < nodePull.size(); ++i) {
        nodePull[i] /= lambda2;
        outStream << i << ": " << nodePull[i] << endl;
    }

    pg.initializeMaxFlow(allNodes, nodePull);
    pg.findMaxFlowEdmondsKarp();
    GraphEdge* e, *eBack;

    // now print out a compact version of the solution
    for(int i = 0; i < beta.size(); ++i) {
        e = pg.nodes[i].back();
        eBack = e->backwards;
        outStream << "Node " << i << " Beta: " << beta[i] << " Deriv: " << nodePullOriginal[i] <<
            " DerivAdj: " << nodePull[i];
        if(i < beta.size()-1) {
            outStream << "t[i,i+1]: " << e->flow << "Back: " << eBack->flow;
        }
        outStream << endl;
    }    
   
    outStream << "================= Whole MaxflowGraph ==================" << endl;
    //    pg.printGraph(cout);
    outStream << "==================== END WHOLE =======================" << endl;
}







