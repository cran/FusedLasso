#include "FusedLasso.h"
#include "SparseMatrix.h"
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <iostream>
#include <string.h>

using namespace std;


extern "C" {
    SEXP copySparseMatrix(SEXP from) {
        SparseMatrix X(from);

        return X.todgCMatrix();
    }

SEXP FusedLassoWrapperL2(SEXP R_XMat, SEXP R_yVec, SEXP R_wObs, SEXP R_beta, SEXP R_wLambda1, SEXP R_graph, SEXP R_maxIterInner, SEXP R_maxIterOuter, SEXP R_accuracy, SEXP R_maxActivateVars, SEXP R_lambda1Vec, SEXP R_lambda2Vec, SEXP R_maxNonZero, SEXP R_addNodes, SEXP R_family, SEXP R_verbose) {
        // first, create a sparse matrix for R_XMat
        SparseMatrix X(R_XMat);                

        // vectors y, wObs, beta
        vector<double> y(REAL(R_yVec), REAL(R_yVec) + LENGTH(R_yVec));
        vector<double> wObs(REAL(R_wObs), REAL(R_wObs) + LENGTH(R_wObs));
        vector<double> beta(REAL(R_beta), REAL(R_beta) + LENGTH(R_beta));
        vector<double> wLambda1(REAL(R_wLambda1), REAL(R_wLambda1) + LENGTH(R_wLambda1));

        // the parameters
        int maxIterInner = INTEGER(R_maxIterInner)[0];
        int maxIterOuter = INTEGER(R_maxIterOuter)[0];
        double accuracy = REAL(R_accuracy)[0];
        int maxActivateVars = INTEGER(R_maxActivateVars)[0];
        int maxNonZero = INTEGER(R_maxNonZero)[0];
        bool verbose = LOGICAL(R_verbose)[0];

        // the graph
        Graph graph;
        if(Rf_isInteger(R_graph) && LENGTH(R_graph) == 1) {
            graph = Graph(INTEGER(R_graph)[0]);
        }
        else if(Rf_isInteger(R_graph) && LENGTH(R_graph) == 2) {
            graph = Graph(INTEGER(R_graph)[0], INTEGER(R_graph)[1]);
        }
        else if(Rf_isInteger(R_graph) && LENGTH(R_graph) == 3) {
            graph = Graph(INTEGER(R_graph)[0], INTEGER(R_graph)[1], INTEGER(R_graph)[2]);
        }
        else { // assumes that it is a list; should check before calling FusedLassoWrapper
            graph = Graph(VECTOR_ELT(R_graph, 0), VECTOR_ELT(R_graph, 1));
        }
        int addNodes = INTEGER(R_addNodes)[0];
        for(int i = 0; i < addNodes; ++i) {
            graph.addNode();
        }

        // the type of the regression
        regEnum regType;
        if(strcmp(CHAR(STRING_ELT(R_family,0)), "binomial")==0) {
            regType = BINOMIAL;
        }
        else if(strcmp(CHAR(STRING_ELT(R_family,0)), "cox")==0) {
            regType = COX;
        }
        else { // default to gaussian
            regType = GAUSSIAN;
        }

        // the lambdas
        vector<double> lambda1Vec(REAL(R_lambda1Vec), REAL(R_lambda1Vec) + LENGTH(R_lambda1Vec));
        vector<double> lambda2Vec(REAL(R_lambda2Vec), REAL(R_lambda2Vec) + LENGTH(R_lambda2Vec));


        // vector for storing if converged or not
        vector<bool> success;
        vector<int> outerIterNum;
        vector<int> innerIterNum;

        FusedLasso fl(X, y, wObs, beta, wLambda1, graph, maxIterInner, maxIterOuter, accuracy, maxActivateVars, 0, 0, regType);
        // adjust the lambda1 for maximum value if necessary

        SparseMatrix res = fl.runAlgorithm(L2, 0, lambda1Vec, lambda2Vec, maxNonZero, success, outerIterNum, innerIterNum, verbose);

        SEXP R_result, R_success, R_names, R_outerIterNum, R_innerIterNum,
             R_lambda1VecNew, R_lambda2VecNew;

        int numRuns = res.p; // number of lambdas performed before stop
        PROTECT(R_success = allocVector(LGLSXP, numRuns));
        PROTECT(R_outerIterNum = allocVector(INTSXP, numRuns));
        PROTECT(R_innerIterNum = allocVector(INTSXP, numRuns));
        PROTECT(R_lambda1VecNew = allocVector(REALSXP, numRuns));
        PROTECT(R_lambda2VecNew = allocVector(REALSXP, numRuns));
        for(int i = 0; i < numRuns; ++i) {
            LOGICAL(R_success)[i] = success[i];
            INTEGER(R_outerIterNum)[i] = outerIterNum[i];
            INTEGER(R_innerIterNum)[i] = innerIterNum[i];
            REAL(R_lambda1VecNew)[i] = lambda1Vec[i];
            REAL(R_lambda2VecNew)[i] = lambda2Vec[i];
        }
        PROTECT(R_result = allocVector(VECSXP, 7));
        SET_VECTOR_ELT(R_result, 0, res.todgCMatrix());
        SET_VECTOR_ELT(R_result, 1, R_success);
        SET_VECTOR_ELT(R_result, 2, R_lambda1VecNew);
        SET_VECTOR_ELT(R_result, 3, R_lambda2VecNew);
        SET_VECTOR_ELT(R_result, 4, R_outerIterNum);
        SET_VECTOR_ELT(R_result, 5, R_innerIterNum);
        SET_VECTOR_ELT(R_result, 6, R_family);
        
        PROTECT(R_names = allocVector(STRSXP, 7));
        SET_STRING_ELT(R_names, 0, mkChar("beta"));
        SET_STRING_ELT(R_names, 1, mkChar("success"));
        SET_STRING_ELT(R_names, 2, mkChar("lambda1"));
        SET_STRING_ELT(R_names, 3, mkChar("lambda2"));
        SET_STRING_ELT(R_names, 4, mkChar("outerIterNum"));
        SET_STRING_ELT(R_names, 5, mkChar("innerIterNum"));
        SET_STRING_ELT(R_names, 6, mkChar("family"));
        setAttrib(R_result, R_NamesSymbol, R_names);
        UNPROTECT(7);

        return R_result;
    }


SEXP FusedLassoWrapperL1(SEXP R_XMat, SEXP R_yVec, SEXP R_wObs, SEXP R_beta, SEXP R_wLambda1, SEXP R_graph, SEXP R_maxIterInner, SEXP R_maxIterOuter, SEXP R_accuracy, SEXP R_maxActivateVars, SEXP R_lambda1Vec, SEXP R_lambda2Vec, SEXP R_maxNonZero, SEXP R_addNodes, SEXP R_family, SEXP R_adjustWithMax, SEXP R_maxFusionLevel, SEXP R_verbose) {
        // first, create a sparse matrix for R_XMat
        SparseMatrix X(R_XMat);                

        // vectors y, wObs, beta
        vector<double> y(REAL(R_yVec), REAL(R_yVec) + LENGTH(R_yVec));
        vector<double> wObs(REAL(R_wObs), REAL(R_wObs) + LENGTH(R_wObs));
        vector<double> beta(REAL(R_beta), REAL(R_beta) + LENGTH(R_beta));
        vector<double> wLambda1(REAL(R_wLambda1), REAL(R_wLambda1) + LENGTH(R_wLambda1));

        // the parameters
        int maxIterInner = INTEGER(R_maxIterInner)[0];
        int maxIterOuter = INTEGER(R_maxIterOuter)[0];
        double accuracy = REAL(R_accuracy)[0];
        int maxActivateVars = INTEGER(R_maxActivateVars)[0];
        int maxNonZero = INTEGER(R_maxNonZero)[0];
        int maxFusionLevel = INTEGER(R_maxFusionLevel)[0];
        bool verbose = LOGICAL(R_verbose)[0];

        // the graph
        Graph graph;
        if(Rf_isInteger(R_graph) && LENGTH(R_graph) == 1) {
            graph = Graph(INTEGER(R_graph)[0]);
        }
        else if(Rf_isInteger(R_graph) && LENGTH(R_graph) == 2) {
            graph = Graph(INTEGER(R_graph)[0], INTEGER(R_graph)[1]);
        }
        else if(Rf_isInteger(R_graph) && LENGTH(R_graph) == 3) {
            graph = Graph(INTEGER(R_graph)[0], INTEGER(R_graph)[1], INTEGER(R_graph)[2]);
        }
        else { // assumes that it is a list; should check before calling FusedLassoWrapper
            graph = Graph(VECTOR_ELT(R_graph, 0), VECTOR_ELT(R_graph, 1));
        }
        int addNodes = INTEGER(R_addNodes)[0];
        for(int i = 0; i < addNodes; ++i) {
            graph.addNode();
        }

        // the type of the regression
        regEnum regType;
        if(strcmp(CHAR(STRING_ELT(R_family,0)), "binomial")==0) {
            regType = BINOMIAL;
        }
        else if(strcmp(CHAR(STRING_ELT(R_family,0)), "cox")==0) {
            regType = COX;
        }
        else { // default to gaussian
            regType = GAUSSIAN;
        }

        // the lambdas
        vector<double> lambda1Vec(REAL(R_lambda1Vec), REAL(R_lambda1Vec) + LENGTH(R_lambda1Vec));
        vector<double> lambda2Vec(REAL(R_lambda2Vec), REAL(R_lambda2Vec) + LENGTH(R_lambda2Vec));


        // vector for storing if converged or not
        vector<bool> success;
        vector<int> outerIterNum;
        vector<int> innerIterNum;
	
	//REprintf("Before Fused lasso\n");

        FusedLasso fl(X, y, wObs, beta, wLambda1, graph, maxIterInner, maxIterOuter, accuracy, maxActivateVars, 0, 0, regType);
        // adjust the lambda1 for maximum value if necessary
        if(LOGICAL(R_adjustWithMax)[0]==true) {
            list<int> exemptVars;
            for(int i = 0; i < wLambda1.size(); ++i) {
                if(wLambda1[i] < 1e-4) {
                    exemptVars.push_back(i);   
                }
            }   
            double maxLambda1 = fl.findMaxLambda1(exemptVars);
            for(int i = 0; i < lambda1Vec.size(); ++i) {
                lambda1Vec[i] *= maxLambda1;
            }
        }


	//REprintf("Before runAlgorithm\n");
        SparseMatrix res = fl.runAlgorithm(L1, maxFusionLevel, lambda1Vec, lambda2Vec, maxNonZero, success, outerIterNum, innerIterNum, verbose);
	//REprintf("After runAlgorithm\n");

        SEXP R_result, R_success, R_names, R_outerIterNum, R_innerIterNum,
             R_lambda1VecNew, R_lambda2VecNew;

	//REprintf("Before copying of results\n");
        int numRuns = res.p; // number of betas performed before stop
        PROTECT(R_success = allocVector(LGLSXP, numRuns));
        PROTECT(R_outerIterNum = allocVector(INTSXP, numRuns));
        PROTECT(R_innerIterNum = allocVector(INTSXP, numRuns));
        PROTECT(R_lambda1VecNew = allocVector(REALSXP, numRuns));
        PROTECT(R_lambda2VecNew = allocVector(REALSXP, numRuns));
        for(int i = 0; i < numRuns; ++i) {
            LOGICAL(R_success)[i] = success[i];
            INTEGER(R_outerIterNum)[i] = outerIterNum[i];
            INTEGER(R_innerIterNum)[i] = innerIterNum[i];
            REAL(R_lambda1VecNew)[i] = lambda1Vec[i];
            REAL(R_lambda2VecNew)[i] = lambda2Vec[i];
        }
        PROTECT(R_result = allocVector(VECSXP, 7));
        SET_VECTOR_ELT(R_result, 0, res.todgCMatrix());
        SET_VECTOR_ELT(R_result, 1, R_success);
        SET_VECTOR_ELT(R_result, 2, R_lambda1VecNew);
        SET_VECTOR_ELT(R_result, 3, R_lambda2VecNew);
        SET_VECTOR_ELT(R_result, 4, R_outerIterNum);
        SET_VECTOR_ELT(R_result, 5, R_innerIterNum);
        SET_VECTOR_ELT(R_result, 6, R_family);
        
        PROTECT(R_names = allocVector(STRSXP, 7));
        SET_STRING_ELT(R_names, 0, mkChar("beta"));
        SET_STRING_ELT(R_names, 1, mkChar("success"));
        SET_STRING_ELT(R_names, 2, mkChar("lambda1"));
        SET_STRING_ELT(R_names, 3, mkChar("lambda2"));
        SET_STRING_ELT(R_names, 4, mkChar("outerIterNum"));
        SET_STRING_ELT(R_names, 5, mkChar("innerIterNum"));
        SET_STRING_ELT(R_names, 6, mkChar("family"));
        setAttrib(R_result, R_NamesSymbol, R_names);
        UNPROTECT(7);

        return R_result;
    }

    
    SEXP FusedLassoWrapperMaxLambdas(SEXP R_XMat, SEXP R_yVec, SEXP R_wObs, SEXP R_wLambda1, SEXP R_graph, SEXP R_maxIterInner, SEXP R_accuracy, SEXP R_maxActivateVars, SEXP R_exemptVars, SEXP R_addNodes, SEXP R_family) {
        // first, create a sparse matrix for R_XMat
        SparseMatrix X(R_XMat);                

        // vectors y, wObs, beta
        vector<double> y(REAL(R_yVec), REAL(R_yVec) + LENGTH(R_yVec));
        vector<double> wObs(REAL(R_wObs), REAL(R_wObs) + LENGTH(R_wObs));
        vector<double> beta(LENGTH(R_wLambda1), 0);
        vector<double> wLambda1(REAL(R_wLambda1), REAL(R_wLambda1) + LENGTH(R_wLambda1));

        // the parameters
        int maxIterInner = INTEGER(R_maxIterInner)[0];
        int maxIterOuter = 100; // actually not needed
        double accuracy = REAL(R_accuracy)[0];
        int maxActivateVars = INTEGER(R_maxActivateVars)[0];

        // the graph
        Graph graph;
        if(Rf_isInteger(R_graph) && LENGTH(R_graph) == 1) {
            graph = Graph(INTEGER(R_graph)[0]);
        }
        else if(Rf_isInteger(R_graph) && LENGTH(R_graph) == 2) {
            graph = Graph(INTEGER(R_graph)[0], INTEGER(R_graph)[1]);
        }
        else if(Rf_isInteger(R_graph) && LENGTH(R_graph) == 3) {
            graph = Graph(INTEGER(R_graph)[0], INTEGER(R_graph)[1], INTEGER(R_graph)[2]);
//            cout << "Graph: " << INTEGER(R_graph)[0] << INTEGER(R_graph)[1] << INTEGER(R_graph)[2] << endl;
        }
        else { // assumes that it is a list; should check before calling FusedLassoWrapper
            graph = Graph(VECTOR_ELT(R_graph, 0), VECTOR_ELT(R_graph, 1));
        }
        int addNodes = INTEGER(R_addNodes)[0];
        for(int i = 0; i < addNodes; ++i) {
            graph.addNode();
        }
	
        // the type of the regression
        regEnum regType;
        if(strcmp(CHAR(STRING_ELT(R_family,0)), "binomial")==0) {
            regType = BINOMIAL;
        }
        else if(strcmp(CHAR(STRING_ELT(R_family,0)), "cox")==0) {
            regType = COX;
        }
        else { // default to gaussian
            regType = GAUSSIAN;
        }

        FusedLasso fl(X, y, wObs, beta, wLambda1, graph, maxIterInner, maxIterOuter, accuracy, maxActivateVars, 0, 0, regType);

        list<int> exemptVars;
        for(int i = 0; i < LENGTH(R_exemptVars); ++i) {
            exemptVars.push_back(INTEGER(R_exemptVars)[i]-1);
        }
        double maxLambda1 = fl.findMaxLambda1(exemptVars);

	// run for several values of lambda1
	double lambda1_here = maxLambda1;
	double maxLambda2= 1e-5;
	double maxLambda2_here;
	for(int i=0; i<5; i++) {
	  lambda1_here /= 4;
	  maxLambda2_here = fl.findMaxLambda2(lambda1_here);
	  // cout << "MaxLambda2_here: " << maxLambda2_here << endl;
	  if(maxLambda2_here > maxLambda2) {
	    maxLambda2 = maxLambda2_here;
	  }
	}
        SEXP R_result, R_names, R_maxLambda1, R_maxLambda2;
        PROTECT(R_maxLambda1 = allocVector(REALSXP, 1));
        REAL(R_maxLambda1)[0] = maxLambda1;
        PROTECT(R_maxLambda2 = allocVector(REALSXP, 1));
        REAL(R_maxLambda2)[0] = maxLambda2;
        PROTECT(R_result = allocVector(VECSXP, 2));
        SET_VECTOR_ELT(R_result, 0, R_maxLambda1);
        SET_VECTOR_ELT(R_result, 1, R_maxLambda2);
        
        PROTECT(R_names = allocVector(STRSXP, 2));
        SET_STRING_ELT(R_names, 0, mkChar("maxLambda1"));
        SET_STRING_ELT(R_names, 1, mkChar("maxLambda2"));
        setAttrib(R_result, R_NamesSymbol, R_names);
        UNPROTECT(4);

        return R_result;
    }
}

