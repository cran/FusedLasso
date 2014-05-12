#include "Lasso.h"
#include "SparseMatrix.h"
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <iostream>

using namespace std;


extern "C" {

    SEXP LassoWrapper(SEXP R_XMat, SEXP R_yVec, SEXP R_wObs, SEXP R_beta, SEXP R_wLambda1, SEXP R_maxIterInner, SEXP R_accuracy, SEXP R_maxActivateVars, SEXP R_lambda1Vec) {
        // first, create a sparse matrix for R_XMat
        SparseMatrix X(R_XMat);                

        // vectors y, wObs, beta
        vector<double> y(REAL(R_yVec), REAL(R_yVec) + LENGTH(R_yVec));
        vector<double> wObs(REAL(R_wObs), REAL(R_wObs) + LENGTH(R_wObs));
        vector<double> beta(REAL(R_beta), REAL(R_beta) + LENGTH(R_beta));
        vector<double> wLambda1(REAL(R_wLambda1), REAL(R_wLambda1) + LENGTH(R_wLambda1));

        // the parameters
        int maxIterInner = INTEGER(R_maxIterInner)[0];
        double accuracy = REAL(R_accuracy)[0];
        int maxActivateVars = INTEGER(R_maxActivateVars)[0];

        // the lambdas
        vector<double> lambda1Vec(REAL(R_lambda1Vec), REAL(R_lambda1Vec) + LENGTH(R_lambda1Vec));

        // vector for storing if converged or not
        vector<bool> success;

        Lasso fl(X, y, wObs, beta, wLambda1, maxIterInner, accuracy, maxActivateVars, 0);
        SparseMatrix res = fl.runAlgorithm(lambda1Vec, success);

        SEXP R_result, R_success, R_names;
        PROTECT(R_success = allocVector(LGLSXP, success.size()));
        for(int i = 0; i < success.size(); ++i) {
            LOGICAL(R_success)[i] = success[i];
        }
        PROTECT(R_result = allocVector(VECSXP, 3));
        SET_VECTOR_ELT(R_result, 0, res.todgCMatrix());
        SET_VECTOR_ELT(R_result, 1, R_success);
        SET_VECTOR_ELT(R_result, 2, R_lambda1Vec);
        
        PROTECT(R_names = allocVector(STRSXP, 3));
        SET_STRING_ELT(R_names, 0, mkChar("beta"));
        SET_STRING_ELT(R_names, 1, mkChar("success"));
        SET_STRING_ELT(R_names, 2, mkChar("lambda1"));
        setAttrib(R_result, R_NamesSymbol, R_names);
        UNPROTECT(3);

        return R_result;
    }

}

