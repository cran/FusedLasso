#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

#include <vector>
#include <iostream>
#include <Rinternals.h>
#include <R.h>
#include <Rdefines.h>

using namespace std;


// Very simple implementation of a sparse matrix
// In the current version, the matrix sparseness structure can not be changed after construction
// It is only usd for the X matrix, which does not need to be changed 

class SparseMatrix {
public:
    int nnz;
    int n;
    int p;
    vector<double> nzVec; // vector of non-zero elements
    vector<int> indexCol; // index of the first element of x of column; length p+1
    vector<int> rowPos; // row position of the ith element of x; length size(x)

public:

    // Constructors
    SparseMatrix();
    // also empty, but with specified number of rows
    SparseMatrix(int n);
    // From a c-Array
    SparseMatrix(const double*x, int n, int p);
    //From a vector
    SparseMatrix(const vector<double>& x, int n, int p);
    // From an already sparse format
    SparseMatrix(const vector<double>& x, const vector<int>& indexCol, const vector<int>& rowPos, int n, int p);
    // create from a dgCMatrix; assumes that it is the right class and does not check
    SparseMatrix(SEXP x);
    // given a set of fusions, generate a new matrix (could be extended to multiplying
    // two matrices
    SparseMatrix createFusedX(const vector<vector<int> > fusedGroups);

    // add a column to the matrix
    // in sparse format
    void addColumn(const vector<double>& newnzVec, const vector<int>& newrowPos);
    // in non-sparse format
    void addColumn(const vector<double>& x);

    // calculate inner product of ith column with itself
    double innerProd(int i) const;
    // calculate inner product of column i and column j
    double innerProd(int i, int j) const;

    // inner product of column i with vector y
    double innerProd(const vector<double>& y, int i) const;

    // WEIGHTED VERSION OF INNER PRODUCTS 
    // calculate inner product of ith column with itself
    double innerProd(int i, const vector<double>& w) const;
    // calculate inner product of column i and column j
    double innerProd(int i, int j, const vector<double>& w) const;

    // inner product of column i with vector y
    double innerProd(const vector<double>& y, int i, const vector<double>& w) const;

   // multiply the sparse matrix by a diagonal matrix from the left
    void multByDiag(const vector<double>& w);

    // add a multiple of column i to the given vector y; done in place
    void addMultOfColumn(vector<double>& y, int i, double mult) const;

    // return the value of a position in the matrix
    // with row i and column j
    double get(int i, int j) const;

    // return the i-th column as a vector
    vector<double> getColumn(int i) const;

    SEXP todgCMatrix() const;

    // print out the current object (used for debugging purposes)
    void print(ostream& outStream) const;

    // print column
    void printColumn(int i, ostream& outStream) const;
};

#endif // SPARSE_MATRIX_H_
