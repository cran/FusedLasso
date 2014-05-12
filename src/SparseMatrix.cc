#include "SparseMatrix.h"
#include "GeneralFunctions.h"
#include <stdlib.h>
#include <R.h>


SparseMatrix::SparseMatrix() {
    p = 0;
    n = 0;
    nnz = 0;
    nzVec.clear();
    indexCol.clear(); indexCol.push_back(0);
    rowPos.clear();
}


SparseMatrix::SparseMatrix(int n) {
    p = 0;
    this->n = n;
    nnz = 0;
    nzVec.clear();
    indexCol.clear(); indexCol.push_back(0);
    rowPos.clear();
}




SparseMatrix::SparseMatrix(const double* x, int n, int p) {
    this->p = p;
    this->n = n;
    this->nnz = 0;
    bool foundFirst;

    for(int i = 0; i < p; ++i) {
        foundFirst = false;
        for(int j = 0; j < n; ++j) {
            if(x[i*n + j] != 0) {
                nnz++;
                nzVec.push_back(x[i*n + j]);
                rowPos.push_back(j);
                if(!foundFirst) {
                    foundFirst = true;
                    indexCol.push_back(nnz-1);
                }
            }
        }
        if(!foundFirst) { // column is empty
            indexCol.push_back(nnz);
        }
    }
    indexCol.push_back(nnz);
}

SparseMatrix::SparseMatrix(const vector<double>& x, int n, int p) {
    this->p = p;
    this->n = n;
    this->nnz = 0;
    bool foundFirst;

    for(int i = 0; i < p; ++i) {
        foundFirst = false;
        for(int j = 0; j < n; ++j) {
            if(x[i*n + j] != 0) {
                nnz++;
                nzVec.push_back(x[i*n + j]);
                rowPos.push_back(j);
                if(!foundFirst) {
                    foundFirst = true;
                    indexCol.push_back(nnz-1);
                }
            }
        }
        if(!foundFirst) { // column is empty
            indexCol.push_back(nnz);
        }
    }
    indexCol.push_back(nnz);
}

SparseMatrix::SparseMatrix(const vector<double>& x, const vector<int>& indexCol, const vector<int>& rowPos, int n, int p) {
    nnz = x.size();
    this-> p = p;
    this-> n = n;
    this->nnz = 0;
    this->nzVec = x;
    this->indexCol = indexCol;    
    this->rowPos = rowPos;
}

// assumes that it is the right class; does not check
SparseMatrix::SparseMatrix(SEXP x) {
    SEXP pslot = GET_SLOT(x, mkString("p"));
    SEXP islot = GET_SLOT(x, mkString("i"));
    SEXP dimslot = GET_SLOT(x, mkString("Dim"));
    SEXP nzVecSlot = GET_SLOT(x, mkString("x"));
    nzVec = vector<double>(REAL(nzVecSlot), REAL(nzVecSlot) + LENGTH(nzVecSlot));
    rowPos = vector<int>(INTEGER(islot), INTEGER(islot) + LENGTH(islot));
    indexCol = vector<int>(INTEGER(pslot), INTEGER(pslot) + LENGTH(pslot));
    n = INTEGER(dimslot)[0];
    p = INTEGER(dimslot)[1];
    nnz = nzVec.size();
}




// for the fusions, create the fusedX matrix
// Implementation at the moment assumes that 
SparseMatrix SparseMatrix::createFusedX(const vector<vector<int> > fusedGroups) {
    // in order to make it more efficient, allocate enough space for the whole
    // X matrix; may be smaller, but for very sparse matrices, will be a good
    // approximation
    SparseMatrix fusedX;
    fusedX.n = n;
    fusedX.p = fusedGroups.size();

    fusedX.nzVec.reserve(nzVec.size());
    fusedX.nzVec.clear();
    fusedX.indexCol.resize(fusedGroups.size() + 1, 0);
    fusedX.rowPos.reserve(rowPos.size());
    fusedX.rowPos.clear();

    vector<int> pos;
    vector<int> endPos;
    int curRow;
    double curVal;

    for(int group = 0; group < fusedGroups.size(); ++group) {
        curRow = 0;
        pos.clear(); pos.resize(fusedGroups[group].size()); // doesn't need to be initialized
        endPos.clear(); endPos.resize(fusedGroups[group].size()); // as immediately below
        // initialize the pos and endPos vector
        for(int i = 0; i < fusedGroups[group].size(); i++) {
            pos[i] = indexCol[fusedGroups[group][i]];
            endPos[i] = indexCol[fusedGroups[group][i] + 1] - 1;
        }

        // save the starting position of this column
        fusedX.indexCol[group] = fusedX.nzVec.size();

        while(curRow < n) {
            curRow = n;
            // find the smallest row available
            for(int i = 0; i < pos.size(); ++i) {
                if(pos[i] <= endPos[i] && rowPos[pos[i]] < curRow) {
                    curRow = rowPos[pos[i]];
                }
            }

            if(curRow == n) {
                continue;
            }

            curVal = 0;
            fusedX.rowPos.push_back(curRow); 
            for(int i = 0; i < pos.size(); ++i) {
                if(pos[i] <= endPos[i] && rowPos[pos[i]] == curRow) {
                    curVal += nzVec[pos[i]];
                    pos[i]++;
                }                
            }
            fusedX.nzVec.push_back(curVal);
        }
    }
    fusedX.indexCol[fusedGroups.size()] = fusedX.nzVec.size();
    fusedX.nnz = fusedX.nzVec.size();

    return fusedX;
}


void SparseMatrix::addColumn(const vector<double>& newnzVec, const vector<int>& newrowPos) {
    nzVec.insert(nzVec.end(), newnzVec.begin(), newnzVec.end());
    rowPos.insert(rowPos.end(), newrowPos.begin(), newrowPos.end());
    nnz += newnzVec.size();
    indexCol.push_back(nnz);
    p += 1;
}


void SparseMatrix::addColumn(const vector<double>& x) {
    for(int j = 0; j < n; ++j) {
        if(x[j] != 0) {
            nnz++;
            nzVec.push_back(x[j]);
            rowPos.push_back(j);
        }
    }
    indexCol.push_back(nnz);
    p += 1;
}


double SparseMatrix::innerProd(int i) const {
    if(i < 0 || i >= p) {
      REprintf("not a valid column %d out of %d\n", i, p);
      error("\n");
    }

    double res = 0;
    for(int pos = indexCol[i]; pos < indexCol[i+1]; ++pos) {
        res += nzVec[pos] * nzVec[pos];
    }
    return res;
}

double SparseMatrix::innerProd(int i, int j) const {
    if(i < 0 || i >= p) {
      REprintf("not a valid column %d out of %d\n", i, p);
      error("\n");
    }

   if(j < 0 || j >= p) {
      REprintf("not a valid column %d out of %d\n", j, p);
      error("\n");
    }

    double res = 0;
    // first check that column j has any entries
    if(indexCol[j] == indexCol[j + 1]) { // no entries
        return 0;
    }

    int posJ = indexCol[j];
    for(int posI = indexCol[i]; posI < indexCol[i + 1]; ++posI) {
        while(rowPos[posJ] < rowPos[posI]) { // increase rows until not less than in I
            posJ++;
        }
        if(rowPos[posJ] == rowPos[posI]) {
            res += nzVec[posJ] * nzVec[posI];
        }
    }

    return res;
}

double SparseMatrix::innerProd(const vector<double>& y, int i) const {
    if(y.size() != n) { // not the right size
      error("y does not have the right size \n");
    }
    if(i < 0 || i >= p) {
      REprintf("not a valid column %d out of %d\n", i, p);
      error("\n");
    }

    double res = 0;
    for(int posI = indexCol[i]; posI < indexCol[i + 1]; ++posI) {
        res += y[rowPos[posI]] * nzVec[posI];
    }

    return res;
}

double SparseMatrix::innerProd(int i, const vector<double>& w) const {
    if(i < 0 || i >= p) {
      REprintf("not a valid column %d out of %d\n", i, p);
      error("\n");
    }

    if(w.size() != n) {
      error("w not the right size\n");
    }

   double res = 0;
    for(int pos = indexCol[i]; pos < indexCol[i+1]; ++pos) {
        res += nzVec[pos] * nzVec[pos] * w[rowPos[pos]];
    }
    return res;
}

double SparseMatrix::innerProd(int i, int j, const vector<double>& w) const {
    if(i < 0 || i >= p) {
      REprintf("not a valid column %d out of %d\n", i, p);
      error("\n");
    }

   if(j < 0 || j >= p) {
      REprintf("not a valid column %d out of %d\n", j, p);
      error("\n");
    }

   if(w.size() != n) {
     error("w not the right size\n");
   }

    double res = 0;
    // first check that column j has any entries
    if(indexCol[j] == indexCol[j + 1]) { // no entries
        return 0;
    }

    int posJ = indexCol[j];
    for(int posI = indexCol[i]; posI < indexCol[i + 1]; ++posI) {
        while(rowPos[posJ] < rowPos[posI]) { // increase rows until not less than in I
            posJ++;
        }
        if(rowPos[posJ] == rowPos[posI]) {
            res += nzVec[posJ] * nzVec[posI] * w[rowPos[posI]];
        }
    }

    return res;
}

double SparseMatrix::innerProd(const vector<double>& y, int i, const vector<double>& w) const {
    if(y.size() != n) { // not the right size
      error("y does not have the right size\n");
    }
    if(i < 0 || i >= p) {
      REprintf("not a valid column %d out of %d\n", i, p);
      error("\n");
    }

    if(w.size() != n) {
      error("w not the right size\n");
    }

    double res = 0;
    for(int posI = indexCol[i]; posI < indexCol[i + 1]; ++posI) {
        res += y[rowPos[posI]] * nzVec[posI] * w[rowPos[posI]];
    }

    return res;
}

void SparseMatrix::multByDiag(const vector<double>& w) {
    // check that w has the right size
    if(w.size() != n) { // not the right size
      error("w does not have the right size\n");
    }

    // no multiply the matrix correctly
    for(int i = 0; i < nnz; ++i) {
        nzVec[i] *= w[rowPos[i]];
    }
}

// need to check later if these checks at the beginning result in a
// speed penalty
void SparseMatrix::addMultOfColumn(vector<double>& y, int i, double mult) const {
    // check that y has the right size
    if(y.size() != n) { // does not have the right size
      error("y does not have the right size\n");
    }
    if(i < 0 || i >= p) {
      REprintf("not a valid column %d out of %d\n", i, p);
      error("\n");
    }

    if(mult == 0) {
        return;
    }

    for(int pos = indexCol[i]; pos < indexCol[i + 1]; ++pos) {
        y[rowPos[pos]] += nzVec[pos] * mult;
    }
}

// also does checking if the column is ok - not good for performance
// but this function is mainly intended for error checking using the
// test suite
double SparseMatrix::get(int i, int j) const {
    if(j < 0 || j >= p) {
      REprintf("not a valid column %d out of %d\n", j, p);
      error("\n");
    }

    for(int pos = indexCol[j]; pos < indexCol[j + 1]; ++pos) {
        
        if(rowPos[pos] == i) {
            return nzVec[pos];
        }
        if(rowPos[pos] > i) { // could not find it
            return 0;
        }
    }
    return 0; // could not find it
}


vector<double> SparseMatrix::getColumn(int i) const {
    vector<double> res(n);

    if(i < 0 || i >= p) {
      REprintf("not a valid column %d out of %d\n", i, p);
      error("\n");
    }

    for(int pos = indexCol[i]; pos < indexCol[i + 1]; ++pos) {
        res[rowPos[pos]] = nzVec[pos];
    }
    return res;
}

SEXP SparseMatrix::todgCMatrix() const {
    SEXP islot, pslot, xslot, dimslot, dgCMatrix, class_def;
    PROTECT(islot = allocVector(INTSXP, nnz));
    PROTECT(pslot = allocVector(INTSXP, p + 1));
    PROTECT(xslot = allocVector(REALSXP, nnz));
    PROTECT(dimslot = allocVector(INTSXP, 2));

    // copy everything over from the SparseMatrix
    INTEGER(dimslot)[0] = n;
    INTEGER(dimslot)[1] = p;
    for(int i = 0; i < indexCol.size(); ++i) {
        INTEGER(pslot)[i] = indexCol[i];
    }
    for(int i = 0; i < nzVec.size(); ++i) {
        REAL(xslot)[i] = nzVec[i];
        INTEGER(islot)[i] = rowPos[i];
    }

    class_def = MAKE_CLASS("dgCMatrix");
    PROTECT(dgCMatrix = NEW_OBJECT(class_def));
    SET_SLOT(dgCMatrix, mkString("i"), islot);
    SET_SLOT(dgCMatrix, mkString("p"), pslot);
    SET_SLOT(dgCMatrix, mkString("x"), xslot);
    SET_SLOT(dgCMatrix, mkString("Dim"), dimslot);

    UNPROTECT(5);
    return(dgCMatrix);
}

void SparseMatrix::print(ostream& outStream) const {
    outStream << "---------------------------------------------" << endl;
    outStream << "nnz: " << nnz << " n: " << n << " p: " << p << endl;
    outStream << "nzVec:" << endl;
    printVector(nzVec, outStream);
    outStream << "indexCol:" << endl;
    printVector(indexCol, outStream);
    outStream << "rowPos:" << endl;
    printVector(rowPos, outStream);
    outStream << "---------------------------------------------" << endl;
}


void SparseMatrix::printColumn(int i, ostream& outStream) const {
    outStream << "Num non zero: " << indexCol[i+1] - indexCol[i] << endl;
    for(int pos = indexCol[i]; pos < indexCol[i + 1]; ++pos) {
        outStream << " Row: " << rowPos[pos] << " Val: " << nzVec[pos];
    }
    outStream << endl;
}
