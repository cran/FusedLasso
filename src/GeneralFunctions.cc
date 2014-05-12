#include "GeneralFunctions.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <set>
#include <limits>
#include <math.h>

using namespace std;

double RelDif(double a, double b)
{
	double c = Abs(a);
	double d = Abs(b);

	d = Max(c, d);

	return d == 0.0 ? 0.0 : Abs(a - b) / d;
}

double RelDifNoAbs(double a, double b)
{
	double c = Abs(a);
	double d = Abs(b);

	d = Max(c, d);

	return d == 0.0 ? 0.0 : (a - b) / d;
}


double maxDiffDoubleVec(const vector<double>& x, const vector<double>& y) {
    double maxDiff = 0;
    if(x.size() != y.size()) {
        return infinite;
    }
    for(int i = 0; i < x.size(); ++i) {
        if(fabs(x[i] - y[i]) > maxDiff) {
            maxDiff = fabs(x[i] - y[i]);
        }
    }
    return maxDiff;
}


void readDoubleMatrix(vector<double>& X, int& n, int& p, string filename) {
    n = 0;
    p = 0;
    double foo;
    ifstream infile;

    X.clear();

    infile.open(filename.c_str());
    //if(!infile) {
    //  REprintf("Couldn't open file %s !\n", filename.c_str());
    //  error("Could not open file!\n");
    //} 
    while(!infile.eof()) {
        if(infile.peek() == '\n') {
            n++;
        }    
        infile >> foo;
        X.push_back(foo);
        p++;
    }
    // check if there is one too many, then delete it
    while(p % n != 0) {
        X.pop_back();
        p--;
    }
    p /= n;
    X = transpose(X, p, n);
}

// does not do any error checking
vector<double> transpose(vector<double>& X, int n, int p) {
    vector<double> res(X.size());

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < p; ++j) {
            res[i * p + j] = X[j * n + i];            
        }
    }

    return res;
}

bool checkIfSubset(const set<int>& x, const vector<int>& fusions) {
    set<int>::const_iterator setIt;
    int group = fusions[*(x.begin())];
    for(setIt = x.begin(); setIt != x.end(); ++setIt) {
        if(fusions[*setIt] != group) {
            return false;
        }
    }
    return true; 
}

int numNonZero(const vector<double>& x) {
    int nonZero = 0;
    for(int i = 0; i < x.size(); ++i) {
        if(x[i]!=0) {
            nonZero++;
        }
    }
    return nonZero;
}


void printVector(const vector<double> x, ostream& outStream) {
    for(int i = 0; i < x.size(); ++i) {
        outStream << i <<":" <<x[i] << ", ";
    }
    outStream << endl;
}

void printVector(const vector<int> x, ostream& outStream) {
    for(int i = 0; i < x.size(); ++i) {
        outStream << i <<":" << x[i] << ", ";
    }
    outStream << endl;
}

void printList(const list<double> x, ostream& outStream) {
    list<double>::const_iterator listIt; 
    for(listIt = x.begin(); listIt != x.end(); ++listIt) {
        outStream << *listIt << ", ";
    }
    outStream << endl;
}

void printList(const list<int> x, ostream& outStream) {
    list<int>::const_iterator listIt; 
    for(listIt = x.begin(); listIt != x.end(); ++listIt) {
        outStream << *listIt << ", ";
    }
    outStream << endl;
}


void printMatrix(vector<double>& X, int n, int p, ostream& outStream) {
    outStream << "n: " << n << "   p: " << p << endl;
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < p; ++j) {
            outStream << X[j*n + i] << ", ";
        }
        outStream << endl;
    }
}


