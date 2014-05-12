#ifndef _GENERALFUNCTIONS_
#define _GENERALFUNCTIONS_

#include <limits>
#include <string>
#include <vector>
#include <list>
#include <set>

using namespace std;

#define Abs(x)    ((x) < 0 ? -(x) : (x))
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))

const double tolerance =1.0e-8;
const double infinite = std::numeric_limits<double>::max();
const int infiniteInt = std::numeric_limits<int>::max();

double RelDif(double a, double b);
double RelDifNoAbs(double a, double b);

double maxDiffDoubleVec(const vector<double>&x, const vector<double>& y);

inline int signum(double x) {return((x>0)-(x<0));};

int numNonZero(const vector<double>& x);

// read in a matrix from file filename
void readDoubleMatrix(vector<double>& X, int& n, int& p, string filename);

// helper function for readDoubleMatrix -> transposes the input
// needed as in R matrices are in column major order
vector<double> transpose(vector<double>& X, int n, int p);

// checks if all positions in x belongs to the same group in fusions
bool checkIfSubset(const set<int>& x, const vector<int>& fusions);

void printVector(const vector<double> x, ostream& outStream); 

void printVector(const vector<int> x, ostream& outStream); 

void printList(const list<double> x, ostream& outStream); 

void printList(const list<int> x, ostream& outStream); 

void printMatrix(vector<double>& X, int n, int p, ostream& outStream);

#endif
