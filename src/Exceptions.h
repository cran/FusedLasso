#ifndef _EXCEPTIONS_
#define _EXCEPTIONS_

#include <exception>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>

using namespace std;

class DidNotConvergeException : public exception {
   
    int iteration;
    double accuracy;
    string message;
public:

    DidNotConvergeException(int iteration, double accuracy) throw () { 
        this->iteration = iteration;
        this->accuracy = accuracy; 
        ostringstream s1;
        s1 << "Did not converge. Iteration: " << iteration << " Accuracy: " << accuracy << endl;
        message = s1.str();
    }
   
    ~DidNotConvergeException() throw () {
    }

    virtual const char* what() const throw ()
    {
        return(message.c_str());
    
    }
};


class DoubleGroupedException : public exception {
    int node;
    string message;

public:
    DoubleGroupedException(int node) throw () {
        this->node = node;
        ostringstream s1;
        s1 << "A node was grouped twice. Node: " << node << endl;
        message = s1.str();
    }

    ~DoubleGroupedException() throw() {
    }

    virtual const char* what() const throw() {
        return(message.c_str());
    }

};

class NodeDoesNotExistException : public exception {
    int node;
    string message;

public:
    NodeDoesNotExistException(int node) throw () {
        this->node = node;
        ostringstream s1;
        s1 << "The node did not exist in MaxFlowGraph. Node: " << node << endl;
        message = s1.str();
    }

    ~NodeDoesNotExistException() throw() {
    }

    virtual const char* what() const throw() {
        return(message.c_str());
    }

};

class OutOfBoundsException : public exception {
    int node;
    string message;

public:
    OutOfBoundsException(int node) throw () {
        this->node = node;
        ostringstream s1;
        s1 << "The node has a group that is too large. Node: " << node << endl;
        message = s1.str();
    }

    ~OutOfBoundsException() throw() {
    }

    virtual const char* what() const throw() {
        return(message.c_str());
    }

};


#endif // _EXCEPTIONS_
