#ifndef GROUP11_H
#define GROUP11_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <numeric>
#include <string>
#include "matrix.h"

using namespace std;

class Schemes
{
protected:
    int N ; //Number of Space Grid Points
    int T; //Value of the Maximum Time Step
    double delt; //Value of the Time Step Increase
    double u; //Speed of Sound
    double L; //Length of Long Tube
    double delx; //Size of the Grid Points
    double f; //Number of Space Grids
    string fileName;
    double setOneBC[2], setTwoBC[2]; //Left and Right Boundary Conditions
    
    vector<double> x; //Array of Space Grid Points
    vector<double> setOneIC, setTwoIC; //Set One Initial Condition results
    vector<double> setOneAS, setTwoAS; //Analytical Results
    vector<double> setOneNS, setTwoNS; //Numerical Results
    vector<double> errorOne, errorTwo; //Arrays for Errors
    double L0, L1, L2;
    
public:
    Schemes(int, int, string);
    void gridPoints(); //Set the Size of Grid Points (Edit)
    void initialConditions(); //Set the value of the Initial Conditions
    virtual void numericalSolution();
    void analyticalSolution(); //Set the values of the Analytical Solutions at a Time Step
    void findError(); //Find the values of the errors
    double L0Norm(vector<double>&x); //Maximum norm
    double L1Norm(vector<double>&x); //Mean Norm
    double L2Norm(vector<double>&x); //Mean Square Norm
    void displayResults(); //Display Norms on Screen
    void createFile(); //Create File to show the solutions and errors
    ~ Schemes(){};
};

/*
 EXPLICIT SCHEME
 */
class Explicit_FTBS : public Schemes
{
private:
	double a; //Coefficient of First Term
	double b; //Coefficient of Second Term
	
public:
    Explicit_FTBS(int, int, string);
    virtual void numericalSolution();
};


/*
 IMPLICIT SCHEME
 */
class Implicit_FTBS : public Schemes
{
private:
    double a, alpha;
    double b, beta;
    
    Matrix Imp;
    
public:
    Implicit_FTBS(int, int, string);
    virtual void numericalSolution();
    void matrixSolution();
    void createMatrix() const;
};

/*
 LAX WENDROFF
 */
class Lax_Wendroff : public Schemes
{
private:
    double a;
    double b;
    double c;
    
public:
    Lax_Wendroff(int, int, string);
    virtual void numericalSolution();
};

/*
 RICHTMEYER SCHEME
 */
class Richtmyer : public Schemes
{
private:
    double c;
    double d;
    double e;
    double g;
    double h; 
public:
    Richtmyer(int, int, string);
    vector< vector <double> > predictorSolution(vector<double>&, vector<double>&);
    virtual void numericalSolution();
};


#endif
