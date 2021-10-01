#include <iostream>
//#include <stdc++.h>
#include <time.h>
#include "Group11.h"

using namespace std;

int main()
{
    clock_t start, end;
    
    
    cout << "------------------------COMPUTATIONAL METHODS ASSIGNMENT------------------------" << endl;
    
    cout << endl;
    
    cout << "Group 11 Members:\tADAM MAIMONA" << endl;
    cout << "\t\t\t\t\tAKINYELE MARVELLOUS" << endl;
    cout << "\t\t\t\t\tBOUJUT CLÉMENT" << endl;

    cout << endl;
    
    cout << "Objective: Examine the Application of Numerical Schemes for the Solution of Partial Differential Equations" << endl;
    
    cout << endl;
    
    
    //PART A:
    cout << "Part A : " << endl;
    cout << "Write a C++ program which solves the above problem on a uniform grid with the prescribed initial and boundary conditions using the following methods : " << endl;
    cout << "In all cases the solution is to be printed and plotted for all x locations at time levels t = 5 and t = 10." << endl;
    
    cout << endl;
    
    cout << "Explicit Forward Time Backward Space Scheme :" << endl;
    cout << "At N = 100, T = 5 and Result File: ExpN100T5.csv" << endl;
    Explicit_FTBS a = Explicit_FTBS (100, 5, "ExpN100T5.csv");
    a.numericalSolution();
    a.findError();
    a.displayResults();
    a.createFile();
    cout << endl;
    
    cout << "Explicit Forward Time Backward Space Scheme :" << endl;
    cout << "At N = 100, T = 10 and Result File: ExpN100T10.csv" << endl;
    
    start = clock();
    
    Explicit_FTBS b = Explicit_FTBS (100, 10, "ExpN100T10.csv");
    b.numericalSolution();
    b.findError();
    b.displayResults();
    b.createFile();
    
    end = clock();
    double bTime = double(end - start) / double(CLOCKS_PER_SEC); //Computation Time for b
    
    cout << endl;
    
    cout << "Implicit Forward Time Backward Space Scheme :" << endl;
    cout << "At N = 100, T = 5 and Result File: ImpN100T5.csv" << endl;
    Implicit_FTBS c = Implicit_FTBS (100, 5, "ImpN100T5.csv");
    c.numericalSolution();
    c.findError();
    c.displayResults();
    c.createFile();
    
    cout << endl;
    
    cout << "Implicit Forward Time Backward Space Scheme :" << endl;
    cout << "At N = 100, T = 10 and Result File: ImpN100T10.csv" << endl;
    
    start = clock();
    
    Implicit_FTBS d = Implicit_FTBS (100, 10, "ImpN100T10.csv");
    d.numericalSolution();
    d.findError();
    d.displayResults();
    d.createFile();

    end = clock();
    double dTime = double(end - start) / double(CLOCKS_PER_SEC); //Computation Time for d
    
    cout << endl;
    
    cout << "Lax Wendroff Scheme :" << endl;
    cout << "At N = 100, T = 5 and Result File: LaxWN100T5.csv" << endl;
    Lax_Wendroff e = Lax_Wendroff (100, 5, "LaxWN100T5.csv");
    e.numericalSolution();
    e.findError();
    e.displayResults();
    e.createFile();
    
    cout << endl;
    
    cout << "Lax Wendroff Scheme :" << endl;
    cout << "At N = 100, T = 10 and Result File: LaxWN100T10.csv" << endl;
    Lax_Wendroff f = Lax_Wendroff (100, 10, "LaxWN100T10.csv");
    f.numericalSolution();
    f.findError();
    f.displayResults();
    f.createFile();
    
    cout << endl;
    
    cout << "Richtmyer Scheme :" << endl;
    cout << "At N = 100, T = 5 and Result File: RichN100T5.csv" << endl;
    Richtmyer g = Richtmyer (100, 5, "RichN100T5.csv");
    g.numericalSolution();
    g.findError();
    g.displayResults();
    g.createFile();
    
    cout << endl;
    
    cout << "Richtmyer Scheme :" << endl;
    cout << "At N = 100, T = 10 and Result File: RichN100T5.csv" << endl;
    Richtmyer h = Richtmyer (100, 10, "RichN100T10.csv");
    h.numericalSolution();
    h.findError();
    h.displayResults();
    h.createFile();
    
    cout << endl;
    cout << endl;
    
    //PART B:
    cout << "Part B : " << endl;
    cout << "Investigate the effect of various number of space grid points N on the accuracy of the solution, as well as the required computation time, of one explicit and one implicit method of your choice, for t = 10 and appropriate deltat of your choice:" << endl;
    cout << "In all cases the solution is to be printed and plotted for all x locations at time levels t = 5 and t = 10." << endl;
    
    //Since Explicit and Implicit Schemes have been done in previous part, the result and computation time would only be displayed in this part.
    cout << "Explicit Forward Time Backward Space Scheme :" << endl;
    cout << "At N = 100, T = 10 and Result File: ExpN100T10.csv" << endl;
    b.displayResults();
    cout << bTime << endl;
    
    cout << endl;
    
    cout << "Implicit Forward Time Backward Space Scheme :" << endl;
    cout << "At N = 100, T = 10 and Result File: ImpN100T10.csv" << endl;
    d.displayResults();
    cout << dTime << endl;
    
    cout << endl;
    
    start = clock();
    
    cout << "Explicit Forward Time Backward Space Scheme :" << endl;
    cout << "At N = 200, T = 10 and Result File: ExpN200T10.csv" << endl;
    Explicit_FTBS i = Explicit_FTBS (200, 10, "ExpN200T10.csv");
    i.numericalSolution();
    i.findError();
    i.displayResults();
    i.createFile();
    
    end = clock();
    
    double iTime = double(end - start) / double(CLOCKS_PER_SEC); //Computation Time for i
    cout << iTime << endl;
    
    cout << endl;
    
    start = clock();
    
    cout << "Explicit Forward Time Backward Space Scheme :" << endl;
    cout << "At N = 400, T = 10 and Result File: ExpN400T10.csv" << endl;
    Explicit_FTBS j = Explicit_FTBS (400, 10, "ExpN400T10.csv");
    j.numericalSolution();
    j.findError();
    j.displayResults();
    j.createFile();
    
    end = clock();
    
    double jTime = double(end - start) / double(CLOCKS_PER_SEC); //Computation Time for j
    cout << jTime << endl;
    
    cout << endl;
    
    start = clock();
    
    cout << "Implicit Forward Time Backward Space Scheme :" << endl;
    cout << "At N = 200, T = 10 and Result File: ImpN200T10.csv" << endl;
    Implicit_FTBS k = Implicit_FTBS (200, 10, "ImpN200T10.csv");
    k.numericalSolution();
    k.findError();
    k.displayResults();
    k.createFile();
    
    end = clock();
    
    double kTime = double(end - start) / double(CLOCKS_PER_SEC); //Computation Time for k
    cout << kTime << endl;
    
    cout << endl;
    
    start = clock();
    
    cout << "Implicit Forward Time Backward Space Scheme :" << endl;
    cout << "At N = 400, T = 10 and Result File: ImpN400T10.csv" << endl;
    Implicit_FTBS l = Implicit_FTBS (400, 10, "ImpN400T10.csv");
    l.numericalSolution();
    l.findError();
    l.displayResults();
    l.createFile();
    
    end = clock();
    
    double lTime = double(end - start) / double(CLOCKS_PER_SEC); //Computation Time for l
    cout << lTime << endl;
    
    return 0;
}
