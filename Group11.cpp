#include "Group11.h"

/*
 Define the function signum which returns -1 for any value less than 0; +1 for any value greater than 0
 and 0 for an input 0.
 */

template <typename T> int sign(T val)
{
    return (T(0) < val) - (val < T(0));
}


/*
 
 BASE CLASS : SCHEMES
 
 ===================================================================================================================================
 This Class contains all the common attributes and methods for all the four methods which are derived classes of this base class.
 It takes in the Number of Space Grid Points, N and the required Time Level, T. This constructor is modified in the derived classes
 to solve for the four schemes. The common methods selected are: gridPoints; initialConditions; analyticalSolutions;
 virtual numericalSolutions and findError. The numericalSolutions was made virtual as the implementation in the derived classes
 would depend on the scheme formulation. In the base class, the numericalSolutions method only prints an information on the screen
 indicating that the base class has no implementation and that an object of a derived class should be instantiated instead.
 ===================================================================================================================================
 
 */

Schemes::Schemes(int N, int T, string fileName)
{
    //The arguments given are equated to the protected data members
    this -> N = N;
    this -> T = T;
    this -> fileName = fileName;
    
    
    //Initialize the Given Conditions
    u = 1.75;
    L = 100;
    delx = L/N;
    delt = 0.1;
    f = N + 1;
    setOneBC[0] = 0; setOneBC [1] = 1;
    setTwoBC[0] = 0; setTwoBC [1] = 0;
    
    x = vector<double>(f);
    
    //Call the required methods in the constructor
    gridPoints();
    initialConditions();
    analyticalSolution();
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 gridPoints: This method is to specify the values of the array of the space grid points which depends on delx. delx is simply gotten
 from the number of space grid points, N in the tube length, L. These values are included in a vector<double> array, x.
 -----------------------------------------------------------------------------------------------------------------------------------
 */

void Schemes::gridPoints()
{
    x[0] = -50; //Initialize the Left Boundary
    for (int i = 1; i < f; i++)
    {
        x[i] = x[i - 1] + delx;
    }
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 initialConditions: This method is to get the values of the initial conditions of both Set One and Set Two using the Initial
 Condition provided which covers for the values of x at -L/2 and L/2 during any time step as well as the solutions of the functions
 at the time step, t = 0
 -----------------------------------------------------------------------------------------------------------------------------------
 */

void Schemes::initialConditions()
{
    for (int i = 0; i < f; i++)
    {
        setOneIC.push_back(0.5 * (sign(x[i]) + 1));
        setTwoIC.push_back((0.5 * exp(-(pow(x[i], 2)))));
    }
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 analyticalSolution: This method is used to get the values of the analytical solution at a specified time step, t using the specified
 formula for both Set One and Set Two. This result would be compared with the Numerical Solution at that same time step to solve for
 the errors and see the behaviour of the the different schemes
 -----------------------------------------------------------------------------------------------------------------------------------
 */

void Schemes::analyticalSolution()
{
    for (int i = 0; i < f; i++)
    {
        setOneAS.push_back(0.5 * (sign(x[i] - (1.75 * T)) + 1));
        setTwoAS.push_back(0.5 * exp(-(pow((x[i] - (1.75 * T)), 2))));
    }
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 numericalSolution: The only virtual method in this class is used to solve for the numerical solutions of the different schemes. It
 was chosen to be virtual because different formulas are utilized in the different schemes. A simple message is printed on the screen
 in the base class in order to inform the user to choose one of the derived classes for the object.\
 -----------------------------------------------------------------------------------------------------------------------------------
 */

void Schemes::numericalSolution()
{
    
    cout << "This should solve for the Numerical Solution in each Scheme " << endl;
    cout << "Kindly choose an appropriate scheme for your object in order to get a solution" << endl;
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 findError: This method is used to find the absolute value of the difference between the analytical solution and numerical solution
 in order to qualitatively and quantitavely examine the schemes
 -----------------------------------------------------------------------------------------------------------------------------------
 */

void Schemes::findError()
{
    
    for (int i = 0; i < f; i++)
    {
        errorOne.push_back(abs(setOneAS[i] - setOneNS[i]));
        errorTwo.push_back(abs(setTwoAS[i] - setTwoNS[i]));
    }
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 L0Norm: This method is used to find the maximum value of the difference between the analytical solution and numerical solution
 in order to qualitatively and quantitavely examine the schemes
 -----------------------------------------------------------------------------------------------------------------------------------
 */
double Schemes::L0Norm(vector <double>&s)
{
    L0 = 0;
    for ( int i=0; i<f; i++ )
    if (L0 < s[i])
    {
        L0 = s[i];
    }
    return L0;
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 L1Norm: This method is used to find the mean of the difference between the analytical solution and numerical solution
 in order to qualitatively and quantitavely examine the schemes
 -----------------------------------------------------------------------------------------------------------------------------------
 */

double Schemes::L1Norm(vector <double>&s)
{
    L1 = 0;
    for ( int i=0; i<f; i++ )
    {
        L0 += s[i];
    }
    return L0/s.size();
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 findError: This method is used to find the mean square of the difference between the analytical solution and numerical solution
 in order to qualitatively and quantitavely examine the schemes
 -----------------------------------------------------------------------------------------------------------------------------------
 */
double Schemes::L2Norm(vector <double>&s)
{
    L2 = 0;
    for ( int i=0; i<f; i++ )
    {
        L2 += pow(s[i], 2);
    }
    return L2/s.size();
}

void Schemes::displayResults()
{
    cout << "\t\t" << left << setw(21) << "SET ONE SOLUTION" << setw(21) << "SET TWO SOLUTION" << endl;
    cout << "--------------------------------------------------" << endl;
    cout << "L0 =\t" << left << setw(21) << L0Norm(errorOne) << setw(21) << L0Norm(errorTwo) << endl;
    cout << "L1 =\t" << left << setw(21) << L1Norm(errorOne) << setw(21) << L1Norm(errorTwo) << endl;
    cout << "L2 =\t" << left << setw(21) << L2Norm(errorOne) << setw(21) << L2Norm(errorTwo) << endl;
    cout << endl;
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 createFile: This method is used to store the results of the Numerical Solution, Analytical Solution and Error in a csv file. The
 file is reset everytime the code is run
 -----------------------------------------------------------------------------------------------------------------------------------
 */

void Schemes::createFile()
{
    fstream fout;
    fout.open(fileName.c_str(), ios::out | ios::trunc);
    
    fout << "SET ONE SOLUTION" << ",,,,," << "SET TWO SOLUTION" << "\n" << endl;
    fout << "x[i]" << ", " << "Analytical Solution" << ", " << "Numerical Solution" << "," << "Error" << ",,";
    fout << "x[i]" << ", " << "Analytical Solution" << ", " << "Numerical Solution" << "," << "Error" << "\n" << endl;
    for (int i = 0; i < f; i++)
    {
        fout << x[i] << ", " << setOneAS[i] << ", " << setOneNS[i] << ", " << errorOne[i] << ",, "
        << x[i] << ", " << setTwoAS[i] << ", " << setTwoNS[i] << ", " << errorTwo[i] << endl;
    }
    
    fout.close();
}

/*
 DERIVED CLASS: EXPLICIT METHOD
 
 ===================================================================================================================================
 This derived class contains two methods with two attributes. The attributes are the coefficients of the known elements. The
 constructor accepts the input variables from the base class: N and T and uses these variables to solve for the Explicit Forward-Time
 Backward-Space method in the virtual numericalSolution method.
 ===================================================================================================================================
 */

Explicit_FTBS::Explicit_FTBS(int N, int T, string fileName) : Schemes (N, T, fileName)
{
    //Coefficients of the knowns
    a = ((u * delt)/(delx));
    b = 1 - ((u * delt)/(delx));
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 numericalSolution: This virtual method in this class is used to solve for the numerical solutions of the Explicit Forward-Time
 Backward-Space scheme. The iterOne and iterTwo vector<double>s are used to duplicate the solution of the present time step so that the
 next time step can be solved for.
 -----------------------------------------------------------------------------------------------------------------------------------
 */

void Explicit_FTBS::numericalSolution()
{
    setOneNS = setOneIC;
    setTwoNS = setTwoIC;
    
    for (int j = 0; j < T/delt; j++)
    {
        vector<double> iterOne = setOneNS, iterTwo = setTwoNS;
        setOneNS.clear(); setTwoNS.clear();
        
        for (int i = 0; i < f; i++)
        {
            if (i == 0) //Check for Boundary Conditions
            {
                setOneNS.push_back(setOneBC[0]);
                setTwoNS.push_back(setTwoBC[0]);
            }
            else if (i == f - 1) //Check for Boundary Conditions
            {
                setOneNS.push_back(setOneBC[1]);
                setTwoNS.push_back(setTwoBC[1]);
            }
            else
            {
                setOneNS.push_back(a * iterOne[i-1] + b * iterOne[i]);
                setTwoNS.push_back(a * iterTwo[i-1] + b * iterTwo[i]);
            }
        }
    }
}

/*
 DERIVED CLASS: IMPLICIT METHOD
 
 ===================================================================================================================================
 This derived class contains two methods with two attributes. The attributes are the coefficients of the known elements. The
 constructor accepts the input variables from the base class: N and T and uses these variables to solve for the Implicit Forward-Time
 Backward-Space method in the virtual numericalSolution method.
 ===================================================================================================================================
 */
Implicit_FTBS::Implicit_FTBS(int N, int T, string fileName) : Schemes(N, T, fileName)
{
    a = (u * delt) / (delx + (u * delt));
    b =  delx / (delx + (u * delt));
    
    alpha = (u * delt)/delx;
    beta = 1 + ((u * delt)/delx);
    
    Imp = Matrix (T/delt, N);
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 numericalSolution: This virtual method is used to solve for the numerical solutions of the Explicit Forward-Time Backward-Space
 scheme. The iterOne and iterTwo vector<double>s are used to duplicate the solution of the present time step so that the next time step
 can be solved for.
 -----------------------------------------------------------------------------------------------------------------------------------
 */

void Implicit_FTBS::numericalSolution()
{
    setOneNS = setOneIC;
    setTwoNS = setTwoIC;
    
    for (int j = 0; j < T/delt; j++)
    {
        vector<double> iterOne = setOneNS, iterTwo = setTwoNS;
        setOneNS.clear(); setTwoNS.clear();
        
        for (int i = 0; i < f; i++)
        {
            if (i == 0) //Check for Boundary Conditions
            {
                setOneNS.push_back(setOneBC[0]);
                setTwoNS.push_back(setTwoBC[0]);
            }
            else if (i == f - 1) //Check for Boundary Conditions
            {
                setOneNS.push_back(setOneBC[1]);
                setTwoNS.push_back(setTwoBC[1]);
            }
            else
            {
                setOneNS.push_back(a * setOneNS[i - 1] + b * iterOne[i]);
                setTwoNS.push_back(a * setTwoNS[i - 1] + b * iterTwo[i]);
            }
        }
    }
}

void Implicit_FTBS::matrixSolution() 
{
    for (int i = 0; i < T/delt; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (j == i) Imp[i][j] = beta;
            else if (j == (i - 1)) Imp[i][j] = alpha;
            else Imp[i][j] = 0;
        }
    }
}

void Implicit_FTBS::createMatrix() const
{
    fstream mout;
    mout.open("matrices.csv", ios::out | ios::trunc);
    
    for (int i = 0; i < T/delt; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mout << Imp[i][j] << ',';
        }
        mout << endl;
    }
    mout.close();
}
/*
 
 DERIVED CLASS: LAX WENDROFF METHOD
 
 ===================================================================================================================================
 This derived class contains a constructor, a method and three attributes. The attributes are the coefficients of known elements. The
 constructor accepts the input variables from the base class: N and T and uses these variables to solve for the Lax-Wendroff Scheme
 in the virtual numericalSolution method.
 ===================================================================================================================================
 
 */
Lax_Wendroff::Lax_Wendroff(int N, int T, string fileName) : Schemes(N, T, fileName)
{
    a = ((pow(delt * u, 2)) + (u * delx * delt))/(2 * pow(delx, 2));
    b = ((pow(delx, 2)) - (pow(delt * u, 2)))/(pow(delx, 2));
    c = ((pow(delt * u, 2)) - (u * delx * delt)) / (2 * pow(delx, 2));
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 numericalSolution: This virtual method is used to solve for the numerical solutions of the Explicit Lax-Wendroff scheme. The iterOne
 and iterTwo vector<double>s are used to duplicate the solution of the present time step so that the next time step can be solved for.
 -----------------------------------------------------------------------------------------------------------------------------------
 */
void Lax_Wendroff::numericalSolution()
{
    setOneNS = setOneIC;
    setTwoNS = setTwoIC;
    
    for (int j = 0; j < T/delt; j++)
    {
        vector<double>iterOne = setOneNS;
        vector<double>iterTwo = setTwoNS;
        setOneNS.clear();
        setTwoNS.clear();
        
        for (int i = 0; i < f; i++)
        {
            if (i == 0) //Check for Boundary Conditions
            {
                setOneNS.push_back(setOneBC[0]);
                setTwoNS.push_back(setTwoBC[0]);
                //setOneNS[f - 1] = setOneBC[1];
            }
            else if (i == f - 1) //Check for Boundary Conditions
            {
                setOneNS.push_back(setOneBC[1]);
                setTwoNS.push_back(setTwoBC[1]);
            }
            else
            {
                setOneNS.push_back((a * iterOne[i - 1]) + (b * iterOne[i]) + (c * iterOne[i + 1]));
                setTwoNS.push_back((a * iterTwo[i - 1]) + (b * iterTwo[i]) + (c * iterTwo[i + 1]));
            }
        }
    }
}

/*
 
 DERIVED CLASS: RICHTMYER METHOD
 
 ===================================================================================================================================
 This derived class contains a constructor, a method and five attributes. The attributes are the coefficients of known elements and
 predictor elements. The constructor accepts the input variables from the base class: N and T and uses these variables to solve for
 the Lax-Wendroff Scheme in the virtual numericalSolution method.
 ===================================================================================================================================
 
 */
Richtmyer::Richtmyer(int N, int T, string fileName) : Schemes(N, T, fileName)
{
    c = (u * delt) / delx;
    d = (2 - c) / 4;
    e = (2 + c) / 4;
    g = c / 2;
    h = -g;
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 predictorSolution: This virtual method is used to solve for the numerical solutions of the Explicit Lax-Wendroff scheme. The
 predictorOne and predictorTwo vector<double>s are used to find the half time steps above of the present time step so that the next time
 step can be solved for.
 -----------------------------------------------------------------------------------------------------------------------------------
 */
vector< vector<double> > Richtmyer::predictorSolution( vector<double> & iterOne, vector<double> &iterTwo)
{
    
    
    vector< vector<double> > results;
    vector <double> predOne;
    vector <double> predTwo;

    for(int i = 0; i < f ; ++i)
    {
        
        if (i == 0 ) //Check for boundary conditions
        {
            
            predOne.push_back(setOneBC[0]);
            predTwo.push_back(setTwoBC[0]);
            
        }
        else if(i == (f-1)) //Check for boundary conditions
        {
            predOne.push_back(setOneBC[1]);
            predTwo.push_back(setTwoBC[1]);
            
        }
        else
        {
            
            predOne.push_back((d*iterOne[i+1]) + (e*iterOne[i-1]));
            predTwo.push_back((d*iterTwo[i+1]) + (e*iterTwo[i-1]));
        }
    }
    results.push_back(predOne);
    results.push_back(predTwo);
    
    return results;
}

/*
 -----------------------------------------------------------------------------------------------------------------------------------
 numericalSolution: This virtual method is used to solve for the numerical solutions of the Explicit Forward-Time Backward-Space
 scheme. This scheme is carried out in two steps: the predictor and the corrector. The predictor vector contains the values at the
 time step: n+1/2. The iterOne and iterTwo vectors are used to duplicate the solution of the present time step so that the next time
 step can be solved for.
 -----------------------------------------------------------------------------------------------------------------------------------
 */
void Richtmyer::numericalSolution()
{
    

    vector<double> predictorOne;
    vector<double> predictorTwo;
    
    setOneNS = setOneIC;
    setTwoNS = setTwoIC;
    
    for(int j = 0; j < (T/delt); j++)
    {
        
        vector<double> iterOne  = setOneNS;
        vector<double> iterTwo  = setTwoNS;
        setOneNS.clear();
        setTwoNS.clear();
        
        //Half time step
        predictorOne = predictorSolution(iterOne, iterTwo)[0];
        predictorTwo = predictorSolution(iterOne, iterTwo)[1];
        
        
        //Next Time Step
        for(int i = 0; i < f; i++)
        {
            if (i == 0 )
            {
                
                setOneNS.push_back(setOneBC[0]);
                setTwoNS.push_back(setTwoBC[0]);
            }
            else if(i == (f-1))
            {
                setOneNS.push_back(setOneBC[1]);
                setTwoNS.push_back(setTwoBC[1]);
                
                
            } else
            {
                setOneNS.push_back((g*predictorOne[i-1]) + (h*predictorOne[i+1]) + (iterOne[i]));
                setTwoNS.push_back((g*predictorTwo[i-1]) + (h*predictorTwo[i+1]) + (iterTwo[i]));
            }
        }
        
        predictorOne.clear();
        predictorTwo.clear();
        
    }
}

