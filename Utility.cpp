#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <gsl/gsl_multifit.h>

#include "FileInput.h"
#include "Global.h"
#include "Utility.h"


#include "pwl_interp_2d_scattered.hpp" 


using namespace std;


void Utility::linspace(double a, double b, size_t n, double *yArray)
{
    double stepLength = (b - a) / (n - 1); 

    if( n == 1)
    {
       yArray[0] = b;
    }
    else
    {
        yArray[0] = a;
        for (size_t i = 1; i <= n-2; ++i)
        {
            yArray[i] = yArray[i-1] + stepLength;
        } 

        yArray[n-1] = b;

    } 

} 


void Utility::trim(string &s)
{

    if(!s.empty())
    {
        size_t position = s.find_first_not_of(" \t\r");

        if(position != string::npos)
        {
            s.erase(0, position); 
            s.erase(s.find_last_not_of(" \t\r") + 1);
        }
        else
        {
            s.erase(s.begin(), s.end());
        } 

    } 

} 


bool Utility::emptyOrBlanksTabs(const string &s)
{
    if (s.empty())
    {
        return true;
    }
    else if(s.find_first_not_of(" \t\r") == string::npos)
    {
        return true;
    }
    else
    {
        return false;
    } 

} 


bool Utility::isOnlyNumbers(string str)
{
    trim(str);
    if(str.empty())
        return false;


    istringstream streamOfstring(str);
    double d;

    while (streamOfstring >> d);

    if(streamOfstring.eof())
        return true;
    else
        return false;
} 


bool Utility::isOnlyOneNumber(string str)
{
    trim(str);
    if(str.empty())
        return false;


    istringstream streamOfstring(str);
    double d;

    streamOfstring >> d;

    if(streamOfstring.eof())
        return true;
    else
        return false;
} 





void Utility::scatterInterpolation2DSpaceLinear(size_t numNode, double *xyNodeArrayPtr,
                                              double *zNodeArrayPtr, size_t numInterpo,
                                              double *xyInterpoArrayPtr, double *zInterpoArrayPtr)
{
    
    int numTriangle; 
    int nodesTriangelArray[3*2*numNode]; 
    int neighborTriangle[3*2*numNode]; 

    r8tris2(numNode, xyNodeArrayPtr, numTriangle, nodesTriangelArray, neighborTriangle);  

    for ( int j = 0; j < numTriangle; j++ )
    {
        for ( int i = 0; i < 3; i++ )
        {
            if ( 0 < neighborTriangle[i+j*3] )
            {
                neighborTriangle[i+j*3] = neighborTriangle[i+j*3] - 1;
            } 
        } 
    } 


    
    double * zArrayPtrTempora; 

    zArrayPtrTempora = pwl_interp_2d_scattered_value(numNode, xyNodeArrayPtr,
                                                     zNodeArrayPtr, numTriangle,
                                                     nodesTriangelArray, neighborTriangle,
                                                     numInterpo, xyInterpoArrayPtr);   

    for(size_t i = 0; i < numInterpo; ++i)
    {
        zInterpoArrayPtr[i] = zArrayPtrTempora[i];
    } 

    delete [] zArrayPtrTempora;

} 









void Utility::estimateLogkAnalyticFormula(const double *tMeasureArrayPtr, const double *logkMeasureArrayPtr,
                                        size_t numMeasure, double *coefficientsArrayPtr,
                                        double &maxRelativeError, double &corresponAbsoluteError,
                                        double &correspondingLogk, double &corresponTemporature)
{
    
    double tMeasureArray[numMeasure];
    for (size_t i = 0; i < numMeasure; ++i)
        tMeasureArray[i] = tMeasureArrayPtr[i] + 273.15;


    
    double chisq;
    gsl_matrix *XMatrix, *covMatrix;
    gsl_vector *yVector, *cVector;

    XMatrix = gsl_matrix_alloc(numMeasure, 6);
    covMatrix = gsl_matrix_alloc(6, 6);

    yVector = gsl_vector_alloc(numMeasure);
    cVector = gsl_vector_alloc(6);

    for (size_t i = 0; i < numMeasure; ++i)
    {
        gsl_matrix_set(XMatrix, i, 0, 1);
        gsl_matrix_set(XMatrix, i, 1, tMeasureArray[i] );
        gsl_matrix_set(XMatrix, i, 2, 1.0 / tMeasureArray[i] );
        gsl_matrix_set(XMatrix, i, 3, log10(tMeasureArray[i]) );
        gsl_matrix_set(XMatrix, i, 4, 1.0 / pow(tMeasureArray[i], 2.0) );
        gsl_matrix_set(XMatrix, i, 5, pow(tMeasureArray[i], 2.0) );

        gsl_vector_set(yVector, i, logkMeasureArrayPtr[i]);
    }

    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(numMeasure, 6);
    gsl_multifit_linear(XMatrix, yVector, cVector, covMatrix, &chisq, work);
    gsl_multifit_linear_free(work);


    
    for (size_t i = 0; i < 6; ++i)
    {
        coefficientsArrayPtr[i] = gsl_vector_get(cVector, i);
    }

    
    double estimatedLogkArray[numMeasure];
    for (size_t i = 0; i < numMeasure; ++i)
        estimatedLogkArray[i] = coefficientsArrayPtr[0] +
                                coefficientsArrayPtr[1] * tMeasureArray[i] +
                                coefficientsArrayPtr[2] / tMeasureArray[i] +
                                coefficientsArrayPtr[3] * log10(tMeasureArray[i]) +
                                coefficientsArrayPtr[4] / pow(tMeasureArray[i], 2.0) +
                                coefficientsArrayPtr[5] * pow(tMeasureArray[i], 2.0);

    double relativeErrorArray[numMeasure]; 
    for (size_t i = 0; i < numMeasure; ++i)
        relativeErrorArray[i] = (estimatedLogkArray[i] - logkMeasureArrayPtr[i]) / logkMeasureArrayPtr[i] * 100;

    size_t position = max_element(relativeErrorArray, relativeErrorArray + numMeasure) - relativeErrorArray;
    maxRelativeError = relativeErrorArray[position];
    corresponAbsoluteError = estimatedLogkArray[position] - logkMeasureArrayPtr[position];
    correspondingLogk = logkMeasureArrayPtr[position];
    corresponTemporature = tMeasureArray[position] - 273.15; 


    
    gsl_matrix_free(XMatrix);
    gsl_matrix_free(covMatrix);
    gsl_vector_free(yVector);
    gsl_vector_free(cVector);

} 


void Utility::obtainVmFromSUPCRTThenStore()
{
    ifstream file( RELATIVEPATH + SPRONSFILE, ios::in); 
    ofstream outputFile(RELATIVEPATH + "vmcriticalpoint.txt", ios::out); 
    string strLine;


    
    outputFile << "\n" << "# minerals that do not undergo phase transitions" << "\n" << endl;

    while(!file.eof())
    {
        getline(file,strLine);

        if(strLine.find("minerals that do not undergo phase transitions") != string::npos)
            break;

    }

    getline(file,strLine);
    getline(file,strLine);

    size_t positionStar = strLine.find("*******************************************************");

    while(positionStar == string::npos)
    {
        string name; 

        name = strLine.substr(0, 21);
        Utility::trim(name);

        outputFile << name << endl;

        for (size_t i = 0; i < 3; ++i)
            getline(file,strLine);

        double vm; 
        if(!Utility::isOnlyOneNumber(strLine.substr(42, 10)))
        {
            vm = 0.0;
        }
        else
        {
            vm = stod(strLine.substr(42, 10)) * 10; 
        }

        outputFile<< fixed << setprecision(4);
        outputFile << "::-Vm\t" << vm << endl;

        for (size_t i = 0; i < 4; ++i)
            getline(file,strLine);

        positionStar = strLine.find("*******************************************************");

    } 


    
    outputFile << "\n" << "# minerals that using Landau theory" << "\n" << endl;

    for (size_t i = 0; i < 3; ++i)
        getline(file, strLine);

    positionStar = strLine.find("*******************************************************");
    while(positionStar == string::npos)
    {
        string name; 

        name = strLine.substr(0, 21);
        Utility::trim(name);

        outputFile << name << endl;

        for (size_t i = 0; i < 3; ++i)
            getline(file,strLine);

        double vm; 
        if(!Utility::isOnlyOneNumber(strLine.substr(42, 10)))
        {
            vm = 0.0;
        }
        else
        {
            vm = stod(strLine.substr(42, 10)) * 10; 
        }

        outputFile<< fixed << setprecision(4);
        outputFile << "::-Vm\t" << vm << endl;

        for (size_t i = 0; i < 5; ++i)
            getline(file,strLine);

        positionStar = strLine.find("*******************************************************");
    } 


    
    outputFile << "\n" << "# minerals that use Bragg-Williams theory" << "\n" << endl;

    for (size_t i = 0; i < 3; ++i)
        getline(file, strLine);

    positionStar = strLine.find("*******************************************************");
    while(positionStar == string::npos)
    {
        string name; 

        name = strLine.substr(0, 21);
        Utility::trim(name);

        outputFile << name << endl;

        for (size_t i = 0; i < 3; ++i)
            getline(file,strLine);

        double vm; 
        if(!Utility::isOnlyOneNumber(strLine.substr(42, 10)))
        {
            vm = 0.0;
        }
        else
        {
            vm = stod(strLine.substr(42, 10)) * 10; 
        }

        outputFile<< fixed << setprecision(4);
        outputFile << "::-Vm\t" << vm << endl;

        for (size_t i = 0; i < 5; ++i)
            getline(file,strLine);

        positionStar = strLine.find("*******************************************************");

    } 


    
    outputFile << "\n" << "# gases" << "\n" << endl;

    for (size_t i = 0; i < 6; ++i)
        getline(file, strLine);

    positionStar = strLine.find("*******************************************************");
    while(positionStar == string::npos)
    {
        string name; 

        name = strLine.substr(0, 21);
        Utility::trim(name);

        outputFile << name << endl;

        for (size_t i = 0; i < 6; ++i)
            getline(file,strLine);

        positionStar = strLine.find("*******************************************************");

    } 


    
    outputFile << "\n" << "# aqueous species" << "\n" << endl;

    for (size_t i = 0; i < 3; ++i)
        getline(file, strLine);

    while(!Utility::emptyOrBlanksTabs(strLine))
    {
        string name; 

        name = strLine.substr(0, 21);
        Utility::trim(name);

        outputFile << name << endl;

        for (size_t i = 0; i < 4; ++i)
            getline(file, strLine);

        double a1, a2, a3, a4, omega;


        if(!Utility::isOnlyOneNumber(strLine.substr(0, 15)))
        {
            a1 = 0.0;
        }
        else
        {
            a1 = stod(strLine.substr(0, 15)) / 4.184; 
        }


        if(!Utility::isOnlyOneNumber(strLine.substr(15, 11)))
        {
            a2 = 0.0;
        }
        else
        {
            a2 = stod(strLine.substr(15, 11)) / 4.184; 
        }


        if(!Utility::isOnlyOneNumber(strLine.substr(26, 11)))
        {
            a3 = 0.0;
        }
        else
        {
            a3 = stod(strLine.substr(26, 11)) / 4.184; 
        }


        if(!Utility::isOnlyOneNumber(strLine.substr(37, 11)))
        {
            a4 = 0.0;
        }
        else
        {
            a4 = stod(strLine.substr(37, 11)) / 4.184; 
        }


        getline(file,strLine);
        if(!Utility::isOnlyOneNumber(strLine.substr(26, 11)))
        {
            omega = 0.0;
        }
        else
        {
            omega = stod(strLine.substr(26, 11)) / 4.184; 
        }


        outputFile << fixed << setprecision(4);
        outputFile << "::-Vm\t" << a1 << "\t" << a2 << "\t" << a3 << "\t" << a4 << "\t" << omega << endl;

        for (size_t i = 0; i < 1; ++i)
            getline(file,strLine);

    } 


    file.close();
    outputFile.close();

} 

