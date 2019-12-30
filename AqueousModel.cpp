#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "ScreenInput.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "AqueousModel.h"
#include "Utility.h"
#include "Global.h"


using namespace std;



double AqueousModel::o2LogkArray[SIZEOFTEMPERATURES] = {0.0};
double AqueousModel::h2oDensityArray[SIZEOFTEMPERATURES] = {0.0};
double AqueousModel::pressureArray[SIZEOFTEMPERATURES] = {0.0};



void AqueousModel::setO2LogkArrayPtr(const double *arrayPtr, size_t sizeOfArray)
{
    for (size_t i = 0; i < sizeOfArray; ++i)
        o2LogkArray[i] = arrayPtr[i];
} 

void AqueousModel::setH2oDensityArrayPtr(const double *arrayPtr, size_t sizeOfArray)
{
    for (size_t i = 0; i < sizeOfArray; ++i)
        h2oDensityArray[i] = arrayPtr[i];
}

void AqueousModel::setPressureArrayPtr(const double *arrayPtr, size_t sizeOfArray)
{
    for (size_t i = 0; i < sizeOfArray; ++i)
        pressureArray[i] = arrayPtr[i];
}

double * AqueousModel::getO2LogkArrayPtr()
{
    return o2LogkArray;
} 

double * AqueousModel::getH2oDensityArrayPtr()
{
    return h2oDensityArray;
} 

double * AqueousModel::getPressureArrayPtr()
{
    return pressureArray;
} 


void AqueousModel::outputPhreeqcKeywords()
{
    ofstream &file = FileOutput::getPhreeqcFileReference();
    file << "LLNL_AQUEOUS_MODEL_PARAMETERS" << endl;
} 


void AqueousModel::outputTemperatures()
{
    size_t sizeOfT = SIZEOFTEMPERATURES; 
    double minT = ScreenInput::getMinTemperature(); 
    double maxT = ScreenInput::getMaxTemperature(); 
    double temperatureArray[sizeOfT]; 

    Utility::linspace(minT,maxT,sizeOfT,temperatureArray); 


    
    ofstream &file = FileOutput::getPhreeqcFileReference();
    file << "-temperatures" << endl;
    file << fixed << setprecision(4);

    for(size_t i = 0; i < sizeOfT; ++i)
    {
        file << setw(10) << temperatureArray[i];

        if ( (i + 1) % 4 == 0)
            file << endl;
    } 
    file << endl;

} 


void AqueousModel::spronsToDprons()
{
    ofstream commandFile("command.txt",ios::out); 
    commandFile << RELATIVEPATH << SPRONSFILE << endl;
    commandFile << "dprons.dat" << endl;

    commandFile.close();

    remove("dprons.dat");


    if(::systemName == WINDOWS)
    {
        system("cpronsbl.exe < command.txt > NUL");
    }
    else
    {
        system("./cpronsbl < command.txt > NUL");
    } 

} 


void AqueousModel::calculateO2LogK()
{

    
    ofstream commandFile("command.txt",ios::out);

    FlagOfDatabaseFormat flag = ScreenInput::getTypeOfDatabaseFormat();

    switch(flag)
    {
    case PHREEQCDAT100:
    {
        commandFile << "n\n"
                    "dprons.dat\n"
                    "3\n"
                    "2\n"
                    "1\n"
                    "1"
                    << endl;

        double minT = ScreenInput::getMinTemperature(); 
        double maxT = ScreenInput::getMaxTemperature(); 
        commandFile << minT
                    << ", "
                    << maxT
                    << ", "
                    << (maxT - minT) / (SIZEOFTEMPERATURES - 1)
                    << endl;

        break;
    }

    case LIQVAP:
    {
        commandFile << "n\n"
                    "dprons.dat\n"
                    "3\n"
                    "2\n"
                    "1\n"
                    "1"
                    << endl;

        double minT = ScreenInput::getMinTemperature(); 
        double maxT = ScreenInput::getMaxTemperature(); 
        commandFile << minT
                    << ", "
                    << maxT
                    << ", "
                    << (maxT - minT) / (SIZEOFTEMPERATURES - 1)
                    << endl;

        break;
    }

    case LLNLDATPSAT:
    {
        commandFile << "n\n"
                    "dprons.dat\n"
                    "3\n"
                    "2\n"
                    "1\n"
                    "1"
                    << endl;

        double minT = ScreenInput::getMinTemperature(); 
        double maxT = ScreenInput::getMaxTemperature(); 
        commandFile << minT
                    << ", "
                    << maxT
                    << ", "
                    << (maxT - minT) / (SIZEOFTEMPERATURES - 1)
                    << endl;

        break;
    }

    case CONSTANTPREESURE:
    {
        commandFile << "n\n"
                    "dprons.dat\n"
                    "3\n"
                    "1\n"
                    "2\n"
                    "n\n"
                    "1\n"
                    "1"
                    << endl;

        double pressure = ScreenInput::getConstantPressure(); 
        commandFile << pressure << ", " << pressure << ", 0" << endl;

        double minT = ScreenInput::getMinTemperature(); 
        double maxT = ScreenInput::getMaxTemperature(); 
        commandFile << minT
                    << ", "
                    << maxT
                    << ", "
                    << (maxT - minT) / (SIZEOFTEMPERATURES - 1)
                    << endl;

        break;
    }

    default:
        cout << "Should never reach here!" << endl;
        break;

    } 

    commandFile << "n\n"
                "2\n"
                "1\n"
                "O2,aq\n"
                "-2 H2O\n"
                " 1 O2,aq\n"
                "0\n"
                "y\n"
                "n\n"
                "logK.out\n"
                "1"
                << endl;

    commandFile.close();


    
    remove("logK.out");


    if(::systemName == WINDOWS)
    {
        system("SUPCRTBL.exe < command.txt > NUL");
    }
    else
    {
        system("./SUPCRTBL < command.txt > NUL");
    } 

} 


void AqueousModel::obtainO2LogkH2ODensityPressure()
{
    ifstream file("logK.out", ios::in); 

    string name; 
    double logkArray[SIZEOFTEMPERATURES] = {0.0}; 
    double densityArray[SIZEOFTEMPERATURES] = {0.0}; 
    double pressureArray[SIZEOFTEMPERATURES] = {0.0}; 

    FileInput objectFileInput;
    objectFileInput.getNameLogkDensityPressure(file, name, logkArray, densityArray, pressureArray, SIZEOFTEMPERATURES);

    setO2LogkArrayPtr(logkArray, SIZEOFTEMPERATURES);
    setH2oDensityArrayPtr(densityArray, SIZEOFTEMPERATURES);
    setPressureArrayPtr(pressureArray, SIZEOFTEMPERATURES);

    file.close();
} 


double AqueousModel::calculateDielectricConstant(double T, double density)
{
    const double a1 = 0.1470333593e2;
    const double a2 = 0.2128462733e3;
    const double a3 = -0.1154445173e3;
    const double a4 = 0.1955210915e2;
    const double a5 = -0.8330347980e2;
    const double a6 = 0.3213240048e2;
    const double a7 = -0.6694098645e1;
    const double a8 = -0.378620204e2;
    const double a9 = 0.6887359646e2;
    const double a10 = -0.2729401652e2;

    double T_r = T/298.15;
    double density_r = density / 1;
    double dielectricConstant;

    dielectricConstant = 1 + (a1 / T_r) * density_r +
        (a2 / T_r + a3 + a4 * T_r) * pow(density_r, 2.0) +
        (a5 / T_r + a6 * T_r + a7 * pow(T_r, 2.0)) * pow(density_r, 3.0) +
        (a8 / pow(T_r, 2.0) + a9 / T_r + a10) * pow(density_r, 4.0);

    return dielectricConstant;
} 


double AqueousModel::calculateDebyHuckelA(double T, double density, double dielectricConstant)
{
    double debyHuckelA;

    debyHuckelA = 1.8246e6 * pow(density, 0.5) / pow( (dielectricConstant * T), 1.5);

    return debyHuckelA;
} 


double AqueousModel::calculateDebyHuckelB(double T, double density, double dielectricConstant)
{
    double debyHuckelB;

    debyHuckelB = 50.29e8 * pow(density, 0.5) / pow( (dielectricConstant * T), 0.5) / 1e8;

    return debyHuckelB;
} 


void AqueousModel::outputebyHuckelAB()
{
    const size_t sizeOfPoints = SIZEOFTEMPERATURES;

    double *densityArrayPtr = getH2oDensityArrayPtr(); 

    double minT = ScreenInput::getMinTemperature(); 
    double maxT = ScreenInput::getMaxTemperature(); 
    double temperatureArray[sizeOfPoints]; 
    Utility::linspace(minT,maxT,sizeOfPoints,temperatureArray); 

    for (size_t i = 0; i < sizeOfPoints; ++i) 
        temperatureArray[i] = temperatureArray[i] + 273.15;


    
    double dielectricConstantArray[sizeOfPoints];

    for (size_t i = 0; i < sizeOfPoints; ++i)
        dielectricConstantArray[i] = calculateDielectricConstant(temperatureArray[i], densityArrayPtr[i]);


    
    double debyeHuckelAArray[sizeOfPoints];

    for (size_t i = 0; i < sizeOfPoints; ++i)
        debyeHuckelAArray[i] = calculateDebyHuckelA(temperatureArray[i], densityArrayPtr[i], dielectricConstantArray[i]);


    
    double debyeHuckelBArray[sizeOfPoints];

    for (size_t i = 0; i < sizeOfPoints; ++i)
        debyeHuckelBArray[i] = calculateDebyHuckelB(temperatureArray[i], densityArrayPtr[i], dielectricConstantArray[i]);


    
    ofstream &file = FileOutput::getPhreeqcFileReference();
    file << "-dh_a" << endl;
    file << fixed << setprecision(4);

    for(size_t i = 0; i < sizeOfPoints; ++i)
    {
        file << setw(10) << debyeHuckelAArray[i];

        if ( (i + 1) % 4 == 0)
            file << endl;
    } 
    file << endl;


    
    file << "-dh_b" << endl;
    file << fixed << setprecision(4);

    for(size_t i = 0; i < sizeOfPoints; ++i)
    {
        file << setw(10) << debyeHuckelBArray[i];

        if ( (i + 1) % 4 == 0)
            file << endl;
    } 
    file << endl;

} 


void AqueousModel::outputBdot()
{
    ifstream file(RELATIVEPATH + BDOT_TP, ios::in); 
    size_t numRow = 0; 


    
    FileInput objectFileInput;
    while(!file.eof())
    {
        string str;
        objectFileInput.getLineExcludingComments(file, str);

        if(!str.empty())
            ++numRow;
    } 
    file.clear();
    file.seekg(0);


    
    


    
    double xymeasureArray[2 * numRow];
    double zmeasureArray[numRow];

    istringstream inputString;
    size_t position = 0;

    while(!file.eof())
    {
        string str;
        objectFileInput.getLineExcludingComments(file, str);

        if(!str.empty())
        {
            inputString.clear();
            inputString.str(str);

            inputString >> xymeasureArray[2 * position];
            inputString >> xymeasureArray[2 * position + 1];
            inputString >> zmeasureArray[position];

            ++position;

        } 

    } 
    file.close();

    size_t sizeOfT = SIZEOFTEMPERATURES; 
    double minT = ScreenInput::getMinTemperature(); 
    double maxT = ScreenInput::getMaxTemperature(); 
    double temperatureArray[sizeOfT]; 

    Utility::linspace(minT,maxT,sizeOfT,temperatureArray); 

    double *pressureArrayPtr = getPressureArrayPtr(); 
    double xyInterpoArray[2 * sizeOfT]; 
    double bdotArray[sizeOfT]; 

    for(size_t i = 0; i < sizeOfT; ++i)
    {
        xyInterpoArray[2 * i] = temperatureArray[i];
        xyInterpoArray[2 * i + 1] = pressureArrayPtr[i];
    } 

    Utility::scatterInterpolation2DSpaceLinear(numRow, xymeasureArray,
                                               zmeasureArray, sizeOfT,
                                               xyInterpoArray, bdotArray); 





    for(size_t i = 0; i < sizeOfT; ++i)
        bdotArray[i] = bdotArray[i] / 100; 


    
    ofstream &fileRefrence = FileOutput::getPhreeqcFileReference();
    fileRefrence << "-bdot" << endl;
    fileRefrence << fixed << setprecision(4);

    for(size_t i = 0; i < sizeOfT; ++i)
    {
        fileRefrence << setw(10) << bdotArray[i];

        if ( (i + 1) % 4 == 0)
            fileRefrence << endl;
    } 
    fileRefrence << endl;


} 


void AqueousModel::outputCo2_coefs()
{
    ofstream &file = FileOutput::getPhreeqcFileReference();
    file << CO2COEFS << endl;
    file << endl;
} 




void AqueousModel::runInstructs()
{
    spronsToDprons();
    calculateO2LogK();
    obtainO2LogkH2ODensityPressure();


    FlagOfDatabaseFormat flag = ScreenInput::getTypeOfDatabaseFormat();
    switch(flag)
    {
    case PHREEQCDAT100:
    case LIQVAP:
        break;

    case LLNLDATPSAT:
    case CONSTANTPREESURE:
        outputPhreeqcKeywords();
        outputTemperatures();
        outputebyHuckelAB();
        outputBdot();
        outputCo2_coefs();

        break;

    default:
        cout << "Should never reach here!" << endl;
        break;

    } 

} 
