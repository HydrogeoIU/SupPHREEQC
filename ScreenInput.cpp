#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "ScreenInput.h"
#include "FileOutput.h"

using namespace std;


FlagOfDatabaseFormat ScreenInput::typeOfDatabaseFormat = LIQVAP;
double ScreenInput::minTemperature = 0.01;
double ScreenInput::maxTemperature = 300;
double ScreenInput::constantPressure = 1;



void ScreenInput::setTypeOfDatabaseFormat(FlagOfDatabaseFormat typeOfDatabaseFormat_)
{
    typeOfDatabaseFormat = typeOfDatabaseFormat_;
} 

void ScreenInput::setMinTemperature(double minTemperature_)
{
    minTemperature = minTemperature_;
} 

void ScreenInput::setMaxTemperature(double maxTemperature_)
{
    maxTemperature = maxTemperature_;
} 

void ScreenInput::setConstantPressure(double constantPressure_)
{
    constantPressure = constantPressure_;
} 

FlagOfDatabaseFormat ScreenInput::getTypeOfDatabaseFormat()
{
    return typeOfDatabaseFormat;
} 

double ScreenInput::getMinTemperature()
{
    return minTemperature;
} 

double ScreenInput::getMaxTemperature()
{
    return maxTemperature;
} 

double ScreenInput::getConstantPressure()
{
    return constantPressure;
} 


void ScreenInput::readTypeOfDatabaseFromScreen()
{





    cout << "Please choose the database for PHREEQC: \n\n"
            "\t1 -- diagenesis.dat (0.01 - 200 oC up to 1000 bar)\n"
            "\t2 -- geothermal.dat (0.01 - 100 oC at 1 bar and 100 - 300 oC at Psat)\n"
            "\t3 -- bl.dat (0.01 - 1000 oC at a constant P up to 5000 bar)\n"
            << endl;


    
    cout << "?";
    int flag; 
    cin >> flag;
    cin.clear();
    cin.sync();
    cout << endl;

    while(flag != 1 && flag != 2 && flag != 3)
    {
        cout << "?";
        cin >> flag;
        cin.clear();
        cin.sync();
        cout << endl;
    } 

    switch(flag)
    {
    case 1:
        break;
    case 2:
        setTypeOfDatabaseFormat(LLNLDATPSAT);
        break;
    case 3:
        setTypeOfDatabaseFormat(CONSTANTPREESURE);
        break;
    default:
        cout << "Should never reach here!" << endl;
        break;
    } 

} 
























































































void ScreenInput::readMinTemperatureFromScreen(double lowerBoundTemp)
{
    cout << "MIN temperature (>= 0.01 oC):" << endl;
    cout << "?";
    double temperature;
    cin >> temperature; 
    cin.clear();
    cin.sync();
    cout << endl;

    cout << fixed << setprecision(2);
    while(temperature < lowerBoundTemp)
    {
        cout << "MIN temperature must be >= " << lowerBoundTemp << ".\n";
        cout << "?";
        cin >> temperature;
        cin.clear();
        cin.sync();
        cout << endl;
    } 
    setMinTemperature(temperature) ;

} 


void ScreenInput::readMaxTemperatureFromScreen(double lowerBoundTemp, double upperBoundTemp)
{
    cout << "MAX temperature (<= 1000 oC):" << endl;
    cout << "?";
    double temperature;
    cin >> temperature; 
    cin.clear();
    cin.sync();
    cout << endl;

    cout << fixed << setprecision(2);
    while(temperature <= lowerBoundTemp || temperature > upperBoundTemp)
    {
        cout << "MAX temperature must be > " << lowerBoundTemp << " and <= " << upperBoundTemp << ".\n";
        cout << "?";
        cin >> temperature;
        cin.clear();
        cin.sync();
        cout << endl;
    } 

    setMaxTemperature(temperature);

} 


void ScreenInput::readConstantPressure(double lowerBoundPressure, double upperBoundPressure)
{
    cout << "Please enter a constant pressure (<= 5000 bar):" << endl;
    cout << "?";
    double pressure;
    cin >> pressure; 
    cin.clear();
    cin.sync();
    cout << endl;

    while(pressure < lowerBoundPressure || pressure > upperBoundPressure)
    {
        cout << "Pressure must be >= " << lowerBoundPressure << " and <= " << upperBoundPressure << ".\n";
        cout << "?";
        cin >> pressure;
        cin.clear();
        cin.sync();
        cout << endl;
    } 
    setConstantPressure(pressure);
} 


void ScreenInput::runInstructs()
{
    readTypeOfDatabaseFromScreen();

    FlagOfDatabaseFormat flag = getTypeOfDatabaseFormat(); 
    switch(flag)
    {





    case LIQVAP:
    {







        setMinTemperature(0.01);
        setMaxTemperature(200);

        break;
    }

    case LLNLDATPSAT:
        setMinTemperature(0.01);
        setMaxTemperature(300);
        break;

    case CONSTANTPREESURE:
    {
        cout << "Please enter the T-P range (oC and bar) for calculating log K (must be in the applicable T-P range of SUPCRTBL. See Fig. 1 in Zhang et al. in review)\n" << endl;

        readMinTemperatureFromScreen(0.01);

        double lowerBoundTemp = getMinTemperature();
        readMaxTemperatureFromScreen(lowerBoundTemp, 1000);

        readConstantPressure(0.01, 5000);

        break;
    }

    default:
    {
        cout << "Should never reach here!" << endl;
        break;
    }

    } 

} 



