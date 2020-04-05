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
            "\t1 -- phreeqc.dat format  (0.01 - 200 oC up to 1000 bar)\n"
            "\t2 -- llnl.dat format (0.01 - 100 oC at 1 bar and 100 - 300 oC at Psat)\n"
            "\t3 -- llnl.dat format (0.01 - 1000 oC at a constant P up to 5000 bar)\n"
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
























































































void ScreenInput::readMinTemperatureFromScreen_Psat(double lowerBoundTemp)
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


void ScreenInput::readMaxTemperatureFromScreen_Psat(double lowerBoundTemp, double upperBoundTemp)
{
    cout << "MAX temperature (<= " << upperBoundTemp << " oC):" << endl;
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


void ScreenInput::readMinTemperatureFromScreen_isobaric(double lowerBoundTemp)
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


void ScreenInput::readMaxTemperatureFromScreen_isobaric(double lowerBoundTemp, double upperBoundTemp)
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


void ScreenInput::readConstantPressure_isobaric(double lowerBoundPressure, double upperBoundPressure)
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
        cout << "Please enter a temperature range(oC) for calculating log K\n" << endl;

        readMinTemperatureFromScreen_Psat(0.01);

        double lowerBoundTemp = max(getMinTemperature(),100.0);
        readMaxTemperatureFromScreen_Psat(lowerBoundTemp, 200);

        break;
    }

    case LLNLDATPSAT:
    {
        cout << "Please enter a temperature range(oC) for calculating log K\n" << endl;

        readMinTemperatureFromScreen_Psat(0.01);

        double lowerBoundTemp = max(getMinTemperature(),100.0);
        readMaxTemperatureFromScreen_Psat(lowerBoundTemp, 300);

        break;
    }

    case CONSTANTPREESURE:
    {
        cout << "Please enter the T-P range (oC and bar) for calculating log K (must be in the applicable T-P range of SUPCRTBL. See Fig. 1 in Zhang et al. in review)\n" << endl;

        readMinTemperatureFromScreen_isobaric(0.01);

        double lowerBoundTemp = getMinTemperature();
        readMaxTemperatureFromScreen_isobaric(lowerBoundTemp, 1000);

        readConstantPressure_isobaric(0.01, 5000);

        break;
    }

    default:
    {
        cout << "Should never reach here!" << endl;
        break;
    }

    } 

} 



