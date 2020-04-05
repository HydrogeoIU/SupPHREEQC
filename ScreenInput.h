#ifndef SCREENINPUT_H
#define SCREENINPUT_H

#include "Global.h"

class ScreenInput
{
public:

    
    static void setTypeOfDatabaseFormat(FlagOfDatabaseFormat);
    static void setMinTemperature(double);
    static void setMaxTemperature(double);
    static void setConstantPressure(double);

    static FlagOfDatabaseFormat getTypeOfDatabaseFormat();
    static double getMinTemperature();
    static double getMaxTemperature();
    static double getConstantPressure();


    
    void readTypeOfDatabaseFromScreen();












    
    void readMinTemperatureFromScreen_Psat(double lowerBoundTemp);


    
    void readMaxTemperatureFromScreen_Psat(double lowerBoundTemp, double upperBoundTemp);


    
    void readMinTemperatureFromScreen_isobaric(double lowerBoundTemp);


    
    void readMaxTemperatureFromScreen_isobaric(double lowerBoundTemp, double upperBoundTemp);


    
    void readConstantPressure_isobaric(double lowerBoundPressure, double upperBoundPressure);


    
    void runInstructs();


private:
    static FlagOfDatabaseFormat typeOfDatabaseFormat;  
    static double minTemperature; 
    static double maxTemperature; 
    static double constantPressure; 

}; 

#endif 
