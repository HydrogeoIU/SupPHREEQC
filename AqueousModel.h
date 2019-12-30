#ifndef AQUEOUSMODEL_H
#define AQUEOUSMODEL_H

#include "Global.h"

class AqueousModel
{
public:

    
    static void setO2LogkArrayPtr(const double *arrayPtr, size_t sizeOfArray);
    static void setH2oDensityArrayPtr(const double *arrayPtr, size_t sizeOfArray);
    static void setPressureArrayPtr(const double *arrayPtr, size_t sizeOfArray);
    static double * getO2LogkArrayPtr();
    static double * getH2oDensityArrayPtr();
    static double * getPressureArrayPtr();


    
    void outputPhreeqcKeywords();


    
    void outputTemperatures();


    
    void spronsToDprons();


    
    void calculateO2LogK();


    
    void obtainO2LogkH2ODensityPressure();


    
    double calculateDielectricConstant(double T, double density);


    
    double calculateDebyHuckelA(double T, double density, double dielectricConstant);


    
    double calculateDebyHuckelB(double T, double density, double dielectricConstant);


    
    void outputebyHuckelAB();


    
    void outputBdot();


    
    void outputCo2_coefs();


    
    void runInstructs();


private:

    static double o2LogkArray[SIZEOFTEMPERATURES]; 
    static double h2oDensityArray[SIZEOFTEMPERATURES]; 
    static double pressureArray[SIZEOFTEMPERATURES]; 


}; 

#endif 
