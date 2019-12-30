#ifndef UTILITY_H
#define UTILITY_H

#include <cstddef>
#include <string>

namespace Utility
{
    
    void linspace(double a, double b, size_t n, double *yArray);


    
    void trim(std::string &s);


    
    bool emptyOrBlanksTabs(const std::string &s);


    
    bool isOnlyNumbers(std::string str);


    
    bool isOnlyOneNumber(std::string str);


    


    
    void scatterInterpolation2DSpaceLinear(size_t numNode, double *xyNodeArrayPtr,
                                           double *zNodeArrayPtr, size_t numInterpo,
                                           double *xyInterpoArrayPtr, double *zInterpoArrayPtr);


    


    
    void estimateLogkAnalyticFormula(const double *tMeasureArrayPtr, const double *logkMeasureArrayPtr,
                                     size_t numMeasure, double *coefficientsArrayPtr,
                                     double &maxRelativeError, double &corresponAbsoluteError,
                                     double &correspondingLogk, double &corresponTemporature);


    
    void obtainVmFromSUPCRTThenStore();



}; 

#endif 
