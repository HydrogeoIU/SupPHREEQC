#ifndef FILEINPUT_H
#define FILEINPUT_H

#include <fstream>
#include <string>

class FileInput
{
public:

    
    bool getNameLogkDensityPressure(std::ifstream &file, std::string &name,
                                   double *logkArrayPtr, double *densityArrayPtr,
                                   double *pressureArrayPtr, size_t numOfPoints);


    
    void getLineExcludingComments(std::ifstream &file, std::string & str);


    
    void getLineExcludingComments_ButDoulbeColon(std::ifstream &file, std::string & str);


    
    void excludeAqOfLine(std::string &str);


    
    double getElementValenceThenErase(std::string &str);


    
    bool getNameReactionFormulaFromRXN(std::ifstream &file, std::string &name,
                                              std::string &reactionFormula, size_t &type);



}; 

#endif 
