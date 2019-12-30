#ifndef FILEOUTPUT_H_INCLUDED
#define FILEOUTPUT_H_INCLUDED

#include <fstream>

#include "Global.h"

class FileOutput
{
public:
    
    static std::ofstream & getPhreeqcFileReference();


    
    void outputTypeInformOfDatabaseFormat();


    
    void outputKineticsScript();


    
    void closePhreeqcFile();


    
    void changeDatabaseFileName();


private:
    static std::ofstream phreeqcFile; 

}; 

#endif 
