#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "ScreenInput.h"
#include "ScreenOutput.h"
#include "FileOutput.h"
#include "AqueousModel.h"
#include "SpeciesMineralGas.h"
#include "Utility.h"
#include "Global.h"

#include <iomanip>
using namespace std;


int main()
{

    
    ScreenOutput objectScreenOutput;
    objectScreenOutput.outputCodeTitleToScreen();


    
    ScreenInput objectScreenInput;
    objectScreenInput.runInstructs();


    
    FileOutput objectFileOutput;
    objectFileOutput.outputTypeInformOfDatabaseFormat();


    
    AqueousModel objectAqueousModel;
    objectAqueousModel.runInstructs();


    
    SpeciesMineralGas objectSpeciesMineralGas;
    objectSpeciesMineralGas.runInstructs();


    
    objectFileOutput.closePhreeqcFile();
    objectFileOutput.changeDatabaseFileName();


    cout << "Done!" << endl;

    if(::systemName == WINDOWS)
    {
        system("pause");
    }
    else
    {
    } 


} 
