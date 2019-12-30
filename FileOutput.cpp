#include <stdio.h>
#include <iostream>

#include "ScreenInput.h"
#include "FileOutput.h"
#include "Global.h"


using namespace std;



ofstream FileOutput::phreeqcFile(NAMEOFPHREEQCDATABASE, ios::out);



ofstream & FileOutput::getPhreeqcFileReference()
{
    return phreeqcFile;
} 


void FileOutput::outputTypeInformOfDatabaseFormat()
{
    ofstream & file = getPhreeqcFileReference();
    FlagOfDatabaseFormat flag = ScreenInput::getTypeOfDatabaseFormat();

    file << "# Developed by Guanru Zhang, Peng Lu, Yilun Zhang, Kevin Tu, and Chen Zhu\n"
            "# All inquiries should be directed to supcrt@indiana.edu\n"
            "# On 12/27/2019\n"
            "# Citation: Zhang GR, Lu P, Zhang YL, Tu K, Zhu C (in review) SUPPHREEQC: A program\n"
            "#    to generate customized PHREEQC thermodynamic databases from SUPCRTBL and\n"
            "#    to extend calculations to elevated pressures and temperatures. Computers &\n"
            "#    Geosciences\n"
            "# Visit https://hydrogeochem.earth.indiana.edu for updates\n\n\n";

    switch(flag)
    {

    case PHREEQCDAT100:
    {
        file << "# PHREEQC.DAT framework (LogKs are along Psat and molar volumes are included)." << endl;
        file << "# Temperature(oC): 0.01 - 100." << endl;
        file << "# Pressure is corrected up to 1000 bar by molar volumes." << endl;

        break;
    }

    case LIQVAP:
    {
        file << "# PHREEQC.DAT framework (LogKs are along Psat and molar volumes are included)." << endl;
        file << "# Temperature(oC): " << ScreenInput::getMinTemperature() << " - " << ScreenInput::getMaxTemperature() << "." << endl;
        file << "# Pressure is corrected up to 1000 bar by molar volumes." << endl;

        break;
    }

    case LLNLDATPSAT:
    {
        file << "# LLNL.DAT framework (LogKs are along Psat)." << endl;
        file << "# Temperature(oC): 0.01 - 300" << endl;

        break;
    }


    case CONSTANTPREESURE:
    {
        file << "# LLNL.DAT framework (LogKs are under varied T and constant P in one-phase region)." << endl;
        file << "# Temperature(oC): " << ScreenInput::getMinTemperature() << " - " << ScreenInput::getMaxTemperature() << endl;
        file << "# Constant pressure(bar): " << ScreenInput::getConstantPressure() << endl;

        break;
    }

    default:
    {
        cout << "Should never reach here!" << endl;
        break;
    }

    } 

    file << "# It's users' responsibility to mind the valid upper limits of temperature and pressure for all data used.";
    file << "\n" << "\n" << endl;

} 


void FileOutput::outputKineticsScript()
{
    ifstream file(RELATIVEPATH + KINETICSSCRIPT_F, ios::in); 
    ofstream &outputFile = getPhreeqcFileReference(); 

    outputFile << endl << endl;

    while(!file.eof())
    {
        string strLine;
        getline(file, strLine);

        outputFile << strLine << endl;
    }

    file.close();
} 


void FileOutput::closePhreeqcFile()
{
    ofstream &outputFile = getPhreeqcFileReference(); 
    outputFile.close();
} 


void FileOutput::changeDatabaseFileName()
{
    remove("PHREEQC-bl.DAT");
    remove("LLNL-bl.DAT");

    remove("diagenesis.dat");
    remove("geothermal.dat");
    remove("bl.dat");

    string oldNameStr = NAMEOFPHREEQCDATABASE;
    string newNameStr = "";

    FlagOfDatabaseFormat flag = ScreenInput::getTypeOfDatabaseFormat();
    switch(flag)
    {




    case LIQVAP:
        newNameStr = "diagenesis.dat";
        break;

    case LLNLDATPSAT:
        newNameStr = "geothermal.dat";
        break;

    case CONSTANTPREESURE:
        newNameStr = "bl.dat";
        break;

    default:
        cout << "Should never reach here!" << endl;
        break;

    } 

    rename(oldNameStr.c_str(), newNameStr.c_str());

}



