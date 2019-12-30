#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <iomanip>

#include "ScreenInput.h"
#include "FileInput.h"
#include "Utility.h"
#include "Global.h"


using namespace std;


bool FileInput::getNameLogkDensityPressure(ifstream &file, string &name,
                                        double *logkArrayPtr, double *densityArrayPtr,
                                        double *pressureArrayPtr, size_t numOfPoints)
{   
    if(file.eof())
        return false;


    
    string temporaryString; 
    GoodBad flag = BAD; 

    while(!file.eof())
    {
        getline(file, temporaryString);

        size_t index = temporaryString.find("************* REACTION");

        if(index != string::npos)
        {
            flag = GOOD;
            break;
        } 

    } 

    if(flag == BAD)
        return false;


    
    while(!file.eof())
    {
        getline(file, temporaryString);

        if (Utility::emptyOrBlanksTabs(temporaryString))
            continue;

        Utility::trim(temporaryString); 
        name = temporaryString;
        break;

    } 


    
    while(!file.eof())
    {
        getline(file, temporaryString);

        size_t index = temporaryString.find("STANDARD STATE PROPERTIES OF THE REACTION AT ELEVATED TEMPERATURES AND PRESSURES");

        if(index != string::npos)
        {
            break;
        } 

    } 

    for(size_t i = 1; i <= 6; ++i)
         getline(file, temporaryString);

    for (size_t i = 0; i < numOfPoints; i = i)
    {
        getline(file, temporaryString);

        if(Utility::emptyOrBlanksTabs(temporaryString))
        {

        }
        else
        {
            string subString;
            subString = temporaryString.substr(24 - 1, 12);
            if(!Utility::isOnlyOneNumber(subString))
            {
                throw domain_error("H2O density is not calculated correctly by SUPCRT!");
            }
            densityArrayPtr[i] = stod(subString); 

            subString = temporaryString.substr(36 - 1, 12);
            if(!Utility::isOnlyOneNumber(subString))
            {
                throw domain_error("logK is not calculated correctly by SUPCRT!");
            }
            logkArrayPtr[i] = stod(subString);


            FlagOfDatabaseFormat flag = ScreenInput::getTypeOfDatabaseFormat();
            switch(flag)
            {
            case PHREEQCDAT100:
            case LIQVAP:
            case LLNLDATPSAT:
                subString = temporaryString.substr(12 - 1, 12);
                break;

            case CONSTANTPREESURE:
                subString = temporaryString.substr(1 - 1, 11);
                break;

            default:
                cout << "Should never reach here!" << endl;
                break;

            } 


            if(!Utility::isOnlyOneNumber(subString))
            {
                throw domain_error("Pressure is not calculated correctly!");
            }
            pressureArrayPtr[i] = stod(subString);

            ++i;

        } 

    } 

    return true;

} 


void FileInput::getLineExcludingComments(ifstream &file, string &str)
{
    GoodBad flag = BAD;


    if(file.eof()) 
    {
        str = "";
        return;
    } 


    while(!file.eof())
    {
        getline(file, str);

        if (str.empty()) 
            continue;

        if (str.find_first_not_of(" \t\r") == string::npos) 
            continue;

        size_t position = str.find_first_not_of(" \t\r");
        if( (str[position]) == '#') 
            continue;

        flag = GOOD;
        break;
    } 


    if(flag == GOOD)
    {
        size_t position = str.find_first_of("#");

        if(position == string::npos)
        {

        }
        else
        {
            str.erase(position);  
        } 

        return;

    }
    else
    {
        str = "";
        return; 
    } 

} 


void FileInput::getLineExcludingComments_ButDoulbeColon(ifstream &file, string & str)
{
    GoodBad flag = BAD;


    if(file.eof()) 
    {
        str = "";
        return;
    } 


    while(!file.eof())
    {
        getline(file, str);

        if (str.empty()) 
            continue;

        if (str.find_first_not_of(" \t\r") == string::npos) 
            continue;

        size_t position = str.find_first_not_of(" \t\r");
        if( (str[position]) == '#') 
            continue;

        flag = GOOD;
        break;
    } 


    if(flag == GOOD)
    {
        size_t position = str.find_first_not_of(" \t\r");
        string subString = str.substr(position, 2);
        if( subString == "::") 
            return;


        position = str.find_first_of("#");

        if(position == string::npos)
        {

        }
        else
        {
            str.erase(position);  
        } 

        return;

    }
    else
    {
        str = "";
        return; 
    } 

} 


void FileInput::excludeAqOfLine(string &str)
{
    if(str.find(",aq") == string::npos)
    {
        return;
    } 

    size_t position = str.find(",aq");
    while(position != string::npos)
    {
        str.replace(position, 3, "   ");
        position = str.find(",aq", position + 3);
    } 

} 


double FileInput::getElementValenceThenErase(string &str)
{
    size_t position = str.find("::valence:");
    string subString;
    subString = str.substr(position + 10);

    stringstream strStream;
    strStream << subString;

    double valence;
    strStream >> valence;

    str.erase(position);

    return valence;

} 


bool FileInput::getNameReactionFormulaFromRXN(ifstream &file, string &name, string &reactionFormula, size_t &type)
{
    
    if(file.eof())
        return false;


    
    string strLine;
    getline(file, strLine);
    Utility::trim(strLine);
    if(strLine.empty())
        return false;


    
    if(strLine.find("Line 1:") != string::npos)
    {
        for (size_t i = 0; i < 16; ++i)
            getline(file, strLine);
    }


    
    stringstream strStream;
    strStream.str(strLine);
    strStream >> name;


    
    size_t numMineralType;
    size_t numSpeciesType;
    size_t numGasType;
    size_t numWaterType;

    getline(file, strLine);
    strStream.clear();
    strStream.str(strLine);

    strStream >> numMineralType;
    strStream >> numSpeciesType;
    strStream >> numGasType;
    strStream >> numWaterType;


    
    unordered_map<string, double> leftRXNFormula;
    unordered_map<string, double> rightRXNFormula;

    for (size_t i = 0; i < numMineralType + numSpeciesType + numGasType + numWaterType; ++i)
    {
        getline(file, strLine);
        strStream.clear();
        strStream.str(strLine);

        double coefficient;
        string speciesPhase;
        strStream >> coefficient;
        strStream >> speciesPhase;

        if(coefficient < 0)
        {
            leftRXNFormula.insert({speciesPhase, -1 * coefficient});
        }
        else
        {
            rightRXNFormula.insert({speciesPhase, coefficient});
        } 

    } 

    strStream.clear();
    strStream.str("");
    strStream.setf(ios::fixed);
    strStream.precision(4);
    if(numGasType > 0)
    {
        type = 3;

        leftRXNFormula.erase(name);

        strStream << name;

        for(auto iterm : leftRXNFormula)
            strStream << " + " << iterm.second << " " << iterm.first;

        strStream << " = ";

        auto iter = rightRXNFormula.begin();
        strStream << iter->second << " " << iter->first;

        ++iter;
        for( ; iter != rightRXNFormula.end(); ++iter)
            strStream << " + " << iter->second << " " << iter->first;

    }
    else if(numMineralType > 0)
    {
        type = 1;

        leftRXNFormula.erase(name);

        strStream << name;

        for(auto iterm : leftRXNFormula)
            strStream << " + " << iterm.second << " " << iterm.first;

        strStream << " = ";

        auto iter = rightRXNFormula.begin();
        strStream << iter->second << " " << iter->first;

        ++iter;
        for( ; iter != rightRXNFormula.end(); ++iter)
            strStream << " + " << iter->second << " " << iter->first;
    }
    else
    {
        type = 2;

        rightRXNFormula.erase(name);

        auto iter = leftRXNFormula.begin();
        strStream << iter->second << " " << iter->first;

        ++iter;
        for( ; iter != leftRXNFormula.end(); ++iter)
            strStream << " + " << iter->second << " " << iter->first;

        strStream << " = " << name;

        for(auto iterm : rightRXNFormula)
            strStream << " + " << iterm.second << " " << iterm.first;
    } 

    reactionFormula = strStream.str();


    
    excludeAqOfLine(name); 
    Utility::trim(name);

    excludeAqOfLine(reactionFormula); 
    Utility::trim(reactionFormula);


    getline(file, strLine);
    return true;

} 



