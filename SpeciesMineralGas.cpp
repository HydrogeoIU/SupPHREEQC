#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>

#include "ScreenInput.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "AqueousModel.h"
#include "SpeciesMineralGas.h"
#include "Utility.h"
#include "Global.h"


using namespace std;


void SpeciesMineralGas::getOutputSolutionMasterSpecies()
{
    ifstream file(RELATIVEPATH + SOLUTION_MASTER_SPECIES_F, ios::in); 
    ofstream &outpuFile = FileOutput::getPhreeqcFileReference(); 

    outpuFile << "SOLUTION_MASTER_SPECIES" << endl; 


    while(!file.eof())
    {

        
        string strLine;
        FileInput objectFileInput;
        objectFileInput.getLineExcludingComments(file, strLine); 
        if(strLine.empty())
            break;

        objectFileInput.excludeAqOfLine(strLine); 

        stringstream strStream;
        strStream << strLine;

        string shortStr;
        strStream >> shortStr;

        string primaryMasterName = ""; 
        string secondaryMasterName = ""; 
        PrimaryMasterChemicalData primaryMasterDataStruct; 
        SecondaryMasterChemicalData secondaryMasterDataStruct; 
        if(shortStr.find("(") == string::npos)
        {
            if(shortStr == NAMEALKALINITYELEMENT) 
            {
                objectFileInput.getElementValenceThenErase(strLine); 
                outpuFile << strLine << endl; 
                continue;
            }

            primaryMasterDataStruct.element = shortStr;  

            strStream >> shortStr;
            primaryMasterName = shortStr; 
            primaryMasterDataStruct.chemicalFormula = primaryMasterName; 

            primaryMasterDataStruct.valence = objectFileInput.getElementValenceThenErase(strLine); 

            primaryMasterDataStruct.reactionFormula = primaryMasterName + " = " + primaryMasterName; 

            primaryMasterMap.insert({primaryMasterName, primaryMasterDataStruct}); 
        }
        else
        {
            secondaryMasterDataStruct.element = shortStr.substr(0, shortStr.find("(")); 

            size_t position1 = shortStr.find("(");
            size_t position2 = shortStr.find(")");
            secondaryMasterDataStruct.valence = stod( shortStr.substr(position1 + 1,  position2 - position1 - 1) ); 

            strStream >> shortStr;
            secondaryMasterName = shortStr; 
            secondaryMasterDataStruct.chemicalFormula = secondaryMasterName; 

            secondaryMasterMap.insert({secondaryMasterName, secondaryMasterDataStruct}); 
        } 


        
        outpuFile << strLine << endl;


    } 

    outpuFile << "\n" << endl;

    file.close();


    
    for (auto iterm : primaryMasterMap)
    {
        for (auto iter = secondaryMasterMap.begin(); iter != secondaryMasterMap.end();  )
        {

            if(iterm.first == iter->first)
            {
                secondaryMasterMap.erase(iter++);
                break;
            }
            else
            {
                ++iter;
            }

        }

    } 

} 


void SpeciesMineralGas::readExtraIdentifiers()
{
    string ExtraIdentifierFileName; 
    FlagOfDatabaseFormat flag = ScreenInput::getTypeOfDatabaseFormat(); 

    switch(flag)
    {

    case PHREEQCDAT100:
    case LIQVAP:
    {
        ExtraIdentifierFileName = EXTRAIDENTIFIER_PHREEQCDAT_F;
        break;
    }

    case LLNLDATPSAT:
    case CONSTANTPREESURE:
    {
        ExtraIdentifierFileName = EXTRAIDENTIFIER_F;
        break;
    }

    default:
    {
        cout << "Should never reach here!" << endl;
        break;
    }

    } 

    ifstream file(RELATIVEPATH + ExtraIdentifierFileName, ios::in); 


    string speciesName; 
    IdentifierData identiferDataTempora; 


    while(!file.eof())
    {
        string strLine;
        FileInput objectFileInput;
        objectFileInput.getLineExcludingComments_ButDoulbeColon(file, strLine); 
        if(strLine.empty())
            break;

        objectFileInput.excludeAqOfLine(strLine); 

        size_t position = strLine.find_first_not_of(" \t\r");
        string subString = strLine.substr(position, 2);
        if(subString != "::")
        {
            stringstream strStream;
            strStream << strLine;
            strStream >> speciesName;
        }
        else
        {
            stringstream strStream;
            strStream << strLine.substr(position + 2);

            strStream >> identiferDataTempora.identifer;

            position = strLine.find(identiferDataTempora.identifer);
            identiferDataTempora.valueAndComment = strLine.substr( position + identiferDataTempora.identifer.size() );


            identiferAddReplaceMultiMap.insert({speciesName, identiferDataTempora}); 

        } 

    } 

    file.close();

} 


void SpeciesMineralGas::readVmCriticalData()
{
    ifstream file(RELATIVEPATH + VMCRITICALDATA_F, ios::in); 
    string speciesName; 
    IdentifierData identiferDataTempora; 


    while(!file.eof())
    {
        string strLine;
        FileInput objectFileInput;
        objectFileInput.getLineExcludingComments_ButDoulbeColon(file, strLine); 
        if(strLine.empty())
            break;

        objectFileInput.excludeAqOfLine(strLine); 

        size_t position = strLine.find_first_not_of(" \t\r");
        string subString = strLine.substr(position, 2);
        if(subString != "::")
        {
            stringstream strStream;
            strStream << strLine;
            strStream >> speciesName;
        }
        else
        {
            stringstream strStream;
            strStream << strLine.substr(position + 2);

            strStream >> identiferDataTempora.identifer;

            position = strLine.find(identiferDataTempora.identifer);
            identiferDataTempora.valueAndComment = strLine.substr( position + identiferDataTempora.identifer.size() );


            identiferAddReplaceMultiMap.insert({speciesName, identiferDataTempora}); 

        } 

    } 

    file.close();
} 


void SpeciesMineralGas::readAbnormalGamma()
{
    ifstream file(RELATIVEPATH + ABNORMALGAMMA_F, ios::in); 

    while(!file.eof())
    {
        string strLine;
        FileInput objectFileInput;
        objectFileInput.getLineExcludingComments(file, strLine); 
        if(strLine.empty())
            break;

        objectFileInput.excludeAqOfLine(strLine); 

        stringstream strStream;
        strStream << strLine;

        string speciesName;
        double gamma;
        strStream >> speciesName; 
        strStream >> gamma; 

        abnormalGammaMap.insert({speciesName, gamma});

    } 

} 


double SpeciesMineralGas::obtainSpeciesValenceFromName(const string &str)
{
    string nameStr = str;

    Utility::trim(nameStr);

    size_t position;
    position = nameStr.find_first_of("+-");

    if(position == string::npos)
    {
        return 0;
    } 

    if(position + 1 == nameStr.size())

    {
        return stod( nameStr.substr(position, 1) + "1");
    }
    else
    {
        return stod( nameStr.substr(position) );
    }

} 


double SpeciesMineralGas::estimateSpeciesGammaFromName(const string &nameStr)
{
    double valence;
    valence = obtainSpeciesValenceFromName(nameStr);


    
    double gamma;
    if(valence < 0)
    {
        gamma = 4.0;
    }
    else if(valence == 0)
    {
        gamma = 3.0;
    }
    else if(valence > 0 && valence <= 4)
    {
        gamma = 3.5 + 0.5 * valence;
    }
    else
    {
        gamma = 6.0;
    } 


    
    for(auto iterm : abnormalGammaMap)
    {
        if(iterm.first == nameStr)
        {
            gamma = iterm.second;
            break;
        } 

    } 


    return gamma;

} 


void SpeciesMineralGas::calculateAllLogk()
{
    
    ofstream commandFile("command.txt",ios::out);
    FlagOfDatabaseFormat flag = ScreenInput::getTypeOfDatabaseFormat(); 

    switch(flag)
    {

    case PHREEQCDAT100:
    {
        commandFile << "n\n"
                    "dprons.dat\n"
                    "3\n"
                    "2\n"
                    "1\n"
                    "1"
                    << endl;

        double minT = ScreenInput::getMinTemperature(); 
        double maxT = ScreenInput::getMaxTemperature(); 
        commandFile << minT
                    << ", "
                    << maxT
                    << ", "
                    << (maxT - minT) / (SIZEOFTEMPERATURES - 1)
                    << endl;

        break;
    }

    case LIQVAP:
    {
        commandFile << "n\n"
                    "dprons.dat\n"
                    "3\n"
                    "2\n"
                    "1\n"
                    "1"
                    << endl;

        double minT = ScreenInput::getMinTemperature(); 
        double maxT = ScreenInput::getMaxTemperature(); 
        commandFile << minT
                    << ", "
                    << maxT
                    << ", "
                    << (maxT - minT) / (SIZEOFTEMPERATURES - 1)
                    << endl;

        break;
    }

    case LLNLDATPSAT:
    {
        commandFile << "n\n"
                    "dprons.dat\n"
                    "3\n"
                    "2\n"
                    "1\n"
                    "1"
                    << endl;

        double minT = ScreenInput::getMinTemperature(); 
        double maxT = ScreenInput::getMaxTemperature(); 
        commandFile << minT
                    << ", "
                    << maxT
                    << ", "
                    << (maxT - minT) / (SIZEOFTEMPERATURES - 1)
                    << endl;

        break;
    }

    case CONSTANTPREESURE:
    {
        commandFile << "n\n"
                    "dprons.dat\n"
                    "3\n"
                    "1\n"
                    "2\n"
                    "n\n"
                    "1\n"
                    "1"
                    << endl;

        double pressure = ScreenInput::getConstantPressure(); 
        commandFile << pressure << ", " << pressure << ", 0" << endl;

        double minT = ScreenInput::getMinTemperature(); 
        double maxT = ScreenInput::getMaxTemperature(); 
        commandFile << minT
                    << ", "
                    << maxT
                    << ", "
                    << (maxT - minT) / (SIZEOFTEMPERATURES - 1)
                    << endl;

        break;
    }

    default:
    {
        cout << "Should never reach here!" << endl;
        break;
    }

    } 

    commandFile << "n\n"
                << "1\n"
                << RELATIVEPATH + RXNFILE << "\n"
                << "logK.out\n"
                << "1"
                << endl;

    commandFile.close();


    
    remove("logK.out");


    if(::systemName == WINDOWS)
    {
        system("SUPCRTBL.exe < command.txt > NUL");
    }
    else
    {
        system("./SUPCRTBL < command.txt > NUL");
    } 

} 


void SpeciesMineralGas::obtainNameReactionFormulaFromRXN()
{
    ifstream file( RELATIVEPATH + RXNFILE , ios::in); 

    while(!file.eof())
    {
        string name; 
        string reactionFormula; 
        size_t typeSubstance = 2; 

        bool flag;
        FileInput objectFileInput;
        flag = objectFileInput.getNameReactionFormulaFromRXN(file, name, reactionFormula, typeSubstance);

        if (!flag)
            break;

        SecondarySpeciesChemicalData secondarySepciesDataStruct; 
        PhaseChemicalData phaseDataStruct;  
        if(typeSubstance == 1)
        {
            phaseDataStruct.reactionFormula = reactionFormula;

            mineralsMap.insert({name, phaseDataStruct});
        }
        else if(typeSubstance == 3)
        {
            phaseDataStruct.reactionFormula = reactionFormula;

            gasesMap.insert({name, phaseDataStruct});
        }
        else
        {
            if(secondaryMasterMap.find(name) != secondaryMasterMap.end())
            {
                secondaryMasterMap[name].reactionFormula = reactionFormula;
            }
            else
            {
                secondarySepciesDataStruct.chemicalFormula = name;
                secondarySepciesDataStruct.reactionFormula = reactionFormula;

                secondarySpeciesMap.insert({name, secondarySepciesDataStruct});

            } 

        } 

    } 


    file.close();

} 


void SpeciesMineralGas::obtainPhaseChemicalFormulaFromSUPCRT()
{
    ifstream file( RELATIVEPATH + SPRONSFILE, ios::in); 
    string strLine;


    
    while(!file.eof())
    {
        getline(file,strLine);

        if(strLine.find("minerals that do not undergo phase transitions") != string::npos)
            break;

    }

    getline(file,strLine);
    getline(file,strLine);

    size_t positionStar = strLine.find("*******************************************************");

    while(positionStar == string::npos)
    {

        stringstream strStream;
        strStream.str(strLine);

        string name; 
        string chemicalFormula; 

        
        
        name = strLine.substr(0, 21);
        Utility::trim(name);
        chemicalFormula = strLine.substr(21);
        Utility::trim(chemicalFormula);

        if (mineralsMap.find(name) != mineralsMap.end())
        {
            mineralsMap[name].chemicalFormula = chemicalFormula;

            string temporaryStr = mineralsMap[name].reactionFormula;
            temporaryStr.replace( temporaryStr.find(name), name.size(), chemicalFormula );
            mineralsMap[name].reactionFormula = temporaryStr;

        } 

        for (size_t i = 0; i < 7; ++i)
            getline(file,strLine);

        positionStar = strLine.find("*******************************************************");

    } 


    
    for (size_t i = 0; i < 3; ++i)
        getline(file, strLine);

    positionStar = strLine.find("*******************************************************");
    while(positionStar == string::npos)
    {
        stringstream strStream;
        strStream.str(strLine);

        string name; 
        string chemicalFormula; 

        
        
        name = strLine.substr(0, 21);
        Utility::trim(name);
        chemicalFormula = strLine.substr(21);
        Utility::trim(chemicalFormula);

        if (mineralsMap.find(name) != mineralsMap.end())
        {
            mineralsMap[name].chemicalFormula = chemicalFormula;

            string temporaryStr = mineralsMap[name].reactionFormula;
            temporaryStr.replace( temporaryStr.find(name), name.size(), chemicalFormula );
            mineralsMap[name].reactionFormula = temporaryStr;

        } 

        for (size_t i = 0; i < 8; ++i)
            getline(file,strLine);

        positionStar = strLine.find("*******************************************************");
    } 


    
    for (size_t i = 0; i < 3; ++i)
        getline(file, strLine);

    positionStar = strLine.find("*******************************************************");
    while(positionStar == string::npos)
    {
        stringstream strStream;
        strStream.str(strLine);

        string name; 
        string chemicalFormula; 

        
        
        name = strLine.substr(0, 21);
        Utility::trim(name);
        chemicalFormula = strLine.substr(21);
        Utility::trim(chemicalFormula);

        if (mineralsMap.find(name) != mineralsMap.end())
        {
            mineralsMap[name].chemicalFormula = chemicalFormula;

            string temporaryStr = mineralsMap[name].reactionFormula;
            temporaryStr.replace( temporaryStr.find(name), name.size(), chemicalFormula );
            mineralsMap[name].reactionFormula = temporaryStr;

        } 

        for (size_t i = 0; i < 8; ++i)
            getline(file,strLine);

        positionStar = strLine.find("*******************************************************");

    } 


    
    for (size_t i = 0; i < 6; ++i)
        getline(file, strLine);

    positionStar = strLine.find("*******************************************************");
    while(positionStar == string::npos)
    {
        stringstream strStream;
        strStream.str(strLine);

        string name; 
        string chemicalFormula; 

        
        
        name = strLine.substr(0, 21);
        Utility::trim(name);
        chemicalFormula = strLine.substr(21);
        Utility::trim(chemicalFormula);

        if (gasesMap.find(name) != gasesMap.end())
        {
            gasesMap[name].chemicalFormula = chemicalFormula;

            string temporaryStr = gasesMap[name].reactionFormula;
            temporaryStr.replace( temporaryStr.find(name), name.size(), chemicalFormula );
            gasesMap[name].reactionFormula = temporaryStr;

        } 

        for (size_t i = 0; i < 6; ++i)
            getline(file,strLine);

        positionStar = strLine.find("*******************************************************");

    } 


    file.close();

} 


void SpeciesMineralGas::calculateStoreLogKAnalytic()
{
    ifstream file("logK.out", ios::in); 

    string name; 
    double logkArray[SIZEOFTEMPERATURES] = {0.0}; 
    double densityArray[SIZEOFTEMPERATURES] = {0.0}; 
    double pressureArray[SIZEOFTEMPERATURES] = {0.0}; 

    while(!file.eof())
    {

        
        bool flag;

        FileInput objectFileInput;
        flag = objectFileInput.getNameLogkDensityPressure(file, name, logkArray, densityArray, pressureArray, SIZEOFTEMPERATURES);

        if (!flag)
            break;

        objectFileInput.excludeAqOfLine(name); 
        Utility::trim(name); 

        double minT = ScreenInput::getMinTemperature(); 
        double maxT = ScreenInput::getMaxTemperature(); 
        double temperatureArray[SIZEOFTEMPERATURES]; 
        Utility::linspace(minT,maxT,SIZEOFTEMPERATURES,temperatureArray); 

        double coefficientsArray[6]; 
        double maxRelativeError; 
        double correspondingLogk; 
        double corresponAbsoluteErrorLogk; 
        double corresponTemporature; 

        Utility::estimateLogkAnalyticFormula(temperatureArray, logkArray,
                                             SIZEOFTEMPERATURES, coefficientsArray,
                                             maxRelativeError, corresponAbsoluteErrorLogk,
                                             correspondingLogk, corresponTemporature);


        
        if (maxRelativeError > 5 && fabs(corresponAbsoluteErrorLogk) > 0.05)
            logkErrorWarn = true; 


        
        if(secondaryMasterMap.find(name) != secondaryMasterMap.end())
        {
            for(size_t i = 0; i < 6; ++i)
                secondaryMasterMap[name].analytic[i] = coefficientsArray[i];

            secondaryMasterMap[name].maxRelativeErrorLogk = maxRelativeError;
            secondaryMasterMap[name].corresponAbsoluteErrorLogk = corresponAbsoluteErrorLogk;
            secondaryMasterMap[name].correspondingLogk = correspondingLogk;
            secondaryMasterMap[name].corresponTemporature = corresponTemporature;
        }
        else if (secondarySpeciesMap.find(name) != secondarySpeciesMap.end())
        {
            for(size_t i = 0; i < 6; ++i)
                secondarySpeciesMap[name].analytic[i] = coefficientsArray[i];

            secondarySpeciesMap[name].maxRelativeErrorLogk = maxRelativeError;
            secondarySpeciesMap[name].corresponAbsoluteErrorLogk = corresponAbsoluteErrorLogk;
            secondarySpeciesMap[name].correspondingLogk = correspondingLogk;
            secondarySpeciesMap[name].corresponTemporature = corresponTemporature;
        }
        else if(mineralsMap.find(name) != mineralsMap.end())
        {
            for(size_t i = 0; i < 6; ++i)
                mineralsMap[name].analytic[i] = coefficientsArray[i];

            mineralsMap[name].maxRelativeErrorLogk = maxRelativeError;
            mineralsMap[name].corresponAbsoluteErrorLogk = corresponAbsoluteErrorLogk;
            mineralsMap[name].correspondingLogk = correspondingLogk;
            mineralsMap[name].corresponTemporature = corresponTemporature;
        }
        else if(gasesMap.find(name) != gasesMap.end())
        {
            for(size_t i = 0; i < 6; ++i)
                gasesMap[name].analytic[i] = coefficientsArray[i];

            gasesMap[name].maxRelativeErrorLogk = maxRelativeError;
            gasesMap[name].corresponAbsoluteErrorLogk = corresponAbsoluteErrorLogk;
            gasesMap[name].correspondingLogk = correspondingLogk;
            gasesMap[name].corresponTemporature = corresponTemporature;
        }
        else
        {

        } 

    } 

    file.close();

} 


void SpeciesMineralGas::calculateStoreGamma()
{

    for(auto iter = primaryMasterMap.begin(); iter != primaryMasterMap.end(); ++iter)
    {
        string nameStr = iter->first;
        primaryMasterMap[iter->first].llnlGamma = estimateSpeciesGammaFromName(nameStr);

    } 


    for(auto iter = secondaryMasterMap.begin(); iter != secondaryMasterMap.end(); ++iter)
    {
        string nameStr = iter->first;
        secondaryMasterMap[iter->first].llnlGamma = estimateSpeciesGammaFromName(nameStr);

    } 


    for(auto iter = secondarySpeciesMap.begin(); iter != secondarySpeciesMap.end(); ++iter)
    {
        string nameStr = iter->first;
        secondarySpeciesMap[iter->first].llnlGamma = estimateSpeciesGammaFromName(nameStr);

    } 

} 


void SpeciesMineralGas::addExtraIdentifier()
{
    for(auto iterm : identiferAddReplaceMultiMap)
    {
        string name = iterm.first; 
        string identifer = iterm.second.identifer; 
        string valueAndComment = iterm.second.valueAndComment; 

        
        if(primaryMasterMap.find(name) != primaryMasterMap.end())
        {
            primaryMasterMap[name].extraIdentiferMap.insert({identifer, valueAndComment});
        }
        else if(secondaryMasterMap.find(name) != secondaryMasterMap.end())
        {
            secondaryMasterMap[name].extraIdentiferMap.insert({identifer, valueAndComment});
        }
        else if(secondarySpeciesMap.find(name) != secondarySpeciesMap.end())
        {
            secondarySpeciesMap[name].extraIdentiferMap.insert({identifer, valueAndComment});
        }
        else if(mineralsMap.find(name) != mineralsMap.end())
        {
            mineralsMap[name].extraIdentiferMap.insert({identifer, valueAndComment});
        }
        else if(gasesMap.find(name) != gasesMap.end())
        {
            gasesMap[name].extraIdentiferMap.insert({identifer, valueAndComment});
        }
        else
        {

        } 

    }

} 


void SpeciesMineralGas::treatO2()
{
    
    secondaryMasterMap["O2"].reactionFormula = O2REACTIONFORMULA;


    
    double minT = ScreenInput::getMinTemperature(); 
    double maxT = ScreenInput::getMaxTemperature(); 
    double temperatureArray[SIZEOFTEMPERATURES]; 
    Utility::linspace(minT,maxT,SIZEOFTEMPERATURES,temperatureArray); 

    double coefficientsArray[6]; 
    double maxRelativeError; 
    double correspondingLogk; 
    double corresponAbsoluteErrorLogk; 
    double corresponTemporature; 

    Utility::estimateLogkAnalyticFormula(temperatureArray, AqueousModel::getO2LogkArrayPtr(),
                                         SIZEOFTEMPERATURES, coefficientsArray,
                                         maxRelativeError, corresponAbsoluteErrorLogk,
                                         correspondingLogk, corresponTemporature);

    for (size_t i = 0; i < 6; ++i)
        secondaryMasterMap["O2"].analytic[i] = coefficientsArray[i];

    secondaryMasterMap["O2"].maxRelativeErrorLogk = maxRelativeError;
    secondaryMasterMap["O2"].corresponAbsoluteErrorLogk = corresponAbsoluteErrorLogk;
    secondaryMasterMap["O2"].correspondingLogk = correspondingLogk;
    secondaryMasterMap["O2"].corresponTemporature = corresponTemporature;

} 


void SpeciesMineralGas::outputSOLUTION_SPECIESPHASES_LIQVAP()
{
    ofstream &file = FileOutput::getPhreeqcFileReference(); 

    file << "SOLUTION_SPECIES" << endl;


    
    file << fixed << setprecision(2);

    for (auto iterm1 : primaryMasterMap)
    {
        file << iterm1.second.reactionFormula << endl;


        if (iterm1.second.extraIdentiferMap.find("-log_k") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-log_k" << iterm1.second.extraIdentiferMap["-log_k"] << endl;

            iterm1.second.extraIdentiferMap.erase("-log_k");
        }
        else
        {
            file << "\t" << "-log_k" << "\t" << iterm1.second.logK << endl;
        } 


        if (iterm1.second.extraIdentiferMap.find("-delta_h") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-delta_h" << iterm1.second.extraIdentiferMap["-delta_h"] << endl;

            iterm1.second.extraIdentiferMap.erase("-delta_h");
        }
        else
        {
            file << "\t" << "-delta_h" << "\t" << iterm1.second.deltaH << endl;
        } 


        for (auto iterm2 : iterm1.second.extraIdentiferMap)
            file << "\t" << iterm2.first << iterm2.second << endl;

    } 
    file << endl;


    
    file << "#**************** secondary master species *********************" << endl;
    for (auto iterm1 : secondaryMasterMap)
    {
        file << iterm1.second.reactionFormula << endl;


        if (iterm1.second.extraIdentiferMap.find("-analytic") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-analytic" << iterm1.second.extraIdentiferMap["-analytic"] << endl;

            iterm1.second.extraIdentiferMap.erase("-analytic");
        }
        else
        {
            file << "\t" << "-analytic";
            file << setprecision(6) << scientific;

            for (size_t i = 0; i < 6; ++i)
                file << "\t" << iterm1.second.analytic[i];

            file << endl;
            file << fixed;

        } 


        for (auto iterm2 : iterm1.second.extraIdentiferMap)
            file << "\t" << iterm2.first << iterm2.second << endl;

    } 
    file << endl;


    
    file << "#**************** secondary species *********************" << endl;
    for (auto iterm1 : secondarySpeciesMap)
    {
        file << iterm1.second.reactionFormula << endl;


        if (iterm1.second.extraIdentiferMap.find("-analytic") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-analytic" << iterm1.second.extraIdentiferMap["-analytic"] << endl;

            iterm1.second.extraIdentiferMap.erase("-analytic");
        }
        else
        {
            file << "\t" << "-analytic";
            file << setprecision(6) << scientific;

            for (size_t i = 0; i < 6; ++i)
                file << "\t" << iterm1.second.analytic[i];

            file << endl;
            file << fixed;

        } 


        for (auto iterm2 : iterm1.second.extraIdentiferMap)
            file << "\t" << iterm2.first << iterm2.second << endl;

    } 
    file << "\n" << endl;


    
    file << "PHASES" << endl;

    for (auto iterm1 : mineralsMap)
    {
        file << iterm1.first << endl;
        file << iterm1.second.reactionFormula << endl;

        if (iterm1.second.extraIdentiferMap.find("-analytic") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-analytic" << iterm1.second.extraIdentiferMap["-analytic"] << endl;

            iterm1.second.extraIdentiferMap.erase("-analytic");
        }
        else
        {
            file << "\t" << "-analytic";
            file << setprecision(6) << scientific;

            for (size_t i = 0; i < 6; ++i)
                file << "\t" << iterm1.second.analytic[i];

            file << endl;
            file << fixed;

        } 


        for (auto iterm2 : iterm1.second.extraIdentiferMap)
            file << "\t" << iterm2.first << iterm2.second << endl;

    } 


    
    for (auto iterm1 : gasesMap)
    {
        file << iterm1.first << endl;
        file << iterm1.second.reactionFormula << endl;

        if (iterm1.second.extraIdentiferMap.find("-analytic") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-analytic" << iterm1.second.extraIdentiferMap["-analytic"] << endl;

            iterm1.second.extraIdentiferMap.erase("-analytic");
        }
        else
        {
            file << "\t" << "-analytic";
            file << setprecision(6) << scientific;

            for (size_t i = 0; i < 6; ++i)
                file << "\t" << iterm1.second.analytic[i];

            file << endl;
            file << fixed;

        } 


        for (auto iterm2 : iterm1.second.extraIdentiferMap)
            file << "\t" << iterm2.first << iterm2.second << endl;

    } 


    
    remove("command.txt");
    remove("dprons.dat");
    remove("zero.dat");
    remove("warning.txt");

} 


void SpeciesMineralGas::outputSOLUTION_SPECIESPHASES_CONSTANTPREESURE()
{
    ofstream &file = FileOutput::getPhreeqcFileReference(); 

    file << "SOLUTION_SPECIES" << endl;


    
    file << fixed << setprecision(2);

    for (auto iterm1 : primaryMasterMap)
    {
        file << iterm1.second.reactionFormula << endl;


        
        if (iterm1.second.extraIdentiferMap.find("-llnl_gamma") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-llnl_gamma" << iterm1.second.extraIdentiferMap["-llnl_gamma"] << endl;

            iterm1.second.extraIdentiferMap.erase("-llnl_gamma");
        }
        else if (iterm1.second.extraIdentiferMap.find("-CO2_llnl_gamma") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-CO2_llnl_gamma" << iterm1.second.extraIdentiferMap["-CO2_llnl_gamma"] << endl;

            iterm1.second.extraIdentiferMap.erase("-CO2_llnl_gamma");
        }
        else
        {
            file << "\t" << "-llnl_gamma" << "\t" << iterm1.second.llnlGamma << endl;
        } 


        if (iterm1.second.extraIdentiferMap.find("-log_k") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-log_k" << iterm1.second.extraIdentiferMap["-log_k"] << endl;

            iterm1.second.extraIdentiferMap.erase("-log_k");
        }
        else
        {
            file << "\t" << "-log_k" << "\t" << iterm1.second.logK << endl;
        } 


        if (iterm1.second.extraIdentiferMap.find("-delta_h") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-delta_h" << iterm1.second.extraIdentiferMap["-delta_h"] << endl;

            iterm1.second.extraIdentiferMap.erase("-delta_h");
        }
        else
        {
            file << "\t" << "-delta_h" << "\t" << iterm1.second.deltaH << endl;
        } 


        for (auto iterm2 : iterm1.second.extraIdentiferMap)
            file << "\t" << iterm2.first << iterm2.second << endl;

    } 
    file << endl;


    
    file << "#**************** secondary master species *********************" << endl;
    for (auto iterm1 : secondaryMasterMap)
    {
        file << iterm1.second.reactionFormula << endl;


        
        if (iterm1.second.extraIdentiferMap.find("-llnl_gamma") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-llnl_gamma" << iterm1.second.extraIdentiferMap["-llnl_gamma"] << endl;

            iterm1.second.extraIdentiferMap.erase("-llnl_gamma");
        }
        else if (iterm1.second.extraIdentiferMap.find("-CO2_llnl_gamma") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-CO2_llnl_gamma" << iterm1.second.extraIdentiferMap["-CO2_llnl_gamma"] << endl;

            iterm1.second.extraIdentiferMap.erase("-CO2_llnl_gamma");
        }
        else
        {
            file << "\t" << "-llnl_gamma" << "\t" << setprecision(2) << iterm1.second.llnlGamma << endl;
        } 


        if (iterm1.second.extraIdentiferMap.find("-analytic") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-analytic" << iterm1.second.extraIdentiferMap["-analytic"] << endl;

            iterm1.second.extraIdentiferMap.erase("-analytic");
        }
        else
        {
            file << "\t" << "-analytic";
            file << setprecision(6) << scientific;

            for (size_t i = 0; i < 6; ++i)
                file << "\t" << iterm1.second.analytic[i];

            file << endl;
            file << fixed;

        } 


        for (auto iterm2 : iterm1.second.extraIdentiferMap)
            file << "\t" << iterm2.first << iterm2.second << endl;

    } 
    file << endl;


    
    file << "#**************** secondary species *********************" << endl;
    for (auto iterm1 : secondarySpeciesMap)
    {
        file << iterm1.second.reactionFormula << endl;


        
        if (iterm1.second.extraIdentiferMap.find("-llnl_gamma") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-llnl_gamma" << iterm1.second.extraIdentiferMap["-llnl_gamma"] << endl;

            iterm1.second.extraIdentiferMap.erase("-llnl_gamma");
        }
        else if (iterm1.second.extraIdentiferMap.find("-CO2_llnl_gamma") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-CO2_llnl_gamma" << iterm1.second.extraIdentiferMap["-CO2_llnl_gamma"] << endl;

            iterm1.second.extraIdentiferMap.erase("-CO2_llnl_gamma");
        }
        else
        {
            file << "\t" << "-llnl_gamma" << "\t" << setprecision(2) << iterm1.second.llnlGamma << endl;
        } 


        if (iterm1.second.extraIdentiferMap.find("-analytic") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-analytic" << iterm1.second.extraIdentiferMap["-analytic"] << endl;

            iterm1.second.extraIdentiferMap.erase("-analytic");
        }
        else
        {
            file << "\t" << "-analytic";
            file << setprecision(6) << scientific;

            for (size_t i = 0; i < 6; ++i)
                file << "\t" << iterm1.second.analytic[i];

            file << endl;
            file << fixed;

        } 


        for (auto iterm2 : iterm1.second.extraIdentiferMap)
            file << "\t" << iterm2.first << iterm2.second << endl;

    } 
    file << "\n" << endl;


    
    file << "PHASES" << endl;

    for (auto iterm1 : mineralsMap)
    {
        file << iterm1.first << endl;
        file << iterm1.second.reactionFormula << endl;

        if (iterm1.second.extraIdentiferMap.find("-analytic") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-analytic" << iterm1.second.extraIdentiferMap["-analytic"] << endl;

            iterm1.second.extraIdentiferMap.erase("-analytic");
        }
        else
        {
            file << "\t" << "-analytic";
            file << setprecision(6) << scientific;

            for (size_t i = 0; i < 6; ++i)
                file << "\t" << iterm1.second.analytic[i];

            file << endl;
            file << fixed;

        } 


        for (auto iterm2 : iterm1.second.extraIdentiferMap)
            file << "\t" << iterm2.first << iterm2.second << endl;

    } 


    
    for (auto iterm1 : gasesMap)
    {
        file << iterm1.first << endl;
        file << iterm1.second.reactionFormula << endl;

        if (iterm1.second.extraIdentiferMap.find("-analytic") != iterm1.second.extraIdentiferMap.end())
        {
            file << "\t" << "-analytic" << iterm1.second.extraIdentiferMap["-analytic"] << endl;

            iterm1.second.extraIdentiferMap.erase("-analytic");
        }
        else
        {
            file << "\t" << "-analytic";
            file << setprecision(6) << scientific;

            for (size_t i = 0; i < 6; ++i)
                file << "\t" << iterm1.second.analytic[i];

            file << endl;
            file << fixed;

        } 


        for (auto iterm2 : iterm1.second.extraIdentiferMap)
            file << "\t" << iterm2.first << iterm2.second << endl;

    } 


    
    remove("command.txt");
    remove("dprons.dat");
    remove("zero.dat");
    remove("warning.txt");

} 


void SpeciesMineralGas::outputWarnInformation()
{
    if(!logkErrorWarn)
        return;

    cout << endl;
    cout << "Warning! Maximum relative error of logK analytic formula fitting is > 5%.\n"
         << "Before using database file, please check warning information in \"warning.txt\"." << endl;
    cout << "\n" << endl;


    ofstream file("warning.txt", ios::out); 

    file << "# Maximum relative error of logK analytic formula fitting is > 5%, and corresponding absolute logK error is > 0.05.\n"
         << "# Detailed information is listed below. Please carefully check these information before using the database.\n\n"
         << endl;

    file << setw(25) << "species or phase,"
         << setw(25) << "maximum relative error,"
         << setw(30) << "corresponding absolute error"
         << "\n" << endl;

    file << fixed << setprecision(2);

    for(auto iterm : secondaryMasterMap)
    {
        double maxRelativeErrorLogk = secondaryMasterMap[iterm.first].maxRelativeErrorLogk;
        double absoluteErrorLogk = secondaryMasterMap[iterm.first].corresponAbsoluteErrorLogk;

        if(maxRelativeErrorLogk > 5 && absoluteErrorLogk > 0.05)
        {
            file << setw(25) << iterm.first;
            file << setw(25) << maxRelativeErrorLogk;
            file << setw(30) << absoluteErrorLogk;
            file << endl;
        } 

    } 


    for(auto iterm : secondarySpeciesMap)
    {
        double maxRelativeErrorLogk = secondarySpeciesMap[iterm.first].maxRelativeErrorLogk;
        double absoluteErrorLogk = secondarySpeciesMap[iterm.first].corresponAbsoluteErrorLogk;

        if(maxRelativeErrorLogk > 5 && absoluteErrorLogk > 0.05)
        {
            file << setw(25) << iterm.first;
            file << setw(25) << maxRelativeErrorLogk;
            file << setw(30) << absoluteErrorLogk;
            file << endl;
        } 

    } 


    for(auto iterm : mineralsMap)
    {
        double maxRelativeErrorLogk = mineralsMap[iterm.first].maxRelativeErrorLogk;
        double absoluteErrorLogk = mineralsMap[iterm.first].corresponAbsoluteErrorLogk;

        if(maxRelativeErrorLogk > 5 && absoluteErrorLogk > 0.05)
        {
            file << setw(25) << iterm.first;
            file << setw(25) << maxRelativeErrorLogk;
            file << setw(30) << absoluteErrorLogk;
            file << endl;
        } 

    } 


    for(auto iterm : gasesMap)
    {
        double maxRelativeErrorLogk = gasesMap[iterm.first].maxRelativeErrorLogk;
        double absoluteErrorLogk = gasesMap[iterm.first].corresponAbsoluteErrorLogk;

        if(maxRelativeErrorLogk > 5 && absoluteErrorLogk > 0.05)
        {
            file << setw(25) << iterm.first;
            file << setw(25) << maxRelativeErrorLogk;
            file << setw(30) << absoluteErrorLogk;
            file << endl;
        } 

    } 

    file.close();

}



void SpeciesMineralGas::runInstructs()
{
    FlagOfDatabaseFormat flag = ScreenInput::getTypeOfDatabaseFormat(); 

    switch(flag)
    {

    case PHREEQCDAT100:
    case LIQVAP:
    {
        getOutputSolutionMasterSpecies();
        readExtraIdentifiers();

        
        readVmCriticalData();

        calculateAllLogk();
        obtainNameReactionFormulaFromRXN();
        obtainPhaseChemicalFormulaFromSUPCRT();
        calculateStoreLogKAnalytic();
        addExtraIdentifier();
        treatO2();

        outputSOLUTION_SPECIESPHASES_LIQVAP();

        FileOutput objectFileOutput;
        objectFileOutput.outputKineticsScript();

        outputWarnInformation();
        break;
    }

    case LLNLDATPSAT:
    {
        getOutputSolutionMasterSpecies();
        readExtraIdentifiers();
        readAbnormalGamma();
        calculateAllLogk();
        obtainNameReactionFormulaFromRXN();
        obtainPhaseChemicalFormulaFromSUPCRT();
        calculateStoreLogKAnalytic();
        calculateStoreGamma();
        addExtraIdentifier();
        treatO2();

        outputSOLUTION_SPECIESPHASES_CONSTANTPREESURE();

        FileOutput objectFileOutput;
        objectFileOutput.outputKineticsScript();

        outputWarnInformation();
        break;
    }

    case CONSTANTPREESURE:
    {
        getOutputSolutionMasterSpecies();
        readExtraIdentifiers();
        readAbnormalGamma();
        calculateAllLogk();
        obtainNameReactionFormulaFromRXN();
        obtainPhaseChemicalFormulaFromSUPCRT();
        calculateStoreLogKAnalytic();
        calculateStoreGamma();
        addExtraIdentifier();
        treatO2();

        outputSOLUTION_SPECIESPHASES_CONSTANTPREESURE();
        outputWarnInformation();
        break;
    }

    default:
    {
        cout << "Should never reach here!" << endl;
        break;
    }

    } 

} 
