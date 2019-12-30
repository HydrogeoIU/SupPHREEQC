#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <cstddef>
#include <map>



enum FlagOfDatabaseFormat {LIQVAP = 1, CONSTANTPREESURE, PHREEQCDAT100, LLNLDATPSAT}; 

enum GoodBad {BAD, GOOD}; 

enum SystemInf {WINDOWS, LINUX}; 



struct PrimaryMasterChemicalData 
{
    PrimaryMasterChemicalData(); 

    std::string chemicalFormula; 
    std::string reactionFormula; 
    double llnlGamma; 
    double logK; 
    double deltaH; 
    std::string element; 
    double valence; 

    std::map<std::string, std::string> extraIdentiferMap; 
};

struct SecondaryMasterChemicalData 
{
    SecondaryMasterChemicalData();

    std::string chemicalFormula; 
    std::string reactionFormula; 
    double llnlGamma; 
    double analytic[6]; 
    double maxRelativeErrorLogk; 
    double corresponAbsoluteErrorLogk; 
    double correspondingLogk; 
    double corresponTemporature; 
    std::string element; 
    double valence; 

    std::map<std::string, std::string> extraIdentiferMap; 
};

struct SecondarySpeciesChemicalData 
{
    SecondarySpeciesChemicalData();

    std::string chemicalFormula; 
    std::string reactionFormula; 
    double llnlGamma; 
    double analytic[6]; 
    double maxRelativeErrorLogk; 
    double corresponAbsoluteErrorLogk; 
    double correspondingLogk; 
    double corresponTemporature; 

    std::map<std::string, std::string> extraIdentiferMap; 
};

struct PhaseChemicalData 
{
    PhaseChemicalData();

    std::string chemicalFormula; 
    std::string reactionFormula; 
    double analytic[6]; 
    double maxRelativeErrorLogk; 
    double corresponAbsoluteErrorLogk; 
    double correspondingLogk; 
    double corresponTemporature; 

    std::map<std::string, std::string> extraIdentiferMap; 
};

struct IdentifierData 
{
    IdentifierData();

    std::string identifer; 
    std::string valueAndComment; 
};



const std::string NAMEOFPHREEQCDATABASE = "PHREEQC-temporary.DAT";

const std::string SPRONSFILE = "spronsbl.dat";

const std::string BDOT_TP = "bdot vs TP.txt";

const std::string SOLUTION_MASTER_SPECIES_F = "SOLUTION_MASTER_SPECIES.txt";

const std::string ABNORMALGAMMA_F = "abnormal gamma.txt";

const std::string EXTRAIDENTIFIER_F = "addreplace.txt";

const std::string EXTRAIDENTIFIER_PHREEQCDAT_F = "addreplace-PHREEQCDAT.txt";

const std::string RXNFILE = "rxn.dat";

const std::string VMCRITICALDATA_F = "vmcriticalpoint.txt";

const std::string KINETICSSCRIPT_F = "kineticsScript.txt";

#if defined(_WIN32)
    const std::string RELATIVEPATH = ".\\data\\"; 

#elif defined(__linux__)
    const std::string RELATIVEPATH = "./data/"; 

#else
#error "OS not supported!"
#endif



const size_t SIZEOFTEMPERATURES = 11; 

const std::string O2REACTIONFORMULA = "2 H2O = O2 + 4 H+ + 4 e-";

const std::string CO2COEFS = {"-co2_coefs        # range: 0-400 oC, 0-500 bar\n"
                              "        -1.0312              0.0012806\n"
                              "          255.9                 0.4445\n"
                              "      -0.001606\n"
                             };

const std::string NAMEALKALINITYELEMENT = "Alkalinity";

const std::string SPECIESALKALINITYELEMENT = "HCO3-";



extern std::map<std::string, PrimaryMasterChemicalData, std::less<std::string>> primaryMasterMap; 

extern std::map<std::string, SecondaryMasterChemicalData, std::less<std::string>> secondaryMasterMap; 

extern std::map<std::string, SecondarySpeciesChemicalData, std::less<std::string>> secondarySpeciesMap; 

extern std::map<std::string, PhaseChemicalData, std::less<std::string>> mineralsMap; 

extern std::map<std::string, PhaseChemicalData, std::less<std::string>> gasesMap; 

extern std::map<std::string, double, std::less<std::string>> abnormalGammaMap; 

extern std::multimap<std::string, IdentifierData, std::less<std::string>> identiferAddReplaceMultiMap; 

extern bool logkErrorWarn; 

extern SystemInf systemName; 


#endif 
