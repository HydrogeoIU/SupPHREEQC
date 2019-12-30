#include "Global.h"


using namespace std;



PrimaryMasterChemicalData::PrimaryMasterChemicalData()
    : chemicalFormula(""),
      reactionFormula(""),
      llnlGamma(0.0),
      logK(0.0),
      deltaH(0.0),
      element(""),
      valence(0.0)
{

} 

SecondaryMasterChemicalData::SecondaryMasterChemicalData()
    : chemicalFormula(""),
      reactionFormula(""),
      llnlGamma(0.0),
      analytic{0.0},
      maxRelativeErrorLogk(0.0),
      corresponAbsoluteErrorLogk(0.0),
      correspondingLogk(0.0),
      corresponTemporature(0.0),
      element(""),
      valence(0.0)
{

} 

SecondarySpeciesChemicalData::SecondarySpeciesChemicalData()
    : chemicalFormula(""),
      reactionFormula(""),
      llnlGamma(0.0),
      analytic{0.0},
      maxRelativeErrorLogk(0.0),
      corresponAbsoluteErrorLogk(0.0),
      correspondingLogk(0.0),
      corresponTemporature(0.0)
{

} 

PhaseChemicalData::PhaseChemicalData()
    : chemicalFormula(""),
      reactionFormula(""),
      analytic{0.0},
      maxRelativeErrorLogk(0.0),
      corresponAbsoluteErrorLogk(0.0),
      correspondingLogk(0.0),
      corresponTemporature(0.0)
{

} 

IdentifierData::IdentifierData()
    : identifer(""),
      valueAndComment("")
{

}



map<string, PrimaryMasterChemicalData, less<string>> primaryMasterMap;

map<string, SecondaryMasterChemicalData, less<string>> secondaryMasterMap;

map<string, SecondarySpeciesChemicalData, less<string>> secondarySpeciesMap;

map<string, PhaseChemicalData, less<string>> mineralsMap;

map<string, PhaseChemicalData, less<string>> gasesMap;

map<string, double, less<string>> abnormalGammaMap;

multimap<string, IdentifierData, less<string>> identiferAddReplaceMultiMap;

bool logkErrorWarn = false;


#if defined(_WIN32)
    SystemInf systemName = WINDOWS;

#elif defined(__linux__)
    SystemInf systemName = LINUX;

#else
#error "OS not supported!"

#endif
