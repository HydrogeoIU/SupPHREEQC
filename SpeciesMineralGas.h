#ifndef SPECIESMINERALGAS_H
#define SPECIESMINERALGAS_H


class SpeciesMineralGas
{
public:

    
    void getOutputSolutionMasterSpecies();


    
    void readExtraIdentifiers();


    
    void readVmCriticalData();


    
    void readAbnormalGamma();


    
    double obtainSpeciesValenceFromName(const std::string &str);


    
    double estimateSpeciesGammaFromName(const std::string &nameStr);


    
    void calculateAllLogk();


    
    void obtainNameReactionFormulaFromRXN();


    
    void obtainPhaseChemicalFormulaFromSUPCRT();


    
    void calculateStoreLogKAnalytic();


    
    void calculateStoreGamma();


    
    void addExtraIdentifier();


    
    void treatO2();


    
    void outputSOLUTION_SPECIESPHASES_LIQVAP();


    
    void outputSOLUTION_SPECIESPHASES_CONSTANTPREESURE();


    
    void outputWarnInformation();



    
    void runInstructs();

}; 


#endif 
