#ifndef GENETICCODE_H
#define GENETICCODE_H

#include <string>
#include <map>
#include <vector>

class GeneticCode {
private:
    static std::map<std::string, char> codonTable;
    static std::map<char, std::string> codonToAminoAcid;
    static void initializeCodonTable();
    static bool tableInitialized;

public:
    static char translateCodon(const std::string& codon);
    static std::string translateSequence(const std::string& dnaSequence);
    static std::string translateSequenceVerbose(const std::string& dnaSequence);
    static bool isStartCodon(const std::string& codon);
    static bool isStopCodon(const std::string& codon);
    static std::string getAminoAcidName(char aminoAcid);
    static std::vector<int> findStartCodons(const std::string& dnaSequence);
    static std::vector<int> findStopCodons(const std::string& dnaSequence);
};

#endif