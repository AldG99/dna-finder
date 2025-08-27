#ifndef DNASEQUENCE_H
#define DNASEQUENCE_H

#include <string>
#include <map>
#include <vector>

class DNASequence {
private:
    std::string sequence;
    bool valid;
    std::map<char, int> nucleotideCount;
    
    void validateSequence();
    void countNucleotides();

public:
    DNASequence();
    DNASequence(const std::string& seq);
    
    void setSequence(const std::string& seq);
    std::string getSequence() const;
    bool isValid() const;
    
    std::string getComplement() const;
    std::string getReverseComplement() const;
    
    double getGCContent() const;
    int getNucleotideCount(char nucleotide) const;
    std::map<char, int> getAllCounts() const;
    double getMolecularWeight() const;
    
    int getLength() const;
    bool isEmpty() const;
    
    static bool isValidNucleotide(char nucleotide);
    static char getComplementNucleotide(char nucleotide);
};

#endif