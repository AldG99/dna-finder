#ifndef SEQUENCEANALYZER_H
#define SEQUENCEANALYZER_H

#include <string>
#include <vector>
#include "DNASequence.h"

struct ORF {
    int start;
    int end;
    int frame;
    std::string sequence;
    std::string protein;
    int length;
    
    ORF(int s, int e, int f, const std::string& seq, const std::string& prot) 
        : start(s), end(e), frame(f), sequence(seq), protein(prot), length(prot.length()) {}
};

class SequenceAnalyzer {
public:
    static std::vector<ORF> findORFs(const std::string& dnaSequence, int minLength = 10);
    static std::vector<ORF> findORFsAllFrames(const std::string& dnaSequence, int minLength = 10);
    static double calculateMeltingTemperature(const std::string& dnaSequence);
    static double calculateComplexity(const std::string& dnaSequence);
    static std::vector<int> findRepeats(const std::string& dnaSequence, int minRepeatLength = 3);
    static std::string generateReport(const DNASequence& seq);
    
private:
    static std::vector<ORF> findORFsInFrame(const std::string& dnaSequence, int frame, int minLength);
    static double nearestNeighborTm(const std::string& sequence);
};

#endif