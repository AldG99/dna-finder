#ifndef CODONANALYZER_H
#define CODONANALYZER_H

#include <string>
#include <map>
#include <vector>
#include <unordered_map>

// Structure to hold codon usage data
struct CodonUsageData {
    std::string codon;
    char aminoAcid;
    int count;
    double frequency;          // Frequency within all codons
    double relativeFrequency;  // Frequency within codons for same amino acid
    double rscu;              // Relative Synonymous Codon Usage
};

// Structure for organism-specific codon preferences
struct OrganismCodonTable {
    std::string organismName;
    std::map<std::string, double> codonFrequencies;  // Codon -> frequency
    std::map<std::string, double> optimalCodons;     // Codon -> weight for CAI
};

// Structure for codon analysis results
struct CodonAnalysisReport {
    int totalCodons;
    double gcContent;
    double caiScore;                              // Codon Adaptation Index
    std::string expressionPrediction;             // High/Medium/Low
    std::vector<CodonUsageData> codonUsage;
    std::map<char, std::vector<CodonUsageData>> codonsByAminoAcid;
    std::vector<std::string> recommendedOptimizations;
    
    // Advanced metrics
    double enc;                                   // Effective Number of Codons
    double codonBias;                            // Codon bias index
    std::map<std::string, int> rareCodonCount;   // Rare codons found
};

class CodonAnalyzer {
public:
    CodonAnalyzer();
    ~CodonAnalyzer();

    // Main analysis functions
    CodonAnalysisReport analyzeCodonUsage(const std::string& sequence, 
                                         const std::string& organism = "E.coli");
    
    // Codon Adaptation Index calculation
    double calculateCAI(const std::string& sequence, const std::string& organism = "E.coli");
    
    // Expression prediction based on codon usage
    std::string predictExpressionLevel(const std::string& sequence, 
                                     const std::string& organism = "E.coli");
    
    // Codon optimization suggestions
    std::vector<std::string> getOptimizationSuggestions(const std::string& sequence,
                                                       const std::string& targetOrganism = "E.coli");
    
    // Optimize sequence for better expression
    std::string optimizeSequence(const std::string& sequence, 
                                const std::string& targetOrganism = "E.coli");
    
    // Utility functions
    std::map<std::string, int> getCodonCounts(const std::string& sequence);
    std::map<char, std::vector<std::string>> getCodonsByAminoAcid();
    std::vector<std::string> getSupportedOrganisms();
    
    // Rare codon analysis
    std::vector<std::string> findRareCodons(const std::string& sequence,
                                          const std::string& organism = "E.coli",
                                          double threshold = 0.05);
    
    // Generate detailed report
    std::string generateCodonReport(const CodonAnalysisReport& report);
    
    // Static utility functions
    static bool isValidCodon(const std::string& codon);
    static std::vector<std::string> sequenceToCodons(const std::string& sequence);
    
private:
    // Organism codon tables
    std::map<std::string, OrganismCodonTable> m_organismTables;
    
    // Genetic code mapping
    std::map<std::string, char> m_geneticCode;
    std::map<char, std::vector<std::string>> m_aminoAcidToCodons;
    
    // Initialize codon tables for different organisms
    void initializeOrganismTables();
    void initializeEColiTable();
    void initializeYeastTable();
    void initializeHumanTable();
    void initializePlantTable();
    
    // Helper functions
    double calculateRSCU(const std::string& codon, 
                        const std::map<std::string, int>& codonCounts);
    double calculateENC(const std::map<std::string, int>& codonCounts);
    double calculateCodonBias(const std::map<std::string, int>& codonCounts,
                             const std::string& organism);
    
    // Expression prediction helpers
    std::string classifyExpression(double caiScore, int rareCodonCount);
    std::vector<std::string> generateOptimizationTips(const CodonAnalysisReport& report);
    
    // Sequence optimization
    std::string chooseOptimalCodon(char aminoAcid, const std::string& organism);
    bool shouldOptimizeCodon(const std::string& codon, const std::string& organism);
};

#endif // CODONANALYZER_H