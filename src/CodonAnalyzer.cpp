#include "CodonAnalyzer.h"
#include "GeneticCode.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sstream>

CodonAnalyzer::CodonAnalyzer() {
    // Initialize genetic code
    m_geneticCode = {
        {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
        {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
        {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
        {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
        
        {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
        {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        
        {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
        {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        
        {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
        {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };
    
    // Initialize amino acid to codons mapping
    for (const auto& pair : m_geneticCode) {
        char aa = pair.second;
        if (aa != '*') {  // Skip stop codons for most analyses
            m_aminoAcidToCodons[aa].push_back(pair.first);
        }
    }
    
    initializeOrganismTables();
}

CodonAnalyzer::~CodonAnalyzer() {
    // Destructor
}

void CodonAnalyzer::initializeOrganismTables() {
    initializeEColiTable();
    initializeYeastTable(); 
    initializeHumanTable();
    initializePlantTable();
}

void CodonAnalyzer::initializeEColiTable() {
    OrganismCodonTable ecoli;
    ecoli.organismName = "E.coli";
    
    // E.coli codon frequencies (per thousand) - from highly expressed genes
    ecoli.codonFrequencies = {
        {"TTT", 22.0}, {"TTC", 16.8}, {"TTA", 13.5}, {"TTG", 13.0},
        {"TCT", 15.2}, {"TCC", 8.8}, {"TCA", 7.8}, {"TCG", 14.4},
        {"TAT", 16.2}, {"TAC", 12.2}, {"TAA", 2.0}, {"TAG", 0.2},
        {"TGT", 5.2}, {"TGC", 6.2}, {"TGA", 1.0}, {"TGG", 15.2},
        
        {"CTT", 11.2}, {"CTC", 10.8}, {"CTA", 3.8}, {"CTG", 52.6},
        {"CCT", 7.2}, {"CCC", 5.8}, {"CCA", 8.8}, {"CCG", 23.0},
        {"CAT", 13.2}, {"CAC", 9.8}, {"CAA", 15.2}, {"CAG", 29.2},
        {"CGT", 38.4}, {"CGC", 22.2}, {"CGA", 3.8}, {"CGG", 5.8},
        
        {"ATT", 30.2}, {"ATC", 25.2}, {"ATA", 4.8}, {"ATG", 27.2},
        {"ACT", 15.2}, {"ACC", 25.2}, {"ACA", 7.2}, {"ACG", 14.8},
        {"AAT", 17.2}, {"AAC", 22.2}, {"AAA", 33.2}, {"AAG", 10.8},
        {"AGT", 15.2}, {"AGC", 16.2}, {"AGA", 2.2}, {"AGG", 1.8},
        
        {"GTT", 18.2}, {"GTC", 20.8}, {"GTA", 11.2}, {"GTG", 26.2},
        {"GCT", 18.8}, {"GCC", 27.2}, {"GCA", 21.2}, {"GCG", 33.8},
        {"GAT", 32.2}, {"GAC", 19.2}, {"GAA", 39.2}, {"GAG", 18.8},
        {"GGT", 24.8}, {"GGC", 29.2}, {"GGA", 8.8}, {"GGG", 11.2}
    };
    
    // Optimal codons for E.coli (weights for CAI calculation)
    ecoli.optimalCodons = {
        {"CTG", 1.0}, {"CGT", 1.0}, {"GCG", 1.0}, {"GAA", 1.0},
        {"TTC", 1.0}, {"GGC", 1.0}, {"CAC", 1.0}, {"ATC", 1.0},
        {"AAG", 1.0}, {"TTG", 1.0}, {"ATG", 1.0}, {"AAC", 1.0},
        {"CCG", 1.0}, {"CAG", 1.0}, {"AGC", 1.0}, {"ACC", 1.0},
        {"GTG", 1.0}, {"TGG", 1.0}, {"TAC", 1.0}
    };
    
    m_organismTables["E.coli"] = ecoli;
}

void CodonAnalyzer::initializeYeastTable() {
    OrganismCodonTable yeast;
    yeast.organismName = "S.cerevisiae";
    
    // Yeast codon frequencies - optimized for this organism
    yeast.codonFrequencies = {
        {"TTT", 26.1}, {"TTC", 18.4}, {"TTA", 28.1}, {"TTG", 27.2},
        {"TCT", 26.2}, {"TCC", 16.8}, {"TCA", 21.8}, {"TCG", 8.8},
        {"TAT", 19.2}, {"TAC", 14.8}, {"TAA", 1.1}, {"TAG", 0.5},
        {"TGT", 8.1}, {"TGC", 4.8}, {"TGA", 0.7}, {"TGG", 10.4},
        
        {"CTT", 12.3}, {"CTC", 5.4}, {"CTA", 14.2}, {"CTG", 10.5},
        {"CCT", 13.5}, {"CCC", 6.8}, {"CCA", 18.2}, {"CCG", 5.3},
        {"CAT", 13.8}, {"CAC", 7.8}, {"CAA", 27.3}, {"CAG", 12.1},
        {"CGT", 6.4}, {"CGC", 2.6}, {"CGA", 3.0}, {"CGG", 1.7},
        
        {"ATT", 30.1}, {"ATC", 17.2}, {"ATA", 17.8}, {"ATG", 20.9},
        {"ACT", 20.3}, {"ACC", 12.7}, {"ACA", 18.2}, {"ACG", 8.0},
        {"AAT", 35.8}, {"AAC", 24.8}, {"AAA", 42.0}, {"AAG", 30.8},
        {"AGT", 14.2}, {"AGC", 9.8}, {"AGA", 21.3}, {"AGG", 9.2},
        
        {"GTT", 22.1}, {"GTC", 11.8}, {"GTA", 12.1}, {"GTG", 10.8},
        {"GCT", 21.2}, {"GCC", 12.6}, {"GCA", 16.2}, {"GCG", 6.2},
        {"GAT", 37.8}, {"GAC", 20.2}, {"GAA", 45.6}, {"GAG", 19.2},
        {"GGT", 24.0}, {"GGC", 9.8}, {"GGA", 10.8}, {"GGG", 6.2}
    };
    
    m_organismTables["yeast"] = yeast;
    m_organismTables["S.cerevisiae"] = yeast;
}

void CodonAnalyzer::initializeHumanTable() {
    OrganismCodonTable human;
    human.organismName = "Human";
    
    // Human codon frequencies
    human.codonFrequencies = {
        {"TTT", 17.2}, {"TTC", 20.4}, {"TTA", 7.2}, {"TTG", 12.8},
        {"TCT", 15.2}, {"TCC", 17.8}, {"TCA", 12.2}, {"TCG", 4.8},
        {"TAT", 12.2}, {"TAC", 15.8}, {"TAA", 0.7}, {"TAG", 0.6},
        {"TGT", 10.2}, {"TGC", 12.8}, {"TGA", 1.3}, {"TGG", 13.2},
        
        {"CTT", 13.2}, {"CTC", 19.8}, {"CTA", 7.2}, {"CTG", 39.8},
        {"CCT", 17.8}, {"CCC", 19.8}, {"CCA", 16.8}, {"CCG", 6.8},
        {"CAT", 10.8}, {"CAC", 15.2}, {"CAA", 12.2}, {"CAG", 34.2},
        {"CGT", 4.8}, {"CGC", 10.8}, {"CGA", 6.2}, {"CGG", 11.8},
        
        {"ATT", 16.2}, {"ATC", 21.2}, {"ATA", 7.2}, {"ATG", 22.2},
        {"ACT", 13.2}, {"ACC", 18.8}, {"ACA", 15.2}, {"ACG", 6.2},
        {"AAT", 17.2}, {"AAC", 19.2}, {"AAA", 24.2}, {"AAG", 32.8},
        {"AGT", 12.2}, {"AGC", 19.2}, {"AGA", 12.2}, {"AGG", 12.2},
        
        {"GTT", 11.2}, {"GTC", 14.8}, {"GTA", 7.2}, {"GTG", 28.2},
        {"GCT", 18.8}, {"GCC", 27.8}, {"GCA", 15.8}, {"GCG", 7.2},
        {"GAT", 22.2}, {"GAC", 25.8}, {"GAA", 29.2}, {"GAG", 40.8},
        {"GGT", 16.8}, {"GGC", 22.2}, {"GGA", 16.2}, {"GGG", 16.2}
    };
    
    m_organismTables["human"] = human;
    m_organismTables["Human"] = human;
}

void CodonAnalyzer::initializePlantTable() {
    OrganismCodonTable plant;
    plant.organismName = "A.thaliana";
    
    // Arabidopsis thaliana codon frequencies
    plant.codonFrequencies = {
        {"TTT", 22.4}, {"TTC", 18.8}, {"TTA", 8.8}, {"TTG", 13.8},
        {"TCT", 18.4}, {"TCC", 14.8}, {"TCA", 13.8}, {"TCG", 11.8},
        {"TAT", 15.2}, {"TAC", 13.8}, {"TAA", 1.2}, {"TAG", 0.8},
        {"TGT", 12.8}, {"TGC", 9.8}, {"TGA", 1.8}, {"TGG", 12.8},
        
        {"CTT", 16.8}, {"CTC", 14.8}, {"CTA", 8.8}, {"CTG", 24.8},
        {"CCT", 16.8}, {"CCC", 13.8}, {"CCA", 17.8}, {"CCG", 8.8},
        {"CAT", 14.8}, {"CAC", 12.8}, {"CAA", 18.8}, {"CAG", 22.8},
        {"CGT", 8.8}, {"CGC", 7.8}, {"CGA", 7.8}, {"CGG", 6.8},
        
        {"ATT", 19.8}, {"ATC", 16.8}, {"ATA", 9.8}, {"ATG", 23.8},
        {"ACT", 16.8}, {"ACC", 15.8}, {"ACA", 16.8}, {"ACG", 9.8},
        {"AAT", 19.8}, {"AAC", 17.8}, {"AAA", 26.8}, {"AAG", 25.8},
        {"AGT", 14.8}, {"AGC", 12.8}, {"AGA", 14.8}, {"AGG", 10.8},
        
        {"GTT", 16.8}, {"GTC", 13.8}, {"GTA", 9.8}, {"GTG", 22.8},
        {"GCT", 22.8}, {"GCC", 18.8}, {"GCA", 17.8}, {"GCG", 9.8},
        {"GAT", 25.8}, {"GAC", 19.8}, {"GAA", 32.8}, {"GAG", 26.8},
        {"GGT", 19.8}, {"GGC", 16.8}, {"GGA", 17.8}, {"GGG", 12.8}
    };
    
    m_organismTables["plant"] = plant;
    m_organismTables["A.thaliana"] = plant;
}

CodonAnalysisReport CodonAnalyzer::analyzeCodonUsage(const std::string& sequence, const std::string& organism) {
    CodonAnalysisReport report;
    
    // Get codon counts
    std::map<std::string, int> codonCounts = getCodonCounts(sequence);
    report.totalCodons = 0;
    for (const auto& pair : codonCounts) {
        report.totalCodons += pair.second;
    }
    
    if (report.totalCodons == 0) {
        report.expressionPrediction = "Cannot analyze empty sequence";
        return report;
    }
    
    // Calculate GC content
    int gcCount = 0;
    for (char c : sequence) {
        if (c == 'G' || c == 'C') gcCount++;
    }
    report.gcContent = (double)gcCount / sequence.length() * 100.0;
    
    // Calculate CAI
    report.caiScore = calculateCAI(sequence, organism);
    
    // Create codon usage data
    for (const auto& pair : codonCounts) {
        CodonUsageData usage;
        usage.codon = pair.first;
        usage.count = pair.second;
        usage.frequency = (double)pair.second / report.totalCodons * 100.0;
        
        if (m_geneticCode.find(pair.first) != m_geneticCode.end()) {
            usage.aminoAcid = m_geneticCode[pair.first];
            usage.rscu = calculateRSCU(pair.first, codonCounts);
        } else {
            usage.aminoAcid = '?';
            usage.rscu = 0.0;
        }
        
        report.codonUsage.push_back(usage);
        report.codonsByAminoAcid[usage.aminoAcid].push_back(usage);
    }
    
    // Sort by frequency
    std::sort(report.codonUsage.begin(), report.codonUsage.end(),
              [](const CodonUsageData& a, const CodonUsageData& b) {
                  return a.frequency > b.frequency;
              });
    
    // Calculate advanced metrics
    report.enc = calculateENC(codonCounts);
    report.codonBias = calculateCodonBias(codonCounts, organism);
    
    // Find rare codons
    std::vector<std::string> rareCodons = findRareCodons(sequence, organism, 0.05);
    for (const std::string& rareCodon : rareCodons) {
        if (codonCounts.find(rareCodon) != codonCounts.end()) {
            report.rareCodonCount[rareCodon] = codonCounts[rareCodon];
        }
    }
    
    // Expression prediction
    report.expressionPrediction = predictExpressionLevel(sequence, organism);
    
    // Optimization recommendations
    report.recommendedOptimizations = getOptimizationSuggestions(sequence, organism);
    
    return report;
}

double CodonAnalyzer::calculateCAI(const std::string& sequence, const std::string& organism) {
    std::vector<std::string> codons = sequenceToCodons(sequence);
    if (codons.empty()) return 0.0;
    
    auto it = m_organismTables.find(organism);
    if (it == m_organismTables.end()) {
        it = m_organismTables.find("E.coli");  // Default to E.coli
    }
    
    const OrganismCodonTable& table = it->second;
    double logSum = 0.0;
    int validCodons = 0;
    
    for (const std::string& codon : codons) {
        if (m_geneticCode.find(codon) != m_geneticCode.end() && m_geneticCode[codon] != '*') {
            char aa = m_geneticCode[codon];
            
            // Find maximum frequency for this amino acid
            double maxFreq = 0.0;
            for (const std::string& synonymCodon : m_aminoAcidToCodons[aa]) {
                auto freqIt = table.codonFrequencies.find(synonymCodon);
                if (freqIt != table.codonFrequencies.end()) {
                    maxFreq = std::max(maxFreq, freqIt->second);
                }
            }
            
            // Calculate relative weight
            auto freqIt = table.codonFrequencies.find(codon);
            if (freqIt != table.codonFrequencies.end() && maxFreq > 0) {
                double weight = freqIt->second / maxFreq;
                logSum += std::log(weight);
                validCodons++;
            }
        }
    }
    
    return validCodons > 0 ? std::exp(logSum / validCodons) : 0.0;
}

std::string CodonAnalyzer::predictExpressionLevel(const std::string& sequence, const std::string& organism) {
    double cai = calculateCAI(sequence, organism);
    std::vector<std::string> rareCodons = findRareCodons(sequence, organism);
    
    return classifyExpression(cai, rareCodons.size());
}

std::string CodonAnalyzer::classifyExpression(double caiScore, int rareCodonCount) {
    std::string level;
    std::string details;
    
    // Classify based on CAI score
    if (caiScore >= 0.8) {
        level = "Alto";
    } else if (caiScore >= 0.6) {
        level = "Medio";
    } else if (caiScore >= 0.4) {
        level = "Bajo";
    } else {
        level = "Muy Bajo";
    }
    
    // Adjust based on rare codons
    if (rareCodonCount > 10) {
        level = "Bajo (muchos codones raros)";
    } else if (rareCodonCount > 5) {
        if (level == "Alto") level = "Medio-Alto";
        else if (level == "Medio") level = "Medio-Bajo";
    }
    
    return level + " (CAI: " + std::to_string(caiScore).substr(0, 4) + ")";
}

std::vector<std::string> CodonAnalyzer::getOptimizationSuggestions(const std::string& sequence, const std::string& targetOrganism) {
    std::vector<std::string> suggestions;
    
    double cai = calculateCAI(sequence, targetOrganism);
    std::vector<std::string> rareCodons = findRareCodons(sequence, targetOrganism);
    
    if (cai < 0.6) {
        suggestions.push_back("• CAI bajo (" + std::to_string(cai).substr(0, 4) + ") - considerar optimización de codones");
    }
    
    if (!rareCodons.empty()) {
        suggestions.push_back("• Encontrados " + std::to_string(rareCodons.size()) + " codones raros - podrían limitar expresión");
        
        if (rareCodons.size() <= 3) {
            suggestions.push_back("• Codones raros encontrados: ");
            for (const std::string& codon : rareCodons) {
                if (m_geneticCode.find(codon) != m_geneticCode.end()) {
                    suggestions.push_back("  - " + codon + " (" + std::string(1, m_geneticCode[codon]) + ")");
                }
            }
        }
    }
    
    if (cai > 0.8 && rareCodons.size() < 3) {
        suggestions.push_back("• Secuencia bien optimizada para " + targetOrganism);
    }
    
    return suggestions;
}

std::map<std::string, int> CodonAnalyzer::getCodonCounts(const std::string& sequence) {
    std::map<std::string, int> counts;
    std::vector<std::string> codons = sequenceToCodons(sequence);
    
    for (const std::string& codon : codons) {
        if (isValidCodon(codon)) {
            counts[codon]++;
        }
    }
    
    return counts;
}

std::vector<std::string> CodonAnalyzer::findRareCodons(const std::string& sequence, const std::string& organism, double threshold) {
    std::vector<std::string> rareCodons;
    std::vector<std::string> codons = sequenceToCodons(sequence);
    
    auto it = m_organismTables.find(organism);
    if (it == m_organismTables.end()) {
        it = m_organismTables.find("E.coli");
    }
    
    const OrganismCodonTable& table = it->second;
    
    for (const std::string& codon : codons) {
        auto freqIt = table.codonFrequencies.find(codon);
        if (freqIt != table.codonFrequencies.end()) {
            if (freqIt->second / 1000.0 < threshold) {  // Convert per thousand to fraction
                if (std::find(rareCodons.begin(), rareCodons.end(), codon) == rareCodons.end()) {
                    rareCodons.push_back(codon);
                }
            }
        }
    }
    
    return rareCodons;
}

double CodonAnalyzer::calculateRSCU(const std::string& codon, const std::map<std::string, int>& codonCounts) {
    if (m_geneticCode.find(codon) == m_geneticCode.end()) {
        return 0.0;
    }
    
    char aa = m_geneticCode[codon];
    if (aa == '*') return 0.0;  // Skip stop codons
    
    // Count total usage for this amino acid
    int totalForAA = 0;
    int synonymousCodons = 0;
    
    for (const std::string& synonymCodon : m_aminoAcidToCodons[aa]) {
        auto it = codonCounts.find(synonymCodon);
        if (it != codonCounts.end()) {
            totalForAA += it->second;
        }
        synonymousCodons++;
    }
    
    if (totalForAA == 0) return 0.0;
    
    auto it = codonCounts.find(codon);
    int codonCount = (it != codonCounts.end()) ? it->second : 0;
    
    double expectedFreq = (double)totalForAA / synonymousCodons;
    return expectedFreq > 0 ? (double)codonCount / expectedFreq : 0.0;
}

double CodonAnalyzer::calculateENC(const std::map<std::string, int>& codonCounts) {
    // Effective Number of Codons - measures codon bias
    // Higher values = less bias
    
    std::map<char, std::vector<int>> aaCodonCounts;
    
    // Group codons by amino acid
    for (const auto& pair : codonCounts) {
        if (m_geneticCode.find(pair.first) != m_geneticCode.end()) {
            char aa = m_geneticCode[pair.first];
            if (aa != '*') {
                aaCodonCounts[aa].push_back(pair.second);
            }
        }
    }
    
    double enc = 0.0;
    int groups = 0;
    
    for (const auto& pair : aaCodonCounts) {
        const std::vector<int>& counts = pair.second;
        int total = 0;
        for (int count : counts) {
            total += count;
        }
        
        if (total > 0) {
            double homozygosity = 0.0;
            for (int count : counts) {
                double freq = (double)count / total;
                homozygosity += freq * freq;
            }
            
            if (homozygosity > 0) {
                enc += 1.0 / homozygosity;
                groups++;
            }
        }
    }
    
    return groups > 0 ? enc / groups : 0.0;
}

double CodonAnalyzer::calculateCodonBias(const std::map<std::string, int>& codonCounts, const std::string& organism) {
    // Simple codon bias metric based on deviation from expected usage
    auto it = m_organismTables.find(organism);
    if (it == m_organismTables.end()) {
        it = m_organismTables.find("E.coli");
    }
    
    const OrganismCodonTable& table = it->second;
    
    double bias = 0.0;
    int totalCodons = 0;
    for (const auto& pair : codonCounts) {
        totalCodons += pair.second;
    }
    
    if (totalCodons == 0) return 0.0;
    
    for (const auto& pair : codonCounts) {
        double observedFreq = (double)pair.second / totalCodons;
        
        auto expectedIt = table.codonFrequencies.find(pair.first);
        if (expectedIt != table.codonFrequencies.end()) {
            double expectedFreq = expectedIt->second / 1000.0;  // Convert per thousand
            bias += std::abs(observedFreq - expectedFreq);
        }
    }
    
    return bias;
}

std::string CodonAnalyzer::generateCodonReport(const CodonAnalysisReport& report) {
    std::stringstream ss;
    
    ss << "=== ANÁLISIS PROFESIONAL DE USO DE CODONES ===\n\n";
    ss << "Organismo objetivo: E.coli (por defecto)\n";
    ss << "Total de codones: " << report.totalCodons << "\n";
    ss << "Contenido GC: " << std::fixed << std::setprecision(1) << report.gcContent << "%\n";
    ss << "Índice de Adaptación de Codones (CAI): " << std::fixed << std::setprecision(3) << report.caiScore << "\n";
    ss << "Predicción de expresión: " << report.expressionPrediction << "\n";
    ss << "Número Efectivo de Codones (ENC): " << std::fixed << std::setprecision(1) << report.enc << "\n\n";
    
    ss << "=== USO DE CODONES (Top 10) ===\n";
    ss << std::left << std::setw(8) << "Codón" << std::setw(4) << "AA" 
       << std::setw(8) << "Count" << std::setw(10) << "Freq%" 
       << std::setw(8) << "RSCU" << "\n";
    ss << std::string(40, '-') << "\n";
    
    for (size_t i = 0; i < std::min((size_t)10, report.codonUsage.size()); i++) {
        const CodonUsageData& usage = report.codonUsage[i];
        ss << std::left << std::setw(8) << usage.codon 
           << std::setw(4) << usage.aminoAcid
           << std::setw(8) << usage.count 
           << std::setw(10) << std::fixed << std::setprecision(1) << usage.frequency
           << std::setw(8) << std::fixed << std::setprecision(2) << usage.rscu << "\n";
    }
    
    if (!report.rareCodonCount.empty()) {
        ss << "\n=== CODONES RAROS DETECTADOS ===\n";
        for (const auto& pair : report.rareCodonCount) {
            ss << "• " << pair.first << ": " << pair.second << " veces\n";
        }
    }
    
    if (!report.recommendedOptimizations.empty()) {
        ss << "\n=== RECOMENDACIONES DE OPTIMIZACIÓN ===\n";
        for (const std::string& rec : report.recommendedOptimizations) {
            ss << rec << "\n";
        }
    }
    
    ss << "\n=== INTERPRETACIÓN PROFESIONAL ===\n";
    ss << "• CAI > 0.8: Excelente para expresión alta\n";
    ss << "• CAI 0.6-0.8: Bueno para expresión moderada\n";
    ss << "• CAI < 0.6: Considerar optimización\n";
    ss << "• ENC > 45: Baja preferencia de codones\n";
    ss << "• ENC < 35: Alta preferencia de codones\n";
    
    return ss.str();
}

std::vector<std::string> CodonAnalyzer::getSupportedOrganisms() {
    std::vector<std::string> organisms;
    for (const auto& pair : m_organismTables) {
        organisms.push_back(pair.first);
    }
    return organisms;
}

bool CodonAnalyzer::isValidCodon(const std::string& codon) {
    return codon.length() == 3 && 
           std::all_of(codon.begin(), codon.end(), [](char c) {
               return c == 'A' || c == 'T' || c == 'C' || c == 'G';
           });
}

std::vector<std::string> CodonAnalyzer::sequenceToCodons(const std::string& sequence) {
    std::vector<std::string> codons;
    
    for (size_t i = 0; i + 2 < sequence.length(); i += 3) {
        std::string codon = sequence.substr(i, 3);
        if (isValidCodon(codon)) {
            codons.push_back(codon);
        }
    }
    
    return codons;
}

std::map<char, std::vector<std::string>> CodonAnalyzer::getCodonsByAminoAcid() {
    return m_aminoAcidToCodons;
}