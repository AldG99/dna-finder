#include "SequenceAnalyzer.h"
#include "GeneticCode.h"
#include <algorithm>
#include <sstream>
#include <cmath>
#include <map>
#include <set>
#include <iomanip>

std::vector<ORF> SequenceAnalyzer::findORFs(const std::string& dnaSequence, int minLength) {
    return findORFsInFrame(dnaSequence, 1, minLength);
}

std::vector<ORF> SequenceAnalyzer::findORFsAllFrames(const std::string& dnaSequence, int minLength) {
    std::vector<ORF> allOrfs;
    
    for (int frame = 1; frame <= 3; frame++) {
        std::vector<ORF> frameOrfs = findORFsInFrame(dnaSequence, frame, minLength);
        allOrfs.insert(allOrfs.end(), frameOrfs.begin(), frameOrfs.end());
    }
    
    DNASequence dna(dnaSequence);
    std::string reverseComplement = dna.getReverseComplement();
    
    for (int frame = 1; frame <= 3; frame++) {
        std::vector<ORF> frameOrfs = findORFsInFrame(reverseComplement, -frame, minLength);
        allOrfs.insert(allOrfs.end(), frameOrfs.begin(), frameOrfs.end());
    }
    
    std::sort(allOrfs.begin(), allOrfs.end(), [](const ORF& a, const ORF& b) {
        return a.length > b.length;
    });
    
    return allOrfs;
}

std::vector<ORF> SequenceAnalyzer::findORFsInFrame(const std::string& dnaSequence, int frame, int minLength) {
    std::vector<ORF> orfs;
    std::string upperSeq = dnaSequence;
    std::transform(upperSeq.begin(), upperSeq.end(), upperSeq.begin(), ::toupper);
    
    int startFrame = abs(frame) - 1;
    
    for (size_t i = startFrame; i <= upperSeq.length() - 3; i += 3) {
        std::string codon = upperSeq.substr(i, 3);
        
        if (GeneticCode::isStartCodon(codon)) {
            std::string protein = "";
            std::string sequence = "";
            size_t j;
            
            for (j = i; j <= upperSeq.length() - 3; j += 3) {
                std::string currentCodon = upperSeq.substr(j, 3);
                sequence += currentCodon;
                char aminoAcid = GeneticCode::translateCodon(currentCodon);
                protein += aminoAcid;
                
                if (GeneticCode::isStopCodon(currentCodon)) {
                    break;
                }
            }
            
            if (protein.length() >= minLength) {
                orfs.push_back(ORF(i, j + 2, frame, sequence, protein));
            }
            
            if (j > upperSeq.length() - 3 && protein.length() >= minLength) {
                orfs.push_back(ORF(i, upperSeq.length() - 1, frame, sequence, protein));
            }
        }
    }
    
    return orfs;
}

double SequenceAnalyzer::calculateMeltingTemperature(const std::string& dnaSequence) {
    if (dnaSequence.length() < 14) {
        int aCount = 0, tCount = 0, cCount = 0, gCount = 0;
        for (char c : dnaSequence) {
            switch (std::toupper(c)) {
                case 'A': aCount++; break;
                case 'T': tCount++; break;
                case 'C': cCount++; break;
                case 'G': gCount++; break;
            }
        }
        return 2 * (aCount + tCount) + 4 * (cCount + gCount);
    } else {
        return nearestNeighborTm(dnaSequence);
    }
}

double SequenceAnalyzer::nearestNeighborTm(const std::string& sequence) {
    std::map<std::string, double> nnTable = {
        {"AA", -7.9}, {"AT", -7.2}, {"AC", -8.4}, {"AG", -7.8},
        {"TA", -7.2}, {"TT", -7.9}, {"TC", -8.2}, {"TG", -8.5},
        {"CA", -8.5}, {"CT", -7.8}, {"CC", -8.0}, {"CG", -10.6},
        {"GA", -8.2}, {"GT", -8.4}, {"GC", -9.8}, {"GG", -8.0}
    };
    
    double enthalpy = 0.0;
    std::string upperSeq = sequence;
    std::transform(upperSeq.begin(), upperSeq.end(), upperSeq.begin(), ::toupper);
    
    for (size_t i = 0; i < upperSeq.length() - 1; i++) {
        std::string dinucleotide = upperSeq.substr(i, 2);
        if (nnTable.find(dinucleotide) != nnTable.end()) {
            enthalpy += nnTable[dinucleotide];
        }
    }
    
    double entropy = -21.6;
    double saltCorrection = 16.6 * std::log10(0.05);
    
    return (enthalpy * 1000) / (entropy + 1.987 * std::log(1e-6)) - 273.15 + saltCorrection;
}

double SequenceAnalyzer::calculateComplexity(const std::string& dnaSequence) {
    std::map<char, int> counts;
    for (char c : dnaSequence) {
        counts[std::toupper(c)]++;
    }
    
    double entropy = 0.0;
    int total = dnaSequence.length();
    
    for (auto& pair : counts) {
        double probability = static_cast<double>(pair.second) / total;
        if (probability > 0) {
            entropy -= probability * std::log2(probability);
        }
    }
    
    return entropy / 2.0;
}

std::vector<int> SequenceAnalyzer::findRepeats(const std::string& dnaSequence, int minRepeatLength) {
    std::vector<int> repeatPositions;
    std::string upperSeq = dnaSequence;
    std::transform(upperSeq.begin(), upperSeq.end(), upperSeq.begin(), ::toupper);
    
    for (size_t i = 0; i < upperSeq.length() - minRepeatLength; i++) {
        for (int len = minRepeatLength; len <= static_cast<int>(upperSeq.length() - i) / 2; len++) {
            std::string pattern = upperSeq.substr(i, len);
            
            if (upperSeq.substr(i + len, len) == pattern) {
                repeatPositions.push_back(i);
                break;
            }
        }
    }
    
    return repeatPositions;
}

std::string SequenceAnalyzer::generateReport(const DNASequence& seq) {
    std::ostringstream report;
    
    report << "=== ANÁLISIS COMPLETO DE SECUENCIA ===\n\n";
    report << "Secuencia: " << seq.getSequence() << "\n";
    report << "Longitud: " << seq.getLength() << " nucleótidos\n";
    report << "Válida: " << (seq.isValid() ? "Sí" : "No") << "\n\n";
    
    if (!seq.isValid()) {
        report << "La secuencia contiene nucleótidos inválidos.\n";
        return report.str();
    }
    
    report << "--- Secuencias Complementarias ---\n";
    report << "Complementaria: " << seq.getComplement() << "\n";
    report << "Reversa Complementaria: " << seq.getReverseComplement() << "\n\n";
    
    report << "--- Composición ---\n";
    auto counts = seq.getAllCounts();
    for (auto& pair : counts) {
        double percentage = (static_cast<double>(pair.second) / seq.getLength()) * 100;
        report << pair.first << ": " << pair.second << " (" << std::fixed << std::setprecision(2) << percentage << "%)\n";
    }
    report << "Contenido GC: " << std::fixed << std::setprecision(2) << seq.getGCContent() << "%\n";
    report << "Peso Molecular: ~" << std::fixed << std::setprecision(0) << seq.getMolecularWeight() << " Da\n\n";
    
    report << "--- Traducción ---\n";
    std::string protein = GeneticCode::translateSequence(seq.getSequence());
    report << "Proteína (+1): " << protein << "\n\n";
    
    report << "--- Análisis Estadístico ---\n";
    report << "Temperatura de Fusión (Tm): " << std::fixed << std::setprecision(1) 
           << calculateMeltingTemperature(seq.getSequence()) << "°C\n";
    report << "Complejidad: " << std::fixed << std::setprecision(3) 
           << calculateComplexity(seq.getSequence()) << "\n\n";
    
    std::vector<ORF> orfs = findORFsAllFrames(seq.getSequence(), 3);
    report << "--- ORFs Encontrados ---\n";
    if (orfs.empty()) {
        report << "No se encontraron ORFs.\n";
    } else {
        for (size_t i = 0; i < orfs.size() && i < 10; i++) {
            const ORF& orf = orfs[i];
            report << "ORF " << (i+1) << ": Frame " << orf.frame 
                   << ", Posición " << orf.start << "-" << orf.end
                   << ", Longitud " << orf.length << " aa\n";
            report << "  Proteína: " << orf.protein << "\n";
        }
    }
    
    return report.str();
}