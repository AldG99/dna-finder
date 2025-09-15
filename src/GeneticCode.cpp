#include "GeneticCode.h"
#include <algorithm>
#include <cctype>

std::map<std::string, char> GeneticCode::codonTable;
std::map<char, std::string> GeneticCode::codonToAminoAcid;
bool GeneticCode::tableInitialized = false;

void GeneticCode::initializeCodonTable() {
    if (tableInitialized) return;
    
    codonTable["TTT"] = 'F'; codonTable["TTC"] = 'F';
    codonTable["TTA"] = 'L'; codonTable["TTG"] = 'L';
    codonTable["TCT"] = 'S'; codonTable["TCC"] = 'S'; codonTable["TCA"] = 'S'; codonTable["TCG"] = 'S';
    codonTable["TAT"] = 'Y'; codonTable["TAC"] = 'Y';
    codonTable["TAA"] = '*'; codonTable["TAG"] = '*';
    codonTable["TGT"] = 'C'; codonTable["TGC"] = 'C';
    codonTable["TGA"] = '*';
    codonTable["TGG"] = 'W';
    
    codonTable["CTT"] = 'L'; codonTable["CTC"] = 'L'; codonTable["CTA"] = 'L'; codonTable["CTG"] = 'L';
    codonTable["CCT"] = 'P'; codonTable["CCC"] = 'P'; codonTable["CCA"] = 'P'; codonTable["CCG"] = 'P';
    codonTable["CAT"] = 'H'; codonTable["CAC"] = 'H';
    codonTable["CAA"] = 'Q'; codonTable["CAG"] = 'Q';
    codonTable["CGT"] = 'R'; codonTable["CGC"] = 'R'; codonTable["CGA"] = 'R'; codonTable["CGG"] = 'R';
    
    codonTable["ATT"] = 'I'; codonTable["ATC"] = 'I'; codonTable["ATA"] = 'I';
    codonTable["ATG"] = 'M';
    codonTable["ACT"] = 'T'; codonTable["ACC"] = 'T'; codonTable["ACA"] = 'T'; codonTable["ACG"] = 'T';
    codonTable["AAT"] = 'N'; codonTable["AAC"] = 'N';
    codonTable["AAA"] = 'K'; codonTable["AAG"] = 'K';
    codonTable["AGT"] = 'S'; codonTable["AGC"] = 'S';
    codonTable["AGA"] = 'R'; codonTable["AGG"] = 'R';
    
    codonTable["GTT"] = 'V'; codonTable["GTC"] = 'V'; codonTable["GTA"] = 'V'; codonTable["GTG"] = 'V';
    codonTable["GCT"] = 'A'; codonTable["GCC"] = 'A'; codonTable["GCA"] = 'A'; codonTable["GCG"] = 'A';
    codonTable["GAT"] = 'D'; codonTable["GAC"] = 'D';
    codonTable["GAA"] = 'E'; codonTable["GAG"] = 'E';
    codonTable["GGT"] = 'G'; codonTable["GGC"] = 'G'; codonTable["GGA"] = 'G'; codonTable["GGG"] = 'G';
    
    codonToAminoAcid['F'] = "Phenylalanine";
    codonToAminoAcid['L'] = "Leucine";
    codonToAminoAcid['S'] = "Serine";
    codonToAminoAcid['Y'] = "Tyrosine";
    codonToAminoAcid['*'] = "Stop";
    codonToAminoAcid['C'] = "Cysteine";
    codonToAminoAcid['W'] = "Tryptophan";
    codonToAminoAcid['P'] = "Proline";
    codonToAminoAcid['H'] = "Histidine";
    codonToAminoAcid['Q'] = "Glutamine";
    codonToAminoAcid['R'] = "Arginine";
    codonToAminoAcid['I'] = "Isoleucine";
    codonToAminoAcid['M'] = "Methionine";
    codonToAminoAcid['T'] = "Threonine";
    codonToAminoAcid['N'] = "Asparagine";
    codonToAminoAcid['K'] = "Lysine";
    codonToAminoAcid['V'] = "Valine";
    codonToAminoAcid['A'] = "Alanine";
    codonToAminoAcid['D'] = "Aspartic acid";
    codonToAminoAcid['E'] = "Glutamic acid";
    codonToAminoAcid['G'] = "Glycine";
    
    tableInitialized = true;
}

char GeneticCode::translateCodon(const std::string& codon) {
    initializeCodonTable();
    
    if (codon.length() != 3) return 'X';
    
    std::string upperCodon = codon;
    std::transform(upperCodon.begin(), upperCodon.end(), upperCodon.begin(), ::toupper);
    
    auto it = codonTable.find(upperCodon);
    return (it != codonTable.end()) ? it->second : 'X';
}

std::string GeneticCode::translateSequence(const std::string& dnaSequence) {
    initializeCodonTable();
    
    std::string protein = "";
    std::string upperSeq = dnaSequence;
    std::transform(upperSeq.begin(), upperSeq.end(), upperSeq.begin(), ::toupper);
    
    for (size_t i = 0; i <= upperSeq.length() - 3; i += 3) {
        std::string codon = upperSeq.substr(i, 3);
        char aminoAcid = translateCodon(codon);
        protein += aminoAcid;
        
        if (aminoAcid == '*') break;
    }
    
    return protein;
}

std::string GeneticCode::translateSequenceVerbose(const std::string& dnaSequence) {
    initializeCodonTable();
    
    std::string result = "";
    std::string upperSeq = dnaSequence;
    std::transform(upperSeq.begin(), upperSeq.end(), upperSeq.begin(), ::toupper);
    
    for (size_t i = 0; i <= upperSeq.length() - 3; i += 3) {
        std::string codon = upperSeq.substr(i, 3);
        char aminoAcid = translateCodon(codon);
        
        result += codon + " -> " + aminoAcid + " (" + getAminoAcidName(aminoAcid) + ")\n";
        
        if (aminoAcid == '*') break;
    }
    
    return result;
}

bool GeneticCode::isStartCodon(const std::string& codon) {
    std::string upperCodon = codon;
    std::transform(upperCodon.begin(), upperCodon.end(), upperCodon.begin(), ::toupper);
    return upperCodon == "ATG";
}

bool GeneticCode::isStopCodon(const std::string& codon) {
    std::string upperCodon = codon;
    std::transform(upperCodon.begin(), upperCodon.end(), upperCodon.begin(), ::toupper);
    return upperCodon == "TAA" || upperCodon == "TAG" || upperCodon == "TGA";
}

std::string GeneticCode::getAminoAcidName(char aminoAcid) {
    initializeCodonTable();
    
    auto it = codonToAminoAcid.find(aminoAcid);
    return (it != codonToAminoAcid.end()) ? it->second : "Unknown";
}

std::vector<int> GeneticCode::findStartCodons(const std::string& dnaSequence) {
    std::vector<int> positions;
    std::string upperSeq = dnaSequence;
    std::transform(upperSeq.begin(), upperSeq.end(), upperSeq.begin(), ::toupper);
    
    for (size_t i = 0; i <= upperSeq.length() - 3; i++) {
        if (upperSeq.substr(i, 3) == "ATG") {
            positions.push_back(i);
        }
    }
    
    return positions;
}

std::vector<int> GeneticCode::findStopCodons(const std::string& dnaSequence) {
    std::vector<int> positions;
    std::string upperSeq = dnaSequence;
    std::transform(upperSeq.begin(), upperSeq.end(), upperSeq.begin(), ::toupper);
    
    for (size_t i = 0; i <= upperSeq.length() - 3; i++) {
        std::string codon = upperSeq.substr(i, 3);
        if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
            positions.push_back(i);
        }
    }
    
    return positions;
}