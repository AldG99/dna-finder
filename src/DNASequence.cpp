#include "DNASequence.h"
#include <algorithm>
#include <cctype>
#include <stdexcept>

DNASequence::DNASequence() : sequence(""), valid(true) {
    nucleotideCount = {{'A', 0}, {'T', 0}, {'C', 0}, {'G', 0}};
}

DNASequence::DNASequence(const std::string& seq) : sequence(""), valid(false) {
    nucleotideCount = {{'A', 0}, {'T', 0}, {'C', 0}, {'G', 0}};
    setSequence(seq);
}

void DNASequence::setSequence(const std::string& seq) {
    sequence = seq;
    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
    validateSequence();
    if (valid) {
        countNucleotides();
    }
}

std::string DNASequence::getSequence() const {
    return sequence;
}

bool DNASequence::isValid() const {
    return valid;
}

void DNASequence::validateSequence() {
    valid = true;
    for (char nucleotide : sequence) {
        if (!isValidNucleotide(nucleotide)) {
            valid = false;
            break;
        }
    }
}

void DNASequence::countNucleotides() {
    nucleotideCount = {{'A', 0}, {'T', 0}, {'C', 0}, {'G', 0}};
    for (char nucleotide : sequence) {
        if (nucleotideCount.find(nucleotide) != nucleotideCount.end()) {
            nucleotideCount[nucleotide]++;
        }
    }
}

std::string DNASequence::getComplement() const {
    if (!valid) return "";
    
    std::string complement = "";
    for (char nucleotide : sequence) {
        complement += getComplementNucleotide(nucleotide);
    }
    return complement;
}

std::string DNASequence::getReverseComplement() const {
    std::string complement = getComplement();
    std::reverse(complement.begin(), complement.end());
    return complement;
}

double DNASequence::getGCContent() const {
    if (!valid || sequence.empty()) return 0.0;
    
    int gcCount = nucleotideCount.at('G') + nucleotideCount.at('C');
    return (static_cast<double>(gcCount) / sequence.length()) * 100.0;
}

int DNASequence::getNucleotideCount(char nucleotide) const {
    char upperNucleotide = std::toupper(nucleotide);
    auto it = nucleotideCount.find(upperNucleotide);
    return (it != nucleotideCount.end()) ? it->second : 0;
}

std::map<char, int> DNASequence::getAllCounts() const {
    return nucleotideCount;
}

double DNASequence::getMolecularWeight() const {
    if (!valid) return 0.0;
    
    const double weights[] = {331.2, 322.2, 307.2, 347.2}; // A, T, C, G
    const char nucleotides[] = {'A', 'T', 'C', 'G'};
    
    double totalWeight = 0.0;
    for (int i = 0; i < 4; i++) {
        totalWeight += nucleotideCount.at(nucleotides[i]) * weights[i];
    }
    
    return totalWeight - (sequence.length() - 1) * 18.01528; // Subtract water molecules
}

int DNASequence::getLength() const {
    return sequence.length();
}

bool DNASequence::isEmpty() const {
    return sequence.empty();
}

bool DNASequence::isValidNucleotide(char nucleotide) {
    char upper = std::toupper(nucleotide);
    return upper == 'A' || upper == 'T' || upper == 'C' || upper == 'G' || 
           upper == 'N' || upper == 'R' || upper == 'Y' || upper == 'K' || 
           upper == 'M' || upper == 'S' || upper == 'W' || upper == 'B' || 
           upper == 'D' || upper == 'H' || upper == 'V';
}

char DNASequence::getComplementNucleotide(char nucleotide) {
    switch (std::toupper(nucleotide)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'N': return 'N';
        case 'R': return 'Y';
        case 'Y': return 'R';
        case 'K': return 'M';
        case 'M': return 'K';
        case 'S': return 'S';
        case 'W': return 'W';
        case 'B': return 'V';
        case 'V': return 'B';
        case 'D': return 'H';
        case 'H': return 'D';
        default: return 'N';
    }
}