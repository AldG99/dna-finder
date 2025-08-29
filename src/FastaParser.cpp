#include "FastaParser.h"
#include <iostream>
#include <algorithm>
#include <cctype>

std::vector<FastaSequence> FastaParser::parseFile(const std::string& filename) {
    std::vector<FastaSequence> sequences;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo " << filename << std::endl;
        return sequences;
    }
    
    std::string line;
    std::string currentHeader = "";
    std::string currentSequence = "";
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            if (!currentHeader.empty() && !currentSequence.empty()) {
                sequences.push_back(FastaSequence(currentHeader, cleanSequence(currentSequence)));
            }
            
            currentHeader = formatHeader(line);
            currentSequence = "";
        } else {
            currentSequence += line;
        }
    }
    
    if (!currentHeader.empty() && !currentSequence.empty()) {
        sequences.push_back(FastaSequence(currentHeader, cleanSequence(currentSequence)));
    }
    
    file.close();
    return sequences;
}

bool FastaParser::writeFile(const std::string& filename, const std::vector<FastaSequence>& sequences) {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: No se pudo crear el archivo " << filename << std::endl;
        return false;
    }
    
    for (const auto& seq : sequences) {
        file << ">" << seq.header << std::endl;
        
        std::string sequence = seq.sequence;
        for (size_t i = 0; i < sequence.length(); i += 80) {
            file << sequence.substr(i, 80) << std::endl;
        }
    }
    
    file.close();
    return true;
}

FastaSequence FastaParser::parseSingleSequence(const std::string& header, const std::string& sequence) {
    return FastaSequence(formatHeader(header), cleanSequence(sequence));
}

bool FastaParser::isValidFastaFile(const std::string& filename) {
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    bool hasHeader = false;
    bool hasSequence = false;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            hasHeader = true;
        } else {
            hasSequence = true;
            for (char c : line) {
                if (!std::isalpha(c) && !std::isspace(c)) {
                    file.close();
                    return false;
                }
            }
        }
    }
    
    file.close();
    return hasHeader && hasSequence;
}

std::string FastaParser::cleanSequence(const std::string& sequence) {
    std::string cleaned = "";
    
    for (char c : sequence) {
        if (std::isalpha(c)) {
            cleaned += std::toupper(c);
        }
    }
    
    return cleaned;
}

std::string FastaParser::formatHeader(const std::string& header) {
    if (header.empty()) return "Unnamed_sequence";
    
    if (header[0] == '>') {
        return header.substr(1);
    }
    
    return header;
}