#ifndef FASTAPARSER_H
#define FASTAPARSER_H

#include <string>
#include <vector>
#include <fstream>

struct FastaSequence {
    std::string header;
    std::string sequence;
    
    FastaSequence(const std::string& h, const std::string& s) : header(h), sequence(s) {}
};

class FastaParser {
public:
    static std::vector<FastaSequence> parseFile(const std::string& filename);
    static bool writeFile(const std::string& filename, const std::vector<FastaSequence>& sequences);
    static FastaSequence parseSingleSequence(const std::string& header, const std::string& sequence);
    static bool isValidFastaFile(const std::string& filename);
    
private:
    static std::string cleanSequence(const std::string& sequence);
    static std::string formatHeader(const std::string& header);
};

#endif