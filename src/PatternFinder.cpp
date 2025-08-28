#include "PatternFinder.h"
#include <algorithm>
#include <cctype>

std::vector<PatternMatch> PatternFinder::findPattern(const std::string& sequence, const std::string& pattern) {
    std::vector<PatternMatch> matches;
    std::string upperSeq = sequence;
    std::string upperPat = pattern;
    
    std::transform(upperSeq.begin(), upperSeq.end(), upperSeq.begin(), ::toupper);
    std::transform(upperPat.begin(), upperPat.end(), upperPat.begin(), ::toupper);
    
    size_t pos = 0;
    while ((pos = upperSeq.find(upperPat, pos)) != std::string::npos) {
        std::string matched = sequence.substr(pos, upperPat.length());
        matches.push_back(PatternMatch(pos, pattern, matched));
        pos++;
    }
    
    return matches;
}

std::vector<PatternMatch> PatternFinder::findPatternWithWildcards(const std::string& sequence, const std::string& pattern) {
    std::vector<PatternMatch> matches;
    std::string upperSeq = sequence;
    std::string upperPat = pattern;
    
    std::transform(upperSeq.begin(), upperSeq.end(), upperSeq.begin(), ::toupper);
    std::transform(upperPat.begin(), upperPat.end(), upperPat.begin(), ::toupper);
    
    for (size_t i = 0; i <= upperSeq.length() - upperPat.length(); i++) {
        if (matchesWithWildcards(upperSeq, upperPat, i)) {
            std::string matched = sequence.substr(i, upperPat.length());
            matches.push_back(PatternMatch(i, pattern, matched));
        }
    }
    
    return matches;
}

std::vector<PatternMatch> PatternFinder::findRestrictionSites(const std::string& sequence) {
    std::vector<PatternMatch> allMatches;
    auto sites = getCommonRestrictionSites();
    
    for (const auto& site : sites) {
        std::vector<PatternMatch> matches = findPatternWithWildcards(sequence, site.second);
        for (auto& match : matches) {
            match.pattern = site.first + " (" + site.second + ")";
            allMatches.push_back(match);
        }
    }
    
    std::sort(allMatches.begin(), allMatches.end(), 
              [](const PatternMatch& a, const PatternMatch& b) {
                  return a.position < b.position;
              });
    
    return allMatches;
}

std::vector<PatternMatch> PatternFinder::findPrimers(const std::string& sequence, const std::string& primer) {
    std::vector<PatternMatch> matches;
    
    std::vector<PatternMatch> forwardMatches = findPattern(sequence, primer);
    for (auto& match : forwardMatches) {
        match.pattern = "Forward: " + match.pattern;
        matches.push_back(match);
    }
    
    std::string reverseComplement = "";
    for (char c : primer) {
        switch (std::toupper(c)) {
            case 'A': reverseComplement = 'T' + reverseComplement; break;
            case 'T': reverseComplement = 'A' + reverseComplement; break;
            case 'C': reverseComplement = 'G' + reverseComplement; break;
            case 'G': reverseComplement = 'C' + reverseComplement; break;
            default: reverseComplement = 'N' + reverseComplement; break;
        }
    }
    
    std::vector<PatternMatch> reverseMatches = findPattern(sequence, reverseComplement);
    for (auto& match : reverseMatches) {
        match.pattern = "Reverse: " + primer + " (RC: " + reverseComplement + ")";
        matches.push_back(match);
    }
    
    return matches;
}

std::vector<PatternMatch> PatternFinder::findAllMatches(const std::string& sequence, const std::vector<std::string>& patterns) {
    std::vector<PatternMatch> allMatches;
    
    for (const std::string& pattern : patterns) {
        std::vector<PatternMatch> matches = findPatternWithWildcards(sequence, pattern);
        allMatches.insert(allMatches.end(), matches.begin(), matches.end());
    }
    
    std::sort(allMatches.begin(), allMatches.end(), 
              [](const PatternMatch& a, const PatternMatch& b) {
                  return a.position < b.position;
              });
    
    return allMatches;
}

std::map<std::string, std::string> PatternFinder::getCommonRestrictionSites() {
    return {
        {"EcoRI", "GAATTC"},
        {"BamHI", "GGATCC"},
        {"HindIII", "AAGCTT"},
        {"XbaI", "TCTAGA"},
        {"SalI", "GTCGAC"},
        {"PstI", "CTGCAG"},
        {"SmaI", "CCCGGG"},
        {"KpnI", "GGTACC"},
        {"SacI", "GAGCTC"},
        {"XhoI", "CTCGAG"},
        {"NdeI", "CATATG"},
        {"NcoI", "CCATGG"},
        {"BglII", "AGATCT"},
        {"ApaI", "GGGCCC"},
        {"NotI", "GCGGCCGC"}
    };
}

bool PatternFinder::matchesWithWildcards(const std::string& sequence, const std::string& pattern, size_t pos) {
    for (size_t i = 0; i < pattern.length(); i++) {
        char seqChar = sequence[pos + i];
        char patChar = pattern[i];
        
        if (patChar != 'N' && patChar != seqChar) {
            switch (patChar) {
                case 'R': if (seqChar != 'A' && seqChar != 'G') return false; break;
                case 'Y': if (seqChar != 'C' && seqChar != 'T') return false; break;
                case 'K': if (seqChar != 'G' && seqChar != 'T') return false; break;
                case 'M': if (seqChar != 'A' && seqChar != 'C') return false; break;
                case 'S': if (seqChar != 'C' && seqChar != 'G') return false; break;
                case 'W': if (seqChar != 'A' && seqChar != 'T') return false; break;
                case 'B': if (seqChar == 'A') return false; break;
                case 'D': if (seqChar == 'C') return false; break;
                case 'H': if (seqChar == 'G') return false; break;
                case 'V': if (seqChar == 'T') return false; break;
                default: return false;
            }
        }
    }
    return true;
}