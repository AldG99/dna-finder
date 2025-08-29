#ifndef PATTERNFINDER_H
#define PATTERNFINDER_H

#include <string>
#include <vector>
#include <map>

struct PatternMatch {
    int position;
    std::string pattern;
    std::string matchedSequence;
    
    PatternMatch(int pos, const std::string& pat, const std::string& seq)
        : position(pos), pattern(pat), matchedSequence(seq) {}
};

class PatternFinder {
public:
    static std::vector<PatternMatch> findPattern(const std::string& sequence, const std::string& pattern);
    static std::vector<PatternMatch> findPatternWithWildcards(const std::string& sequence, const std::string& pattern);
    static std::vector<PatternMatch> findRestrictionSites(const std::string& sequence);
    static std::vector<PatternMatch> findPrimers(const std::string& sequence, const std::string& primer);
    static std::vector<PatternMatch> findAllMatches(const std::string& sequence, const std::vector<std::string>& patterns);
    
    static std::map<std::string, std::string> getCommonRestrictionSites();
    
private:
    static bool matchesWithWildcards(const std::string& sequence, const std::string& pattern, size_t pos);
    static char wildcardToRegex(char wildcard);
};

#endif