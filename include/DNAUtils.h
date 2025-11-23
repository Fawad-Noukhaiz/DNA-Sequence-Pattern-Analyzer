#pragma once
#include <string>
#include <bitset>

using namespace std;

class DNAUtils {
public:
    static double gcContent(const string& seq);
    static bool containsSRY(const string& seq);
    static string reverseComplement(const string& seq);
    static bool isValidDNA(const string& seq);
    static bool quickValidation(const string& seq);
    static size_t sequenceHash(const string& seq);
};