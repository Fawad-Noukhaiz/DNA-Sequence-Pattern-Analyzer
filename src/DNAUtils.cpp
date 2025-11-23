#include "DNAUtils.h"
#include "PatternSearch.h"
#include <algorithm>
#include <iostream>
#include <functional>

using namespace std;

double DNAUtils::gcContent(const string& seq) {
    if (seq.empty()) return 0.0;

    size_t gc = 0;
    for (char c : seq) {
        if (c == 'G' || c == 'C') gc++;
    }

    return (gc * 100.0) / seq.size();
}

bool DNAUtils::containsSRY(const string& seq) {
    string marker = "TCCAGTTTTGTTACAGGG";
    auto found = PatternSearch::kmp(seq, marker);
    return !found.empty();
}

string DNAUtils::reverseComplement(const string& seq) {
    string result;
    result.reserve(seq.size());

    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        switch (*it) {
        case 'A': result += 'T'; break;
        case 'T': result += 'A'; break;
        case 'C': result += 'G'; break;
        case 'G': result += 'C'; break;
        case 'N': result += 'N'; break;
        default:  result += 'N'; break;
        }
    }

    return result;
}

bool DNAUtils::isValidDNA(const string& seq) {
    for (char c : seq) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
            return false;
        }
    }
    return true;
}

bool DNAUtils::quickValidation(const string& seq) {
    static const bitset<256> validChars = []() {
        bitset<256> bs;
        bs['A'] = bs['C'] = bs['G'] = bs['T'] = bs['N'] = true;
        bs['a'] = bs['c'] = bs['g'] = bs['t'] = bs['n'] = true;
        return bs;
        }();

    for (char c : seq) {
        if (!validChars[static_cast<unsigned char>(c)]) {
            return false;
        }
    }
    return true;
}

size_t DNAUtils::sequenceHash(const string& seq) {
    size_t hash = 0;
    const size_t prime = 31;

    for (char c : seq) {
        int value = 0;
        switch (toupper(c)) {
        case 'A': value = 1; break;
        case 'C': value = 2; break;
        case 'G': value = 3; break;
        case 'T': value = 4; break;
        default: value = 5; break;
        }
        hash = hash * prime + value;
    }
    return hash;
}