#include "PatternSearch.h"
#include <unordered_map>
#include <cmath>
#include <chrono>

using namespace std;

static vector<int> buildLPS(const string& pat) {
    vector<int> lps(pat.size(), 0);
    int len = 0;

    for (int i = 1; i < pat.size(); ) {
        if (pat[i] == pat[len]) {
            lps[i++] = ++len;
        }
        else {
            if (len != 0) len = lps[len - 1];
            else lps[i++] = 0;
        }
    }
    return lps;
}

vector<int> PatternSearch::kmp(const string& text, const string& pat) {
    vector<int> result;
    if (pat.empty() || text.empty() || pat.size() > text.size())
        return result;

    auto lps = buildLPS(pat);
    int i = 0, j = 0;

    while (i < text.size()) {
        if (text[i] == pat[j]) {
            i++; j++;
            if (j == pat.size()) {
                result.push_back(i - j);
                j = lps[j - 1];
            }
        }
        else {
            if (j != 0) j = lps[j - 1];
            else i++;
        }
    }
    return result;
}

vector<int> PatternSearch::boyerMoore(const string& text, const string& pat) {
    vector<int> result;
    if (pat.empty() || text.empty() || pat.size() > text.size())
        return result;

    unordered_map<char, int> badChar;
    int patLen = pat.size();
    int textLen = text.size();

    for (int i = 0; i < patLen; i++) {
        badChar[pat[i]] = i;
    }

    int shift = 0;
    while (shift <= (textLen - patLen)) {
        int j = patLen - 1;

        while (j >= 0 && pat[j] == text[shift + j]) {
            j--;
        }

        if (j < 0) {
            result.push_back(shift);
            shift += (shift + patLen < textLen) ? patLen - badChar[text[shift + patLen]] : 1;
        }
        else {
            int badCharShift = j - badChar[text[shift + j]];
            shift += max(1, badCharShift);
        }
    }
    return result;
}

vector<int> PatternSearch::naiveSearch(const string& text, const string& pat) {
    vector<int> result;
    if (pat.empty() || text.empty() || pat.size() > text.size())
        return result;

    for (int i = 0; i <= text.size() - pat.size(); i++) {
        bool found = true;
        for (int j = 0; j < pat.size(); j++) {
            if (text[i + j] != pat[j]) {
                found = false;
                break;
            }
        }
        if (found) result.push_back(i);
    }
    return result;
}

vector<int> PatternSearch::rabinKarp(const string& text, const string& pat) {
    vector<int> result;
    if (pat.empty() || text.empty() || pat.size() > text.size())
        return result;

    const int prime = 101;
    const int base = 256;
    int patLen = pat.size();
    int textLen = text.size();

    int patHash = 0;
    int textHash = 0;
    int h = 1;

    for (int i = 0; i < patLen - 1; i++) {
        h = (h * base) % prime;
    }

    for (int i = 0; i < patLen; i++) {
        patHash = (base * patHash + pat[i]) % prime;
        textHash = (base * textHash + text[i]) % prime;
    }

    for (int i = 0; i <= textLen - patLen; i++) {
        if (patHash == textHash) {
            bool match = true;
            for (int j = 0; j < patLen; j++) {
                if (text[i + j] != pat[j]) {
                    match = false;
                    break;
                }
            }
            if (match) result.push_back(i);
        }

        if (i < textLen - patLen) {
            textHash = (base * (textHash - text[i] * h) + text[i + patLen]) % prime;
            if (textHash < 0) textHash += prime;
        }
    }

    return result;
}

vector<string> PatternSearch::getAlgorithmNames() {
    return {
        "KMP (Knuth-Morris-Pratt)",
        "Boyer-Moore",
        "Rabin-Karp",
        "Naive Search"
    };
}