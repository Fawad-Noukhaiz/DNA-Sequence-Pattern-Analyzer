#pragma once
#include <string>
#include <vector>

using namespace std;

class PatternSearch {
public:
    static vector<int> kmp(const string& text, const string& pat);
    static vector<int> boyerMoore(const string& text, const string& pat);
    static vector<int> rabinKarp(const string& text, const string& pat);
    static vector<int> naiveSearch(const string& text, const string& pat);
    static vector<string> getAlgorithmNames();
};