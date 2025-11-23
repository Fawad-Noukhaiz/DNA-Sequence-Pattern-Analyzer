#pragma once
#include <string>
#include <unordered_map>
#include <vector>
#include <queue>
#include <functional>

using namespace std;

class KmerAnalyzer {
public:
    static unordered_map<string, int> count(const string& seq, int k);
    static vector<pair<string, int>> topKmers(const unordered_map<string, int>& kmers, int n);
    static vector<pair<string, int>> topKmersHeap(const unordered_map<string, int>& kmers, int n);
};