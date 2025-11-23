#include "KmerAnalyzer.h"
#include <algorithm>
#include <queue>

using namespace std;

unordered_map<string, int> KmerAnalyzer::count(const string& seq, int k) {
    unordered_map<string, int> kmerCounts;

    if (k <= 0 || k > static_cast<int>(seq.size())) {
        return kmerCounts;
    }

    for (size_t i = 0; i <= seq.size() - k; i++) {
        string kmer = seq.substr(i, k);
        kmerCounts[kmer]++;
    }

    return kmerCounts;
}

vector<pair<string, int>> KmerAnalyzer::topKmers(
    const unordered_map<string, int>& kmers, int n)
{
    vector<pair<string, int>> sortedKmers(kmers.begin(), kmers.end());

    sort(sortedKmers.begin(), sortedKmers.end(),
        [](const pair<string, int>& a, const pair<string, int>& b) {
            return a.second > b.second;
        });

    if (n < static_cast<int>(sortedKmers.size())) {
        sortedKmers.resize(n);
    }

    return sortedKmers;
}

vector<pair<string, int>> KmerAnalyzer::topKmersHeap(
    const unordered_map<string, int>& kmers, int n)
{
    using KmerPair = pair<string, int>;

    auto comparator = [](const KmerPair& a, const KmerPair& b) {
        return a.second > b.second;
        };

    priority_queue<KmerPair, vector<KmerPair>, decltype(comparator)> minHeap(comparator);

    for (const auto& entry : kmers) {
        if (minHeap.size() < n) {
            minHeap.push(entry);
        }
        else if (entry.second > minHeap.top().second) {
            minHeap.pop();
            minHeap.push(entry);
        }
    }

    vector<KmerPair> result;
    while (!minHeap.empty()) {
        result.push_back(minHeap.top());
        minHeap.pop();
    }

    reverse(result.begin(), result.end());
    return result;
}