#pragma once
#include <string>
#include <vector>

using namespace std;

struct KmerNode {
    string kmer;
    int count;
    KmerNode* left;
    KmerNode* right;

    KmerNode(const string& k, int c)
        : kmer(k), count(c), left(nullptr), right(nullptr) {
    }
};

class KmerBST {
private:
    KmerNode* root;

    KmerNode* insert(KmerNode* node, const string& kmer, int count);
    void clear(KmerNode* node);
    void inOrderTraversal(KmerNode* node, vector<pair<string, int>>& result) const;
    KmerNode* findKmer(KmerNode* node, const string& kmer) const;

public:
    KmerBST();
    ~KmerBST();

    void insert(const string& kmer, int count);
    bool contains(const string& kmer) const;
    int getCount(const string& kmer) const;
    vector<pair<string, int>> getAllKmers() const;
    void clear();
};