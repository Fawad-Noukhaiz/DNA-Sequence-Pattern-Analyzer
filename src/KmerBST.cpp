#include "KmerBST.h"
#include <algorithm>

using namespace std;

KmerBST::KmerBST() : root(nullptr) {}

KmerBST::~KmerBST() {
    clear();
}

KmerNode* KmerBST::insert(KmerNode* node, const string& kmer, int count) {
    if (!node) {
        return new KmerNode(kmer, count);
    }

    if (kmer < node->kmer) {
        node->left = insert(node->left, kmer, count);
    }
    else if (kmer > node->kmer) {
        node->right = insert(node->right, kmer, count);
    }
    else {
        node->count += count;
    }

    return node;
}

void KmerBST::insert(const string& kmer, int count) {
    root = insert(root, kmer, count);
}

KmerNode* KmerBST::findKmer(KmerNode* node, const string& kmer) const {
    if (!node || node->kmer == kmer) {
        return node;
    }

    if (kmer < node->kmer) {
        return findKmer(node->left, kmer);
    }
    else {
        return findKmer(node->right, kmer);
    }
}

bool KmerBST::contains(const string& kmer) const {
    return findKmer(root, kmer) != nullptr;
}

int KmerBST::getCount(const string& kmer) const {
    KmerNode* node = findKmer(root, kmer);
    return node ? node->count : 0;
}

void KmerBST::inOrderTraversal(KmerNode* node,
    vector<pair<string, int>>& result) const
{
    if (!node) return;

    inOrderTraversal(node->left, result);
    result.emplace_back(node->kmer, node->count);
    inOrderTraversal(node->right, result);
}

vector<pair<string, int>> KmerBST::getAllKmers() const {
    vector<pair<string, int>> result;
    inOrderTraversal(root, result);
    return result;
}

void KmerBST::clear(KmerNode* node) {
    if (!node) return;
    clear(node->left);
    clear(node->right);
    delete node;
}

void KmerBST::clear() {
    clear(root);
    root = nullptr;
}