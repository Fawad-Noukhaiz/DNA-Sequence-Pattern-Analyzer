#pragma once
#include <string>

using namespace std;

class SequenceLoader {
public:
    static bool loadFASTA(const string& filename, string& outSeq, string& outHeader);
};