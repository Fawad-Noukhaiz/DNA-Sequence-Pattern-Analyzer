#include "SequenceLoader.h"
#include <fstream>
#include <iostream>
#include <cctype>
#include <filesystem>

using namespace std;

bool SequenceLoader::loadFASTA(const string& filename, string& outSeq, string& outHeader) {
    try {
        filesystem::path p(filename);

        if (!filesystem::exists(p)) {
            cerr << "Error: File does not exist: " << p << '\n';
            return false;
        }

        uintmax_t filesize = 0;
        try {
            filesize = filesystem::file_size(p);
        }
        catch (...) {}

        ifstream file(p, ios::binary);
        if (!file.is_open()) {
            cerr << "Error: Could not open file: " << p << '\n';
            return false;
        }

        outSeq.clear();
        outHeader.clear();

        if (filesize > 0 && filesize < (uintmax_t)(4ULL * 1024 * 1024 * 1024)) {
            outSeq.reserve(static_cast<size_t>(min<uintmax_t>(filesize, (uintmax_t)(size_t(-1)))));
        }

        cout << "Loading FASTA file: " << p << " ...\n";

        string line;
        bool headerRead = false;
        size_t invalidChars = 0;
        size_t lineNum = 0;

        while (getline(file, line)) {
            lineNum++;

            if (!line.empty() && line.back() == '\r') line.pop_back();

            if (line.empty()) continue;

            if (line[0] == '>') {
                if (headerRead) {
                    cerr << "Warning: Multiple sequences in file. Loading only first sequence.\n";
                    break;
                }
                outHeader = line.substr(1);
                headerRead = true;
            }
            else {
                for (unsigned char uc : line) {
                    if (isspace(uc)) continue;
                    char base = static_cast<char>(toupper(uc));
                    if (base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == 'N') {
                        outSeq.push_back(base);
                    }
                    else {
                        invalidChars++;
                        outSeq.push_back('N');
                    }
                }
            }

            if (lineNum % 1000000 == 0) {
                cout << "  Processed " << lineNum << " lines, "
                    << (outSeq.size() / 1000000) << " Mbp\n";
            }
        }

        file.close();

        if (!headerRead) {
            cerr << "Error: No header line found\n";
            return false;
        }

        if (outSeq.empty()) {
            cerr << "Error: No sequence data found\n";
            return false;
        }

        if (invalidChars > 0) {
            cerr << "Warning: Replaced " << invalidChars << " invalid characters with 'N'\n";
        }

        cout << "Successfully loaded " << outSeq.size() << " base pairs\n";
        return true;
    }
    catch (const exception& e) {
        cerr << "Exception while loading file: " << e.what() << '\n';
        return false;
    }
}