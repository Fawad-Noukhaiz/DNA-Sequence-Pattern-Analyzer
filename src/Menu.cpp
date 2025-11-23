#include "Menu.h"
#include "SequenceLoader.h"
#include "PatternSearch.h"
#include "KmerAnalyzer.h"
#include "DNAUtils.h"
#include "OperationHistory.h"
#include "KmerBST.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>

using namespace std;

static void showMenu(bool loaded, const string& header) {
    cout << "\n==== DNA Analyzer ====\n";
    if (!loaded) {
        cout << "[No FASTA file loaded]\n";
    }
    else {
        cout << "[Loaded: " << header << "]\n";
    }

    cout << "1) Load FASTA file\n";
    cout << "2) Pattern Search\n";
    cout << "3) K-mer Count\n";
    cout << "4) GC Content\n";
    cout << "5) SRY Detection\n";
    cout << "6) Sequence Info\n";
    cout << "7) Operation History\n";
    cout << "8) Validate Sequence\n";
    cout << "9) Toggle K-mer Algorithm\n";
    cout << "10) Exit\n";
    cout << "Choose: ";
}

void Menu::run(int argc, char* argv[]) {
    string sequence;
    string header;
    bool loaded = false;
    OperationHistory history;
    bool useHeapForKmers = false;

    if (argc > 1) {
        string path = argv[1];
        if (SequenceLoader::loadFASTA(path, sequence, header)) {
            loaded = true;
            history.addOperation("Load FASTA", path);
        }
    }

    while (true) {
        showMenu(loaded, header);

        int choice;
        if (!(cin >> choice)) {
            cin.clear();
            cin.ignore(99999, '\n');
            continue;
        }

        switch (choice) {
        case 1: {
            cout << "Enter FASTA filename: ";
            string path;
            cin >> path;

            if (SequenceLoader::loadFASTA(path, sequence, header)) {
                loaded = true;
                history.addOperation("Load FASTA", path);
            }
            break;
        }

        case 2: {
            if (!loaded) {
                cout << "Please load a FASTA first.\n";
                break;
            }

            cout << "Enter pattern to search: ";
            string pat;
            cin >> pat;

            for (char& c : pat) c = toupper(c);

            cout << "\nSelect search algorithm:\n";
            vector<string> algorithms = PatternSearch::getAlgorithmNames();
            for (size_t i = 0; i < algorithms.size(); i++) {
                cout << "  " << (i + 1) << ") " << algorithms[i] << "\n";
            }
            cout << "  " << (algorithms.size() + 1) << ") Compare all algorithms\n";
            cout << "Choose: ";

            int algoChoice;
            if (!(cin >> algoChoice)) {
                cin.clear();
                cin.ignore(99999, '\n');
                algoChoice = 1;
            }

            vector<int> positions;
            string algoName;
            auto start = chrono::high_resolution_clock::now();
            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> duration;

            if (algoChoice >= 1 && algoChoice <= static_cast<int>(algorithms.size())) {
                cout << "\nSearching for pattern '" << pat << "' using "
                    << algorithms[algoChoice - 1] << "...\n";

                start = chrono::high_resolution_clock::now();
                switch (algoChoice) {
                case 1: positions = PatternSearch::kmp(sequence, pat); break;
                case 2: positions = PatternSearch::boyerMoore(sequence, pat); break;
                case 3: positions = PatternSearch::rabinKarp(sequence, pat); break;
                case 4: positions = PatternSearch::naiveSearch(sequence, pat); break;
                }
                end = chrono::high_resolution_clock::now();
                duration = end - start;

                algoName = algorithms[algoChoice - 1];

                cout << "\nFound " << positions.size() << " matches in "
                    << fixed << setprecision(6) << duration.count() << " seconds.\n";

                if (!positions.empty()) {
                    int toShow = min(10, static_cast<int>(positions.size()));
                    cout << "First " << toShow << " positions:\n";
                    for (int i = 0; i < toShow; i++) {
                        cout << "  " << positions[i] << "\n";
                    }
                    if (positions.size() > 10) {
                        cout << "  ... and " << (positions.size() - 10) << " more\n";
                    }
                }

                history.addOperation("Pattern Search",
                    "Pattern: " + pat + ", Algorithm: " + algoName +
                    ", Found: " + to_string(positions.size()) +
                    ", Time: " + to_string(duration.count()) + "s");
            }
            else if (algoChoice == static_cast<int>(algorithms.size() + 1)) {
                cout << "\n=== Comparing All Search Algorithms ===\n";
                cout << "Pattern: '" << pat << "' in sequence of length "
                    << sequence.size() << "\n\n";

                vector<pair<string, vector<int>>> results;
                vector<pair<string, double>> timings;

                for (size_t i = 0; i < algorithms.size(); i++) {
                    cout << "Running " << algorithms[i] << "...";
                    cout.flush();

                    start = chrono::high_resolution_clock::now();
                    vector<int> algoPositions;
                    switch (i) {
                    case 0: algoPositions = PatternSearch::kmp(sequence, pat); break;
                    case 1: algoPositions = PatternSearch::boyerMoore(sequence, pat); break;
                    case 2: algoPositions = PatternSearch::rabinKarp(sequence, pat); break;
                    case 3: algoPositions = PatternSearch::naiveSearch(sequence, pat); break;
                    }
                    end = chrono::high_resolution_clock::now();
                    duration = end - start;

                    results.push_back({ algorithms[i], algoPositions });
                    timings.push_back({ algorithms[i], duration.count() });

                    cout << " found " << algoPositions.size() << " matches, "
                        << fixed << setprecision(6) << duration.count() << " seconds\n";
                }

                bool consistent = true;
                size_t expected = results[0].second.size();
                for (size_t i = 1; i < results.size(); i++) {
                    if (results[i].second.size() != expected) {
                        consistent = false;
                        break;
                    }
                }

                cout << "\n=== Results Summary ===\n";
                if (!consistent) {
                    cout << "WARNING: Algorithms found different numbers of matches!\n";
                }

                sort(timings.begin(), timings.end(),
                    [](const pair<string, double>& a, const pair<string, double>& b) {
                        return a.second < b.second;
                    });

                cout << "\nPerformance Ranking:\n";
                for (size_t i = 0; i < timings.size(); i++) {
                    cout << "  " << (i + 1) << ". " << timings[i].first << ": "
                        << fixed << setprecision(6) << timings[i].second << " seconds\n";
                }

                history.addOperation("Pattern Search Compare",
                    "Pattern: " + pat + ", Best: " + timings[0].first +
                    ", Time: " + to_string(timings[0].second) + "s");
            }
            else {
                cout << "Invalid algorithm choice. Using KMP algorithm.\n";

                start = chrono::high_resolution_clock::now();
                positions = PatternSearch::kmp(sequence, pat);
                end = chrono::high_resolution_clock::now();
                duration = end - start;

                cout << "\nFound " << positions.size() << " matches in "
                    << fixed << setprecision(6) << duration.count() << " seconds.\n";

                history.addOperation("Pattern Search",
                    "Pattern: " + pat + ", Algorithm: KMP" +
                    ", Found: " + to_string(positions.size()) +
                    ", Time: " + to_string(duration.count()) + "s");
            }
            break;
        }

        case 3: {
            if (!loaded) {
                cout << "Please load a FASTA first.\n";
                break;
            }

            cout << "Enter k (recommended: 3-15 for large chromosomes): ";
            int k;
            cin >> k;

            auto kmers = KmerAnalyzer::count(sequence, k);

            if (!kmers.empty()) {
                cout << "\nTop 10 most frequent " << k << "-mers:\n";

                vector<pair<string, int>> top;
                if (useHeapForKmers) {
                    top = KmerAnalyzer::topKmersHeap(kmers, 10);
                    cout << "[Using Heap Algorithm]\n";
                }
                else {
                    top = KmerAnalyzer::topKmers(kmers, 10);
                    cout << "[Using Sorting Algorithm]\n";
                }

                for (size_t i = 0; i < top.size(); i++) {
                    cout << setw(3) << (i + 1) << ". "
                        << top[i].first << " : " << top[i].second << " times\n";
                }

                cout << "\nBuilding K-mer BST for demonstration...\n";
                KmerBST bst;
                for (const auto& kmer : kmers) {
                    bst.insert(kmer.first, kmer.second);
                }
                cout << "BST contains 'ATG': " << (bst.contains("ATG") ? "Yes" : "No") << "\n";
                if (k >= 3 && bst.contains("ATG")) {
                    cout << "Count of 'ATG' in BST: " << bst.getCount("ATG") << "\n";
                }
            }

            history.addOperation("K-mer Analysis",
                "k=" + to_string(k) + ", Algorithm=" +
                (useHeapForKmers ? "Heap" : "Sorting"));
            break;
        }

        case 4: {
            if (!loaded) {
                cout << "Please load a FASTA first.\n";
                break;
            }

            cout << "Calculating GC content...\n";
            double gc = DNAUtils::gcContent(sequence);
            cout << "\nGC Content: " << fixed << setprecision(2) << gc << "%\n";

            if (gc < 40) cout << "(Low GC content)\n";
            else if (gc > 55) cout << "(High GC content)\n";
            else cout << "(Moderate GC content)\n";

            history.addOperation("GC Content", "Result: " + to_string(gc) + "%");
            break;
        }

        case 5: {
            if (!loaded) {
                cout << "Please load a FASTA first.\n";
                break;
            }

            cout << "Searching for SRY marker...\n";
            bool yes = DNAUtils::containsSRY(sequence);

            if (yes) {
                cout << "\nSRY gene marker found.\nLikely MALE.\n";
            }
            else {
                cout << "\nSRY gene marker not found.\nLikely FEMALE or no Y chromosome.\n";
            }

            history.addOperation("SRY Detection", yes ? "Yes" : "No");
            break;
        }

        case 6: {
            if (!loaded) {
                cout << "Please load a FASTA first.\n";
                break;
            }

            cout << "\n=== Sequence Information ===\n";
            cout << "Header: " << header << "\n";
            cout << "Length: " << sequence.size() << " bp";
            cout << " (" << sequence.size() / 1000000.0 << " Mbp)\n";

            size_t countA = 0, countC = 0, countG = 0, countT = 0, countN = 0;
            for (char c : sequence) {
                switch (c) {
                case 'A': countA++; break;
                case 'C': countC++; break;
                case 'G': countG++; break;
                case 'T': countT++; break;
                case 'N': countN++; break;
                }
            }

            cout << "\nBase composition:\n";
            cout << "  A: " << countA << " (" << (countA * 100.0 / sequence.size()) << "%)\n";
            cout << "  C: " << countC << " (" << (countC * 100.0 / sequence.size()) << "%)\n";
            cout << "  G: " << countG << " (" << (countG * 100.0 / sequence.size()) << "%)\n";
            cout << "  T: " << countT << " (" << (countT * 100.0 / sequence.size()) << "%)\n";
            if (countN > 0) {
                cout << "  N: " << countN << " (" << (countN * 100.0 / sequence.size()) << "%)\n";
            }

            history.addOperation("Sequence Info",
                "Length: " + to_string(sequence.size()));
            break;
        }

        case 7:
            history.displayRecent();
            break;

        case 8: {
            if (!loaded) {
                cout << "Please load a FASTA first.\n";
                break;
            }

            cout << "Validating sequence...\n";
            bool isValid = DNAUtils::isValidDNA(sequence);
            bool quickValid = DNAUtils::quickValidation(sequence);
            size_t hash = DNAUtils::sequenceHash(sequence);

            cout << "\nValidation Results:\n";
            cout << "Full Validation: " << (isValid ? "VALID" : "INVALID") << "\n";
            cout << "Quick Validation: " << (quickValid ? "VALID" : "INVALID") << "\n";
            cout << "Sequence Hash: " << hash << "\n";

            history.addOperation("Sequence Validation",
                "Full: " + string(isValid ? "Valid" : "Invalid") +
                ", Quick: " + string(quickValid ? "Valid" : "Invalid"));
            break;
        }

        case 9:
            useHeapForKmers = !useHeapForKmers;
            cout << "K-mer algorithm set to: "
                << (useHeapForKmers ? "HEAP (Priority Queue)" : "SORTING") << "\n";
            history.addOperation("Toggle Algorithm",
                useHeapForKmers ? "Heap" : "Sorting");
            break;

        case 10:
            cout << "Goodbye!\n";
            return;

        default:
            cout << "Invalid option.\n";
        }
    }
}