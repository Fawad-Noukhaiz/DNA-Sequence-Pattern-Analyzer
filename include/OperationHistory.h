#pragma once
#include <string>
#include <ctime>

using namespace std;

struct OperationNode {
    string operation;
    string parameters;
    time_t timestamp;
    OperationNode* next;

    OperationNode(const string& op, const string& params)
        : operation(op), parameters(params), next(nullptr) {
        timestamp = time(nullptr);
    }
};

class OperationHistory {
private:
    OperationNode* head;
    OperationNode* tail;
    int size;
    const int MAX_HISTORY = 50;

public:
    OperationHistory();
    ~OperationHistory();

    void addOperation(const string& operation, const string& parameters = "");
    void displayRecent(int count = 10) const;
    void clear();
    int getSize() const { return size; }
};