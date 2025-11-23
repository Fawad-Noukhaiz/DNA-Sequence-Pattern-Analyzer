#include "OperationHistory.h"
#include <iostream>
#include <iomanip>

using namespace std;

OperationHistory::OperationHistory() : head(nullptr), tail(nullptr), size(0) {}

OperationHistory::~OperationHistory() {
    clear();
}

void OperationHistory::addOperation(const string& operation, const string& parameters) {
    OperationNode* newNode = new OperationNode(operation, parameters);

    if (!head) {
        head = tail = newNode;
    }
    else {
        tail->next = newNode;
        tail = newNode;
    }

    size++;

    if (size > MAX_HISTORY) {
        OperationNode* toDelete = head;
        head = head->next;
        delete toDelete;
        size--;
    }
}

void OperationHistory::displayRecent(int count) const {
    if (size == 0) {
        cout << "No operations recorded.\n";
        return;
    }

    cout << "\nRecent Operations (showing last " << min(count, size) << "):\n";
    cout << "--------------------------------------------------\n";

    OperationNode* current = head;
    int startFrom = max(0, size - count);

    for (int i = 0; i < startFrom && current; i++) {
        current = current->next;
    }

    int index = startFrom + 1;
    while (current) {
        tm timeinfo;
        localtime_s(&timeinfo, &current->timestamp);

        cout << setw(2) << index << ". ["
            << put_time(&timeinfo, "%H:%M:%S") << "] "
            << current->operation;

        if (!current->parameters.empty()) {
            cout << " - " << current->parameters;
        }
        cout << "\n";

        current = current->next;
        index++;
    }
}

void OperationHistory::clear() {
    while (head) {
        OperationNode* temp = head;
        head = head->next;
        delete temp;
    }
    head = tail = nullptr;
    size = 0;
}