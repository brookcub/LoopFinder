#include "Utils.h"

std::vector<unsigned> strideSubsample(const std::vector<unsigned>& v, const unsigned consecutiveElements, const unsigned rejectWindow) {
    // expects that v is in ascending order
    std::vector<unsigned> subsampledElements;

    unsigned consecutiveCount = 0;
    unsigned lastSelectedValue = 0;

    for (unsigned i = 0; i < v.size(); i++) {
        unsigned elementValue = v[i];

        // Keep the first consecutiveElements Elements
        if (consecutiveCount < consecutiveElements) {
            subsampledElements.push_back(elementValue);
            consecutiveCount++;
            lastSelectedValue = elementValue;
        } else {
            // Check if the candidate is outside the reject window
            if (elementValue > lastSelectedValue + rejectWindow) {
                subsampledElements.push_back(elementValue);
                lastSelectedValue = elementValue;
                consecutiveCount = 1;
            }
        }
    }

    return subsampledElements;
}

void keepBestLoops(std::vector<std::pair<std::pair<unsigned int, unsigned int>, double>>& loops, unsigned keep) {
    struct Compare {
        bool operator()(const std::pair<std::pair<unsigned int, unsigned int>, double>& a, const std::pair<std::pair<unsigned int, unsigned int>, double>& b) {
            if (std::abs(a.second - b.second) < 0.03) {
                // If qualities are very close, prefer longer loops
                return (a.first.second - a.first.first) > (b.first.second - b.first.first);
            }

            return a.second > b.second; // Compare doubles by their values
        }
    };

    std::priority_queue<std::pair<std::pair<unsigned int, unsigned int>, double>, std::vector<std::pair<std::pair<unsigned int, unsigned int>, double>>, Compare> minHeap;

    for (const auto& entry : loops) {
        if (minHeap.size() < keep) {
            minHeap.push(entry);
        } else if (entry.second > minHeap.top().second) {
            minHeap.pop();
            minHeap.push(entry);
        }
    }

    loops.clear();

    while (!minHeap.empty()) {
        loops.push_back(minHeap.top());
        minHeap.pop();
    }

    // Reverse the vector to maintain the highest values first
    std::reverse(loops.begin(), loops.end());
}
