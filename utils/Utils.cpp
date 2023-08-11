/*************************************************************

  LoopFinder: A spectral hybrid loopfindng algorithm

  Brooklyn Rose Ludlow // August 2023 // MIT License

**************************************************************/

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


std::vector<std::pair<unsigned, unsigned>> determineWindowsToCompare(double* spectralCentroids, unsigned numWindows, unsigned analysisBufferWindows, unsigned centroidComparisonRange) {
    std::vector<std::pair<unsigned, unsigned>> windowsToCompare;

    // Create a map to store centroid buckets and the corresponding window indices
    std::unordered_map<unsigned, std::vector<unsigned>> buckets;

    for (unsigned i = analysisBufferWindows; i < numWindows - analysisBufferWindows; i++) {
        double centroid = spectralCentroids[i];
        int bucketIndex = static_cast<int>(centroid / centroidComparisonRange);
        buckets[bucketIndex].push_back(i);
    }

    // Add pairs of window indices to output vector
    for (auto it = buckets.begin(); it != buckets.end(); ++it) {
        unsigned bucketIndex = it->first;
        const std::vector<unsigned>& windowIndices = it->second;

        // Add all pairs of windows indices within the bucket
        for (unsigned i = 0; i < windowIndices.size(); i++) {
            unsigned indexA = windowIndices[i];
            for (unsigned j = i + 1; j < windowIndices.size(); j++) {
                unsigned indexB = windowIndices[j];
                // smaller index always goes first
                if (indexA < indexB) {
                    windowsToCompare.push_back({indexA, indexB});
                } else {
                    windowsToCompare.push_back({indexB, indexA});
                }
            }
        }

        // Also add pairs between this bucket and the next largest bucket, if it exists
        // This is to ensure that we don't miss any pairs of windows that are close in frequency but cross a bucket boundary
        auto nextBucketIt = buckets.find(bucketIndex + 1);

        if (nextBucketIt != buckets.end()) {
            unsigned nextBucketIndex = nextBucketIt->first;
            const std::vector<unsigned>& nextWindowIndices = nextBucketIt->second;

            for (unsigned i = 0; i < windowIndices.size(); i++) {
                unsigned indexA = windowIndices[i];
                for (unsigned j = 0; j < nextWindowIndices.size(); j++) {
                    unsigned indexB = nextWindowIndices[j];
                    if (std::abs(spectralCentroids[indexA] - spectralCentroids[indexB]) > centroidComparisonRange) {
                        continue;
                    }

                    if (indexA < indexB) {
                        windowsToCompare.push_back({indexA, indexB});
                    } else {
                        windowsToCompare.push_back({indexB, indexA});
                    }
                }
            }
        }
    }

    return windowsToCompare;

}
