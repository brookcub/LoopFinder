/*************************************************************

  LoopFinder: A spectral hybrid loopfindng algorithm

  Brooklyn Rose Ludlow // August 2023 // MIT License

**************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <queue>
#include <unordered_map>
#include <vector>

std::vector<unsigned> strideSubsample(const std::vector<unsigned>& v, const unsigned consecutiveElements, const unsigned rejectWindow);

void keepBestLoops(std::vector<std::pair<std::pair<unsigned int, unsigned int>, double>>& loops, unsigned keep);

std::vector<std::pair<unsigned, unsigned>> determineWindowsToCompare(double* spectralCentroids, unsigned numWindows, unsigned analysisBufferWindows, unsigned centroidComparisonRange);

#endif // UTILS_H
