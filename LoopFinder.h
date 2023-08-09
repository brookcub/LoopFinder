/*************************************************************

  LoopFinder: A spectral hybrid loopfindng algorithm

  Brooklyn Rose Ludlow // August 2023 // MIT License

**************************************************************/

#ifndef LOOPFINDER_H
#define LOOPFINDER_H

#include <cmath>
#include <utility>
#include <vector>

std::pair<unsigned, unsigned> LoopFinderAlgorithm(
    const double lpcmData[],
    const unsigned dataSize,
    const unsigned sampleRate,
    const unsigned sustainStart,
    const unsigned sustainEnd,
    const unsigned minLoopDurationMs
);

#endif
