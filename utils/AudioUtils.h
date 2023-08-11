/*************************************************************

  LoopFinder: A spectral hybrid loopfindng algorithm

  Brooklyn Rose Ludlow // August 2023 // MIT License

**************************************************************/

#ifndef AUDIOUTILS_H
#define AUDIOUTILS_H

#include <cmath>
#include <vector>

double* copyLPCMData(const double* lpcmData, const unsigned startPos, const unsigned size);

std::vector<unsigned> findZeroCrossings(const double* lpcmData, const unsigned startPos, const unsigned endPos);

void curveSpectrum(double* spectrum, const unsigned windowSize);

double normalizeSpectrumRMSE(double x);
double normalizeSpectrumJSDivergence(double x);

double normalizeSampleWindowMAD(double x);

double computeSpectralCentroid(const double* powerSpectrum, unsigned sampleRate, unsigned windowSize);

#endif // AUDIOUTILS_H
