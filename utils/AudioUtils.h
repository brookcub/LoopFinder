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

#endif // AUDIOUTILS_H
