/*************************************************************

  LoopFinder: A spectral hybrid loopfindng algorithm

  Brooklyn Rose Ludlow // August 2023 // MIT License

**************************************************************/

#ifndef STATUTILS_H
#define STATUTILS_H

#include <cmath>

double BartlettWeightedAverage(double* data, unsigned dataSize);
double RMSE(const double* dataA, const double* dataB, unsigned dataSize);
double JensenShannonDivergence(const double* dataA, const double* dataB, unsigned dataSize);
double MeanAbsoluteDifference(const double* dataA, const double* dataB, int dataSize);

#endif // STATUTILS_H
