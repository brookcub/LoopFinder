/*************************************************************

  LoopFinder: A spectral hybrid loopfindng algorithm

  Brooklyn Rose Ludlow // August 2023 // MIT License

**************************************************************/

#include "StatUtils.h"

double BartlettWeightedAverage(double* data, unsigned dataSize) {
    double totalWeightedSum = 0.0;
    double totalWeight = 0.0;

    for (int i = 0; i < dataSize; i++) {
        double weight = 1.0 / (1.0 + std::abs(i - static_cast<int>(dataSize) / 2));
        totalWeightedSum += data[i] * weight;
        totalWeight += weight;
    }

    if (totalWeight == 0.0) {
        return 0.0;  // Handle the case where all weights are zero
    }
    return totalWeightedSum / totalWeight;
}

double MeanAbsoluteDifference(const double* dataA, const double* dataB, int dataSize) {
    double difference = 0;
    for (int i = 0; i < dataSize; i++) {
      difference += fabs(dataA[i] - dataB[i]);
    }
    return difference / dataSize;
}

double RMSE(const double* dataA, const double* dataB, unsigned dataSize) {
    double squaredErrorSum = 0.0;

    // Calculate the squared error for each corresponding pair of elements
    for (unsigned i = 0; i < dataSize; i++) {
        double error = dataA[i] - dataB[i];
        squaredErrorSum += error * error;
    }

    double meanSquaredError = squaredErrorSum / dataSize;
    double rootMeanSquaredError = std::sqrt(meanSquaredError);
    return rootMeanSquaredError;
}

double JensenShannonDivergence(const double* dataA, const double* dataB, unsigned dataSize) {
    double klDivergence1 = 0.0;
    double klDivergence2 = 0.0;
    double epsilon = 1e-7; // Small epsilon value to prevent taking log(0)

    // Calculate the Kullback-Leibler divergences between each data point and their average
    for (unsigned i = 0; i < dataSize; i++) {
        double avg = (dataA[i] + dataB[i]) / 2.0;
        if (dataA[i] > 0.0 && avg > 0.0) {
            klDivergence1 += dataA[i] * std::log2((dataA[i] + epsilon) / (avg + epsilon));
        }
        if (dataB[i] > 0.0 && avg > 0.0) {
            klDivergence2 += dataB[i] * std::log2((dataB[i] + epsilon) / (avg + epsilon));
        }
    }

    return 0.5 * (klDivergence1 + klDivergence2);
}
