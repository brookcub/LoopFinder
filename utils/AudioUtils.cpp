#include "AudioUtils.h"

double* copyLPCMData(const double* lpcmData, const unsigned startPos, const unsigned size) {
    double* window = new double[size];
    for (unsigned i = 0; i < size; i++) {
        window[i] = lpcmData[i + startPos];
    }
    return window;
}

std::vector<unsigned> findZeroCrossings(const double* lpcmData, const unsigned startPos, const unsigned endPos) {
    // Return all zero crossings in a segment of LPCM data
    std::vector<unsigned> zeroCrossings;
    for (unsigned i = startPos; i < endPos - 1; i++) {
        if ((lpcmData[i] <= 0 && lpcmData[i + 1] >= 0) || (lpcmData[i] >= 0 && lpcmData[i + 1] <= 0)) {
            zeroCrossings.push_back(i);
        }
    }
    return zeroCrossings;
}

void curveSpectrum(double* spectrum, const unsigned windowSize) {
    for (unsigned i = 0; i < windowSize / 2; i++) {
        spectrum[i] = 10 * log10(spectrum[i]);
    }
}

double normalizeSpectrumRMSE(double x) {
    // this is sketch but it basically does what I want

    // Ensure x is within the specified range
    double minValue = 6.0;
    double maxValue = 40.0;
    x = std::max(minValue, std::min(maxValue, x));

    // Normalize using a power function
    x = 1.0 - (x - minValue) / (maxValue - minValue);
    x = std::pow(x, 3);

    // Curve with sigmoid
    //double center = 0.7;
    //double steepness = 10.0;
    //x =  (1.0 / (1.0 + exp(-steepness * (x - center + 0.5))) - 0.5) * 2.0;

    return x;
}

double normalizeSpectrumJSDivergence(double x) {
    return normalizeSpectrumRMSE(x / 10.0 + 6);
}

double normalizeSampleWindowMAD(double x) {
    const double minValue = 0.0;
    const double maxValue = 1.0;

    x *= 9;
    x = std::pow(x, 1.5) * 1.5;
    x = std::max(minValue, std::min(maxValue, 1.0 - x));

    return x;
}
