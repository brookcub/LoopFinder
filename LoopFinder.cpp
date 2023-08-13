/*************************************************************

  LoopFinder: A spectral hybrid loopfindng algorithm

  Brooklyn Rose Ludlow // August 2023 // MIT License

**************************************************************/

#include "LoopFinder.h"
#include "AudioUtils.h"
#include "StatUtils.h"
#include "Utils.h"
#include "FFT.h"

#include <iostream>

//     _   _              _ _   _
//    /_\ | |__ _ ___ _ _(_) |_| |_  _ __
//   / _ \| / _` / _ \ '_| |  _| ' \| '  \
//  /_/ \_\_\__, \___/_| |_|\__|_||_|_|_|_|
//          |___/

// Meaningful Parameters by Section
//
// 1 ) Compute Power Spectra
//     - FFTWindowSize
//     - FFTWindowType
//     - FFTWindowOverlap
// 2 ) Spectral Similarity Calculation
//     - windowComparisonRange
//     - spectralSimilarityThreshold
// 3 ) Local Sample Correlation Calculation
//     - spectralWindowsToCheck
//     - localSampleRange
// 4 ) Metric Integration
//     - spectralWeight
//     - sampleWeight
//     - qualityThreshold

// Naming conventions:
// "Pos" refers specifically to sample position indices in the raw lpcmData
// "Index" refers to any other indices (zeroCrossing, windows, spectra, etc)
// This avoids some confusion when we're keeping track of many indices

std::pair<unsigned, unsigned> LoopFinderAlgorithm(
  const double lpcmData[],
  const unsigned lpcmDataSize,
  const unsigned sampleRate,
  const unsigned sustainStartPos,
  const unsigned sustainEndPos,
  const unsigned minLoopDurationMs
) {

  /////////////////////////////////////////////////////////////////////////////
  //  INTERNAL PARAMETERS AND SETUP
  /////////////////////////////////////////////////////////////////////////////
  const std::pair<unsigned, unsigned> nullReturn = std::make_pair(0, 0);

  // Finding candidates and subsampling
  const unsigned consecutiveCandidates = 15;
  const unsigned rejectWindowMs = 20;

  // FFT
  const unsigned FFTWindowSize = 2048; // 1024 should be good. could try 2048 w overlaps.any bigger would obscure the features we're interested in at high pitches
  const unsigned FFTWindowType = 3; // 3 = Hanning. 3 is great, 4, 6, and 9 look same as 3
  const unsigned FFTWindowOverlap = FFTWindowSize / 2;

  // Spectral Similarity
  const unsigned windowComparisonRange = 8;
  const double spectralSimilarityThreshold = 0.5;

  // Local Sample Correlation
  const unsigned spectralWindowsToCheck = 50;
  const unsigned localSampleRange = 800;

  // Metric Integration
  const double spectralWeight = 0.5;
  const double sampleWeight = 0.5;
  const double qualityThreshold = 0.85;

  // General useful things
  const unsigned FFTWindowSizeMinusOverlap = FFTWindowSize - FFTWindowOverlap;
  const unsigned analysisBuffer = std::max(FFTWindowSizeMinusOverlap * windowComparisonRange, localSampleRange);
  const unsigned analysisBufferWindows = std::max(windowComparisonRange, localSampleRange / FFTWindowSizeMinusOverlap);
  const unsigned earliestLoopablePos = sustainStartPos + analysisBuffer;
  const unsigned latestLoopablePos = sustainEndPos - analysisBuffer;

  if ((sustainEndPos - sustainStartPos) < FFTWindowSize) {
    throw std::invalid_argument("Sustain range is smaller than window size");
  }

  /////////////////////////////////////////////////////////////////////////////
  //  WINDOW AND FFT AUDIO DATA
  /////////////////////////////////////////////////////////////////////////////

  const unsigned numWindows = (sustainEndPos - sustainStartPos) / FFTWindowSizeMinusOverlap - 1;

  // Allocate memory for windows, powerSpectra, and spectralCentroids
  double** windows = new double*[numWindows];
  double** powerSpectra = new double*[numWindows];
  double* spectralCentroids = new double[numWindows];

  // Calculate power spectra FFT for each window of audio
  for (unsigned i = 0; i < numWindows; i++) {
      unsigned windowStartPos = i * FFTWindowSizeMinusOverlap + sustainStartPos;
      windows[i] = copyLPCMData(lpcmData, windowStartPos, FFTWindowSize);

      // Apply window function to the lpcm window data
      WindowFunc(FFTWindowType, FFTWindowSize, windows[i]);

      // Perform the FFT
      powerSpectra[i] = new double[FFTWindowSize / 2];
      PowerSpectrum(FFTWindowSize, windows[i], powerSpectra[i]);

      // Apply perceptual curving to the PowerSpectra
      curveSpectrum(powerSpectra[i], FFTWindowSize);

      // Compute spectral centroid
      spectralCentroids[i] = computeSpectralCentroid(powerSpectra[i], sampleRate, FFTWindowSize);
  }

  /////////////////////////////////////////////////////////////////////////////
  //  FIND SPECTRALLY SIMILAR WINDOWS
  /////////////////////////////////////////////////////////////////////////////

  const unsigned numComparisonWindows = windowComparisonRange * 2 + 1;
  const unsigned minLoopDurationSamples = sampleRate * minLoopDurationMs / 1000;
  const unsigned minLoopDurationWindows = 1 +  (minLoopDurationSamples - 1) / (FFTWindowSizeMinusOverlap);
  const unsigned FFTCenterOffset = FFTWindowSize / 4;
  const unsigned centroidComparisonRange = 250; // todo: make this a function of spectralsimilaritythreshold. 250 is very very roughly 0.85 normalized RMSE

  // reduce our comparison pool by grouping windows by spectral centroid in O(n) time
  std::vector<std::pair<unsigned, unsigned>> loopWindows = determineWindowsToCompare(
    spectralCentroids,
    numWindows,
    analysisBufferWindows,
    centroidComparisonRange
  );

  std::vector<std::pair<std::pair<unsigned, unsigned>, double >> spectralLoopCandidates;

  // memoize duplicate spectral similarity calculations
  double* memoizationTable = new double[numWindows * numWindows];
  std::fill_n(memoizationTable, numWindows * numWindows, -1.0);

  // Loop over pairs, storing pairs that are within the spectral similarity threshold
  for (unsigned i = 0; i < loopWindows.size(); i++) {
    unsigned startWindowIndex = loopWindows[i].first;
    unsigned endWindowIndex = loopWindows[i].second;

    // skip windows that are too close to meet minimum loop length
    if (endWindowIndex - startWindowIndex < minLoopDurationWindows) {
      continue;
    }

    // Calculate spectral similarity for windows in comparison range
    double* spectralResults = new double[numComparisonWindows];

    for (unsigned w = 0; w < numComparisonWindows; w++) {
      unsigned windowAIndex = startWindowIndex + w - windowComparisonRange;
      unsigned windowBIndex = endWindowIndex + w - windowComparisonRange;

      // Check if the similarity value is already calculated and stored
      if (memoizationTable[windowAIndex * numWindows + windowBIndex] != -1.0) {
        spectralResults[w] = memoizationTable[windowAIndex * numWindows + windowBIndex];
      } else {
        spectralResults[w] = RMSE(powerSpectra[windowAIndex], powerSpectra[windowBIndex], FFTWindowSize / 2);
        memoizationTable[windowAIndex * numWindows + windowBIndex] = spectralResults[w];
      }
    }

    // Average and normalize the results
    double averageRMSE = BartlettWeightedAverage(spectralResults, numComparisonWindows);
    double normalRMSE = normalizeSpectrumRMSE(averageRMSE);

    delete[] spectralResults;

    // Add good loops to candidate vector
    if (normalRMSE >= spectralSimilarityThreshold) {
      spectralLoopCandidates.push_back(
        std::make_pair(
          loopWindows[i],
          normalRMSE
        )
      );
    }
  }

  // Free memory for windows and powerSpectra, we're done with them
  for (unsigned i = 0; i < numWindows; i++) {
      delete[] windows[i];
      delete[] powerSpectra[i];
  }

  delete[] windows;
  delete[] powerSpectra;
  delete[] memoizationTable;
  delete[] spectralCentroids;

  // Now we have a vector spectralLoopCandidates with the structure:
  // ( (startWindowIndex, endWindowIndex), spectralSimilarity )

  if (spectralLoopCandidates.empty()) {
    return nullReturn;
    throw std::logic_error("No regions meet spectral criteria. Try reducing minLoopDurationMs or spectralSimilarityThreshold.");
  }

  /////////////////////////////////////////////////////////////////////////////
  //  LOCAL SAMPLE CORRELATION
  /////////////////////////////////////////////////////////////////////////////

  unsigned rejectWindow = sampleRate * rejectWindowMs / 1000;

  const unsigned localWindowSize = localSampleRange * 2 + 1;
  std::vector<std::pair<std::pair<unsigned, unsigned>, double>> finalLoops;

  // loop over top X window pairs, perform sample-per-sample correlation, integrate the results
  keepBestLoops(spectralLoopCandidates, spectralWindowsToCheck); // filter for best candidates

  for (unsigned i = 0; i < spectralLoopCandidates.size(); i++) {
    unsigned startWindowIndex = spectralLoopCandidates[i].first.first;
    unsigned endWindowIndex = spectralLoopCandidates[i].first.second;
    double spectralSimilarity = spectralLoopCandidates[i].second;

    unsigned startWindowStartPos = (startWindowIndex + 1) * FFTWindowSizeMinusOverlap + sustainStartPos - FFTCenterOffset;
    unsigned startWindowEndPos = startWindowStartPos + FFTWindowSizeMinusOverlap - 1;

    unsigned endWindowStartPos = (endWindowIndex + 1) * FFTWindowSizeMinusOverlap + sustainStartPos - FFTCenterOffset;
    unsigned endWindowEndPos = endWindowStartPos + FFTWindowSizeMinusOverlap - 1;

    //unsigned startWindowCenterPos = FFTWindowSizeMinusOverlap * (startWindowIndex + 1) + sustainStartPos;
    //unsigned endWindowCenterPos = FFTWindowSizeMinusOverlap * (endWindowIndex + 1) + sustainStartPos;

    // Get zero crossings for local window analysis
    std::vector<unsigned> startZeroCrossings = findZeroCrossings(lpcmData, startWindowStartPos, startWindowEndPos);
    std::vector<unsigned> endZeroCrossings  = findZeroCrossings(lpcmData, endWindowStartPos, endWindowEndPos); // todo - it would be better to center this on the window

    // subsample the zero crossings to reduce the number of comparisons
    std::vector<unsigned> subsampledStartZeroCrossings = strideSubsample(startZeroCrossings, 1, FFTWindowSizeMinusOverlap / 3); // spread out across window
    std::vector<unsigned> subsampledEndZeroCrossings = strideSubsample(endZeroCrossings, consecutiveCandidates, rejectWindow); // clustered

    // iterate over each pair of zero crossings, and find the best match
    for (unsigned j = 0; j < subsampledStartZeroCrossings.size(); j++) {
      unsigned startZeroCrossing = startZeroCrossings[j];
      double* startPointSampleWindow = copyLPCMData(lpcmData, startZeroCrossing - localSampleRange, localWindowSize);
      WindowFunc(2, localWindowSize, startPointSampleWindow); // 1, 2, or 6 work well
      for (unsigned k = 0; k < subsampledEndZeroCrossings.size(); k++) {
        unsigned endZeroCrossing = endZeroCrossings[k];

        // skip if the loop is too short
        if (endZeroCrossing - startZeroCrossing < minLoopDurationSamples) {
          continue;
        }

        // Compute mean absolute difference of the start and end windows
        double* endPointSampleWindow = copyLPCMData(lpcmData, endZeroCrossing - localSampleRange, localWindowSize);
        WindowFunc(2, localWindowSize, endPointSampleWindow);

        double mad = MeanAbsoluteDifference(startPointSampleWindow, endPointSampleWindow, localWindowSize);
        double sampleSimilarity = normalizeSampleWindowMAD(mad);
        delete[] endPointSampleWindow;

        // Integrate the spectral and sample results
        double integratedQuality = (spectralWeight * spectralSimilarity + sampleWeight * sampleSimilarity) / (spectralWeight + sampleWeight);

        // Add good loops to final vector
        if (integratedQuality >= qualityThreshold) {
          finalLoops.push_back(
            std::make_pair(
              std::make_pair(
                startZeroCrossing,
                endZeroCrossing
              ),
              integratedQuality
            )
          );
        }
      }

      delete[] startPointSampleWindow;
    }
  }

  if (!finalLoops.empty()) {
    keepBestLoops(finalLoops, 1); // get best loop
    std::cout << "The top loop candidate is: " << finalLoops[0].first.first << ", " << finalLoops[0].first.second << " with quality: " << finalLoops[0].second << " and length: " << 1.0 * (finalLoops[0].first.second - finalLoops[0].first.first) / sampleRate << std::endl;
    return std::make_pair(finalLoops[0].first.first, finalLoops[0].first.second);
  }

  return nullReturn;
}
