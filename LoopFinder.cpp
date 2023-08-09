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
// 1 ) Initial zero-crossing candidate selection
//     - N/A
// 2 ) Subsampling of candidates
//     - consecutiveCandidates
//     - rejectWindowMs
// 3 ) Compute Power Spectra
//     - FFTWindowSize
//     - FFTWindowType
// 4 ) Spectral Similarity Calculation
//     - windowComparisonRange
//     - spectralSimilarityThreshold
// 5 ) Local Sample Correlation Calculation
//     - spectralCandidatesToCheck
//     - localSampleRange
// 6 ) Metric Integration
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

  // Spectral Similarity
  const unsigned windowComparisonRange = 8;
  const double spectralSimilarityThreshold = 0.5;

  // Local Sample Correlation
  const unsigned spectralCandidatesToCheck = 6000;
  const unsigned localSampleRange = 800;

  // Metric Integration
  const double spectralWeight = 0.5;
  const double sampleWeight = 0.5;
  const double qualityThreshold = 0.85;

  // General useful things
  const unsigned analysisBuffer = std::max(FFTWindowSize * windowComparisonRange, localSampleRange);
  const unsigned earliestLoopablePos = sustainStartPos + analysisBuffer;
  const unsigned latestLoopablePos = sustainEndPos - analysisBuffer;

  if ((sustainEndPos - sustainStartPos) < FFTWindowSize) {
    throw std::invalid_argument("Sustain range is smaller than window size");
  }

  /////////////////////////////////////////////////////////////////////////////
  //  FINDING CANDIDATES
  /////////////////////////////////////////////////////////////////////////////

  std::vector<unsigned> zeroCrossings  = findZeroCrossings(lpcmData, earliestLoopablePos, latestLoopablePos);

  if (zeroCrossings.size() < 2) {
    throw std::invalid_argument("Sustain region must contain at least 2 zero crossings to loop.");
  }

  // zeroCrossings now contains all zero crossings in the sustain section in ascending order

  /////////////////////////////////////////////////////////////////////////////
  //  SUBSAMPLE CANDIDATES (for performance)
  /////////////////////////////////////////////////////////////////////////////

  // Subsample with this algorithm
  //  1) Keep X consecutive candidates
  //  2) Then reject all candidates for Y milliseconds afterwards
  //  - Repeat until all candidates have been sampled
  // This gives us clusters of points across the audio.
  // . . . .        . . . .        . . . .        . . . .

  unsigned rejectWindow = sampleRate * rejectWindowMs / 1000;
  std::vector<unsigned> subsampledZeroCrossings = strideSubsample(zeroCrossings, consecutiveCandidates, rejectWindow);

  zeroCrossings.clear();

  /////////////////////////////////////////////////////////////////////////////
  //  WINDOW AND FFT AUDIO DATA
  /////////////////////////////////////////////////////////////////////////////

  const unsigned numWindows = (sustainEndPos - sustainStartPos - FFTWindowSize) / (FFTWindowSize / 2) + 1;

  // Allocate memory for windows and powerSpectra
  double** windows = new double*[numWindows];
  double** powerSpectra = new double*[numWindows];

  // Calculate power spectra FFT for each window of audio
  for (unsigned i = 0; i < numWindows; i++) {
      unsigned windowStartPos = i * (FFTWindowSize / 2) + sustainStartPos;
      windows[i] = copyLPCMData(lpcmData, windowStartPos, FFTWindowSize);

      // Apply window function to the lpcm window data
      WindowFunc(FFTWindowType, FFTWindowSize, windows[i]);

      // Perform the FFT
      powerSpectra[i] = new double[FFTWindowSize / 2];
      PowerSpectrum(FFTWindowSize, windows[i], powerSpectra[i]);

      // Apply perceptual curving to the PowerSpectra
      curveSpectrum(powerSpectra[i], FFTWindowSize);
  }

  /////////////////////////////////////////////////////////////////////////////
  //  COMPUTE SPECTRAL SIMILARITY
  /////////////////////////////////////////////////////////////////////////////

  const unsigned numComparisonWindows = windowComparisonRange * 2 + 1;
  const unsigned minLoopDurationSamples = sampleRate * minLoopDurationMs / 1000;
  const unsigned FFTCenterOffset = FFTWindowSize / 4;

  std::vector<std::pair<std::pair<unsigned, unsigned>, double > > spectralLoopCandidates;

  double* memoizationTable = new double[numWindows * numWindows];
  std::fill_n(memoizationTable, numWindows * numWindows, -1.0);

  // loop over each pair of subsampledZeroCrossings
  for (unsigned i = 0; i < subsampledZeroCrossings.size() - 1; i++) {
    unsigned loopStartPos = subsampledZeroCrossings[i];
    unsigned startWindowIndex = 2 * (loopStartPos - sustainStartPos + FFTCenterOffset) / FFTWindowSize - 1;

    for (unsigned j = subsampledZeroCrossings.size() - 1; j > i + 1; j--) {
      unsigned loopEndPos = subsampledZeroCrossings[j];

      if (loopEndPos - loopStartPos < minLoopDurationSamples) {
        continue; // skip if the loop is too short
      }

      unsigned endWindowIndex = 2 * (loopEndPos - sustainStartPos + FFTCenterOffset) / FFTWindowSize - 1;

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
            std::make_pair((loopStartPos), (loopEndPos)),
            normalRMSE
          )
        );
      }

    } // end candidate loop
  } // start candidate loop

  // Free memory for windows and powerSpectra, we're done with them
  for (unsigned i = 0; i < numWindows; i++) {
      delete[] windows[i];
      delete[] powerSpectra[i];
  }

  delete[] windows;
  delete[] powerSpectra;
  delete[] memoizationTable;

  // We're also done with subsampledZeroCrossings
  subsampledZeroCrossings.clear();

  // Now we have a vector spectralLoopCandidates with the structure:
  // ( (loopStartPos, loopEndPos), spectralSimilarity )

  if (spectralLoopCandidates.empty()) {
    return nullReturn;
    throw std::logic_error("No loop candidates meet spectral criteria. Try reducing minLoopDurationMs or spectralSimilarityThreshold.");
  }

  std::cout << "There are " << spectralLoopCandidates.size() << " spectralLoopCandidates" << std::endl;

  /////////////////////////////////////////////////////////////////////////////
  //  LOCAL SAMPLE CORRELATION
  /////////////////////////////////////////////////////////////////////////////

  const unsigned localWindowSize = localSampleRange * 2 + 1;
  std::vector<std::pair<std::pair<unsigned, unsigned>, double>> finalLoops;

  // loop over top X loops, perform sample-per-sample correlation, integrate the results
  keepBestLoops(spectralLoopCandidates, spectralCandidatesToCheck); // filter for best candidates

  for (unsigned i = 0; i < spectralLoopCandidates.size(); i++) {
    unsigned loopStartPos = spectralLoopCandidates[i].first.first;
    unsigned loopEndPos = spectralLoopCandidates[i].first.second;

    // Compute mean absolute difference of the start and end windows
    double* startPointSampleWindow = copyLPCMData(lpcmData, loopStartPos - localSampleRange, localWindowSize);
    double* endPointSampleWindow = copyLPCMData(lpcmData, loopEndPos - localSampleRange, localWindowSize);
    WindowFunc(2, localWindowSize, startPointSampleWindow); // 1, 2, or 6 work well
    WindowFunc(2, localWindowSize, endPointSampleWindow);

    double mad = MeanAbsoluteDifference(startPointSampleWindow, endPointSampleWindow, localWindowSize);
    double sampleSimilarity = normalizeSampleWindowMAD(mad);

    delete[] startPointSampleWindow;
    delete[] endPointSampleWindow;

    // Integrate the spectral and sample results
    double spectralSimilarity = spectralLoopCandidates[i].second;
    double integratedQuality = (spectralWeight * spectralSimilarity + sampleWeight * sampleSimilarity) / (spectralWeight + sampleWeight);

    // Add good loops to final vector
    if (integratedQuality >= qualityThreshold) {
      finalLoops.push_back(
        std::make_pair(
          std::make_pair(
            (loopStartPos),
            (loopEndPos)
          ),
          integratedQuality
        )
      );
    }
  }

  if (!finalLoops.empty()) {
    keepBestLoops(finalLoops, 1); // get best loop
    std::cout << "The top loop candidate is: " << finalLoops[0].first.first << ", " << finalLoops[0].first.second << " with quality: " << finalLoops[0].second << " and length: " << 1.0 * (finalLoops[0].first.second - finalLoops[0].first.first) / sampleRate << std::endl;
    return std::make_pair(finalLoops[0].first.first, finalLoops[0].first.second);
  }

  return nullReturn;
}
