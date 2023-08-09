/*************************************************************

  LoopFinder: A spectral hybrid loopfindng algorithm

  Brooklyn Rose Ludlow // August 2023 // MIT License

**************************************************************/

#include <iostream>
#include <filesystem>
#include <LoopFinder.h>
#include <sndfile.hh>

namespace fs = std::filesystem;

struct WAVInfo {
  int format;
  unsigned int sampleRate;
  int channels;
  int minorFormat;
  unsigned int dataSize;
  double* lpcmData;
};

WAVInfo readWAVFile(std::string filePath) {
  WAVInfo info;

  SndfileHandle sfHandle;
  sfHandle = SndfileHandle(filePath, SFM_READ);

  if (!sfHandle) {
    throw std::runtime_error("Error opening the WAV file.");
  }
  if (sfHandle.channels() != 1 || (sfHandle.format() & SF_FORMAT_TYPEMASK) != SF_FORMAT_WAV) {
    throw std::invalid_argument("Unsupported audio format. The file must be a mono (1 channel) WAV file.");
  }

  // Get header data from file
  info.format = sfHandle.format();
  info.sampleRate = sfHandle.samplerate();
  info.channels = sfHandle.channels();
  info.minorFormat = sfHandle.format() & SF_FORMAT_SUBMASK;
  info.dataSize = sfHandle.frames();

  // Read audio data into array
  info.lpcmData = new double[info.dataSize];
  sfHandle.read(info.lpcmData, info.dataSize);

  return info;
}

void CrossFadeAtLoopPoint(double* lpcmData, unsigned dataSize, std::pair<int, int> loopPoints, unsigned fadeLengthSamples) {

  if (fadeLengthSamples > loopPoints.first){
    throw std::invalid_argument("Crossfade length exceeds audio data before loop start.");
  }

  // linear from 0 to 1
  for (unsigned i = 0; i < fadeLengthSamples; i++){
    double fraction = i * 1.0 / (fadeLengthSamples - 1);
    unsigned fadePos = loopPoints.second + i - fadeLengthSamples  + 1;
    unsigned crossfadePos = loopPoints.first + i - fadeLengthSamples + 1;
    lpcmData[fadePos] = (1.0 - fraction) * lpcmData[fadePos] + fraction * lpcmData[crossfadePos];
  }

  // hacky off-by-one prevention
  if (dataSize > loopPoints.second + 1) {
    lpcmData[loopPoints.second + 1] = lpcmData[loopPoints.first + 1];
  }
  if (dataSize > loopPoints.second + 2) {
    lpcmData[loopPoints.second + 2] = lpcmData[loopPoints.first + 2];
  }

}

void SaveLoopPoints(std::string filePath, WAVInfo info, std::pair<int, int> loopPoints) {

  SndfileHandle sfh = SndfileHandle(filePath, SFM_WRITE, info.format, info.channels, info.sampleRate);

  // Set loop
  SF_INSTRUMENT metadata;
  metadata.loop_count = 1;
  metadata.loops[0].mode = SF_LOOP_FORWARD;
  metadata.loops[0].start = loopPoints.first;
  metadata.loops[0].end = loopPoints.second;
  metadata.loops[0].count = 1;

  // Write loop
  int success = sfh.command(SFC_SET_INSTRUMENT, &metadata, sizeof(metadata));

  // Small crossfade to hide boundary clicks
  CrossFadeAtLoopPoint(info.lpcmData, info.dataSize, loopPoints, 512);

  // Write data
  sfh.write(info.lpcmData, info.dataSize);

  if (success == SF_TRUE) {
    std::cout << "Loop written to file: (" << loopPoints.first << ", " << loopPoints.second << ")" << std::endl;
  } else {
    throw std::runtime_error("Error writing loop to file");
  }
}

void processFile(const fs::path& filePath, unsigned sustainStartPos, unsigned minLoopLengthMs) {
  std::cout << "Processing file: " << filePath << std::endl;
  std::string filePathStr = filePath.string();

  // Read WAV file
  WAVInfo info = readWAVFile(filePathStr);

  // Run algorithm
  std::pair<int, int> optimalLoop = LoopFinderAlgorithm(
    info.lpcmData,
    info.dataSize,
    info.sampleRate,
    sustainStartPos,   // sustainStart
    info.dataSize,    // sustainEnd
    minLoopLengthMs  // minLoopMS
  );

  if (optimalLoop.first != 0 || optimalLoop.second != 0) {
    // Write WAV file with new loop point
    SaveLoopPoints(filePathStr, info, optimalLoop);
  }

  // Clean up and free allocated memory
  delete[] info.lpcmData;
}

void processDirectory(const fs::path& dirPath, unsigned sustainStartPos, unsigned minLoopLengthMs) {
  std::cout << "Processing directory: " << dirPath << std::endl;

  for (const auto& entry : fs::directory_iterator(dirPath)) {
    if (fs::is_regular_file(entry)) {
      processFile(entry.path(), sustainStartPos, minLoopLengthMs);
    } else if (fs::is_directory(entry)) {
      processDirectory(entry.path(), sustainStartPos, minLoopLengthMs);
    }
  }
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <path> <arg1> <arg2>" << std::endl;
    return 1;
  }

  fs::path inputPath = argv[1];
  unsigned int sustainStartPos = std::stoul(argv[2]);
  unsigned int minLoopLengthMs = std::stoul(argv[3]);

  if (!fs::exists(inputPath)) {
    std::cerr << "Path does not exist: " << inputPath << std::endl;
    return 1;
  }

  if (fs::is_regular_file(inputPath)) {
    processFile(inputPath, sustainStartPos, minLoopLengthMs);
  } else if (fs::is_directory(inputPath)) {
    processDirectory(inputPath, sustainStartPos, minLoopLengthMs);
  } else {
    std::cerr << "Unsupported path type: " << inputPath << std::endl;
    return 1;
  }
  return 0;
}
