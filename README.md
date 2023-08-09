# LoopFinder


**LoopFinder** is a tool that detects perceptually seamless loop points in digital audio files. It uses a hybrid algorithm of spectral difference analysis and local waveform correlation. It is particularly well suited to complex, semi-periodic waveforms such as a sustained note from a clarinet or a multi-oscillator synth patch.

## Algorithm

The core algorithm of LoopFinder works by breaking audio down into very short pieces, called windows. The FFT power spectrum, which describes how much energy is contained at each frequency, is calculated for each window. The spectra of different windows are compared to identify non-contiguous segments of audio that are spectrally similar.

 > Spectrally similar windows are identified in different parts of the audio data:
 >
 > ..............................**[a][a][a]**..........................**[a][a][a]**..............................

The most similar segments are then aligned at different permutations of their zero-crossings to find the closest correlation of their local waveforms. The results of the spectral and waveform analyses are integrated to create an overall loop quality metric.

See [LoopFinder.cpp](LoopFinder.cpp) for the full explanation and implementation.

## License

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

LoopFinder is open-source software licensed under the GNU General Public License version 3 (GPLv3). The core algorithm and individual source files, apart from `FFT.cpp` and `FFT.h`, are licensed under the MIT license. See the [LICENSE](LICENSE) file for more details.
