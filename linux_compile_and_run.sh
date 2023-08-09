g++ main.cpp LoopFinder.cpp utils/*.cpp FFT/FFT.cpp -o loopfinder -I. -Iutils -IFFT -lsndfile -std=c++17

./loopfinder
