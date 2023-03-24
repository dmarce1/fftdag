#pragma once

#include "math.hpp"

#define FFT_REAL 1
#define FFT_DCT2 2
#define FFT_DCT3 4
#define FFT_INV 8
#define FFT_DCT1 16
#define FFT_DCT4 32

std::vector<math_vertex> fft(std::vector<math_vertex> xin, int N, int opts );
std::vector<cmplx> fft_modsplit(std::vector<cmplx> xin, int N, int opts = 0);

