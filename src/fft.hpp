#pragma once

#include "math.hpp"

#define FFT_REAL 1
#define FFT_INV 2
#define FFT_DCT1 4
#define FFT_DCT2 8
#define FFT_DCT3 16
#define FFT_DCT4 32
#define FFT_DST1 64
#define FFT_DST2 128
#define FFT_DST3 256
#define FFT_DST4 512

std::vector<math_vertex> fft(std::vector<math_vertex> xin, int N, int opts );
std::vector<cmplx> fft_modsplit(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_cooley_tukey(int N1, int N2, std::vector<cmplx> xin, int opts);
std::vector<cmplx> fft_radix6(std::vector<cmplx> xin, int N, int opts);
void fft_reset();
void print_fft_bests();
