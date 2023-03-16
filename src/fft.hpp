#pragma once

#include "math.hpp"

#define FFT_REAL 1
#define FFT_SYM 2

std::vector<math_vertex> fft_prime_power(int R, std::vector<math_vertex> xin, int N);
std::vector<math_vertex> fft_radix4(std::vector<math_vertex> xin, int N);
std::vector<math_vertex> fft_singleton(std::vector<math_vertex> xin, int N);
std::vector<math_vertex> fft_prime_factor(int N1, int N2, std::vector<math_vertex> xin);
std::vector<math_vertex> fft(std::vector<math_vertex> xin, int N, int opts = 0);
std::vector<math_vertex> fft_raders(std::vector<math_vertex> xin, int N);
