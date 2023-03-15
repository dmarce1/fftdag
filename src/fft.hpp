#pragma once

#include "math.hpp"

std::vector<math_vertex> fft_radix2(std::vector<math_vertex> xin, int N);
std::vector<math_vertex> fft_prime_power(int R, std::vector<math_vertex> xin, int N);
std::vector<math_vertex> fft_radix4(std::vector<math_vertex> xin, int N);
std::vector<math_vertex> fft_singleton(std::vector<math_vertex> xin, int N);
std::vector<math_vertex> fft_prime_factor(int N1, int N2, int R1, int R2, std::vector<math_vertex> xin);
