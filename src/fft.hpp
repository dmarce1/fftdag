#pragma once

#include "math.hpp"

std::vector<math_vertex> fft_radix2(std::vector<math_vertex> xin, int N);
std::vector<math_vertex> fft_prime_power(int R, std::vector<math_vertex> xin, int N);
