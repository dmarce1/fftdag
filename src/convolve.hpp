#pragma once
#include <array>

std::vector<cmplx> convolve(std::vector<cmplx> x, std::vector<std::complex<double>> h, int opts);
std::vector<cmplx> convolve_fast(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> convolve_fft(std::vector<cmplx> x, std::vector<std::complex<double>> h);

