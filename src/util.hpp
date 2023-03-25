/*
 * util.hpp
 *
 *  Created on: Mar 5, 2023
 *      Author: dmarce1
 */

#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <vector>
#include <set>
#include <unordered_map>
#include <map>
#include <functional>
#include <vector>
#include <complex>

inline bool close2(double a, double b ) {
	return std::abs(a - b) < 1.0e-12;
}

std::vector<std::vector<int>> nchoosek(int n, int k);
bool is_prime(int n);
int greatest_prime_factor(int N);
std::map<int, int> prime_factorization(int N);
int mod_pow(int a, int b, int m);
int mod_inv(int a, int m);
int generator(int N);
std::vector<int> raders_ginvq(int N);
const std::vector<int> raders_gq(int N);
void fftw(std::vector<std::complex<double>>& x);
const std::vector<std::complex<double>> twiddles(int N);
const std::vector<std::complex<double>> raders_four_twiddle(int N);
const std::vector<std::complex<double>> raders_four_twiddle(int N, int M);
bool are_coprime(int a, int b);

#endif
