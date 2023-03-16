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

inline bool close2(double a, double b ) {
	return std::abs(a - b) < 1.0e-12;
}

std::vector<std::vector<int>> nchoosek(int n, int k);
bool is_prime(int n);
int greatest_prime_factor(int N);
std::map<int, int> prime_factorization(int N);

#endif
