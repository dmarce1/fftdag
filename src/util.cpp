#include <vector>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <map>
#include <cmath>

__int128 factorial(__int128 k) {
	if (k <= 1) {
		return 1;
	} else {
		return k * factorial(k - 1);
	}
}

std::vector<std::vector<int>> nchoosek(int n, int k) {
	std::vector<std::vector<int>> rc;
	std::vector<int> combo(k);
	std::iota(combo.begin(), combo.end(), 0);
	if (n == k) {
		rc.push_back(combo);
	} else {
		bool done = false;
		while (!done) {
			rc.push_back(combo);
			int dim = k - 1;
			while (combo[dim] == dim + n - k) {
				dim--;
				if (dim < 0) {
					done = true;
					break;
				}
			}
			if (!done) {
				combo[dim]++;
				for (int i = dim + 1; i < k; i++) {
					combo[i] = combo[i - 1] + 1;
				}
			}
		}
	}
	rc.push_back(combo);
	return rc;
}

bool is_prime(int n) {
	static thread_local std::unordered_map<int, bool> values;
	auto i = values.find(n);
	if (i == values.end()) {
		bool v = true;
		if (n == 1) {
			v = false;
		} else {
			int kmax = sqrt(n);
			for (auto i = 2; i <= kmax; i++) {
				if (n % i == 0) {
					v = false;
					break;
				}
				if (i > kmax) {
					break;
				}
			}
		}
		values[n] = v;
		i = values.find(n);
	}
	return i->second;
}

int greatest_prime_factor(int N) {
	static thread_local std::unordered_map<int, int> values;
	auto i = values.find(N);
	if (i == values.end()) {
		int v;
		for (int n = 2; n <= N; n++) {
			if (N % n == 0 && is_prime(n)) {
				v = n;
			}
		}
		values[N] = v;
		i = values.find(N);
	}
	return i->second;
}

std::map<int, int> prime_factorization(int N) {
	std::map<int, int> map;
	while( N != 1 ) {
		int k = greatest_prime_factor(N);
		if( map.find(k) == map.end()) {
			map[k] = 0;
		}
		map[k]++;
		N /= k;
	}
	return map;
}
