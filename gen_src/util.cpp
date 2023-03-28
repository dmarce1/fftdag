#include <vector>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <map>
#include <cmath>
#include <set>
#include <complex>
#include <fftw3.h>

int mod_pow(int a, int b, int m) {
	int rc = 1;
	int apow = a;
	if (a < 0 && b % 2 == 1) {
		return -mod_pow(-a, b, m);
	}
	while (b) {
		if (b & 1) {
			rc = ((rc % m) * (apow % m)) % m;
		}
		b >>= 1;
		apow = ((apow % m) * (apow % m)) % m;
	}
	return rc;
}

static int mod_inv(int a, int m) {
	return mod_pow(a, m - 2, m);
}

static int generator(int N) {
	static thread_local std::unordered_map<int, int> values;
	auto i = values.find(N);
	if (i == values.end()) {
		for (int g = 2;; g++) {
			std::set<int> I;
			bool fail = false;
			for (int m = 0; m < N - 2; m++) {
				int n = mod_pow(g, m, N);
				if (I.find(n) == I.end()) {
					I.insert(n);
				} else {
					fail = true;
					break;
				}
			}
			if (!fail) {
				values[N] = g;
				i = values.find(N);
				break;
			}
		}
	}
	return i->second;
}

std::vector<int> raders_ginvq(int N) {
	const int g = generator(N);
	std::vector<int> ginvq;
	for (int q = 0; q < N - 1; q++) {
		ginvq.push_back(mod_inv(mod_pow(g, q, N), N));
	}
	return ginvq;
}

const std::vector<int> raders_gq(int N) {
	const int g = generator(N);
	std::vector<int> gq;
	for (int q = 0; q < N - 1; q++) {
		gq.push_back(mod_pow(g, q, N));
	}
	return gq;
}

void fftw(std::vector<std::complex<double>>& x) {
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, fftw_complex*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);
		plans[N] = fftw_plan_dft_1d(N, in[N], out[N], FFTW_FORWARD, FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n][0] = x[n].real();
		i[n][1] = x[n].imag();
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		x[n].real(o[n][0]);
		x[n].imag(o[n][1]);
	}

}

const std::vector<std::complex<double>> twiddles(int N) {
	std::vector<std::complex<double>> tw(N);
	for (int k = 0; k < N; k++) {
		tw[k] = std::polar(1.0, -2 * M_PI * k / N);
	}
	return tw;
}

const std::vector<std::complex<double>> raders_four_twiddle(int N) {
	std::vector<std::complex<double>> b(N - 1);
	const auto tws = twiddles(N);
	const auto ginvq = raders_ginvq(N);
	for (int q = 0; q < N - 1; q++) {
		b[q] = tws[ginvq[q]];
	}
	//fftw(b);
	return b;
}

const std::vector<std::complex<double>> raders_four_twiddle(int N, int M) {
	std::vector<std::complex<double>> b(M, 0);
	const auto tws = twiddles(N);
	const auto ginvq = raders_ginvq(N);
	b[0] = tws[ginvq[0]] * (1.0 / M);
	for (int q = 1; q < N - 1; q++) {
		b[M + q - (N - 1)] = b[q] = tws[ginvq[q]];
	}
//	fftw(b);
	return b;
}

__int128 factorial(__int128 k) {
	if (k <= 1) {
		return 1;
	} else {
		return k * factorial(k - 1);
	}
}

std::vector<std::complex<double>> chirp_z_filter(int N) {
	std::vector<std::complex<double>> z(N);
	for( int n = 0; n < N; n++) {
		z[n] = std::polar(1.0, M_PI * n * n / N);
	}
	return z;
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
	while (N != 1) {
		int k = greatest_prime_factor(N);
		if (map.find(k) == map.end()) {
			map[k] = 0;
		}
		map[k]++;
		N /= k;
	}
	return map;
}

bool are_coprime(int a, int b) {
	auto afacs = prime_factorization(a);
	auto bfacs = prime_factorization(b);
	for (auto i : afacs) {
		for (auto j : bfacs) {
			if (i.first == j.first) {
				return false;
			}
		}
	}
	return true;
}
