#include <complex>
#include <vector>
#include <unordered_map>
#include <memory>
#include <chrono>
#include <numeric>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include "sfft.hpp"
#include "util.hpp"
#include "types.hpp"

void fft_complex_simd4(fft_simd4* X, int N) {
	sfft_complex((double*) X, ((double*) X) + sizeof(fft_simd4), 2 * sizeof(fft_simd4), 2 * sizeof(fft_simd4), N);
}

void * operator new(std::size_t n) {
	void* memptr;
	posix_memalign(&memptr, 32, n);
	return memptr;
}

void operator delete(void * p) {
	free(p);
}

void *operator new[](std::size_t s) {
	void* memptr;
	posix_memalign(&memptr, 32, s);
	return memptr;
}

void operator delete[](void *p) {
	free(p);
}

#define SIMD_SIZE 4

int round_down(int i, int m) {
	return m * (i / m);
}

const std::vector<complex<double>> scalar_twiddles(int N) {
	using entry_type = std::shared_ptr<std::vector<complex<double>>>;
	static std::unordered_map<int, entry_type> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return *(iter->second);
	} else {
		std::vector<complex<double>> tw(N);
		for (int k = 0; k < N; k++) {
			tw[k].real() = cos(-2 * M_PI * k / N);
			tw[k].imag() = sin(-2 * M_PI * k / N);
		}
		cache[N] = std::make_shared<std::vector<complex<double>>>(std::move(tw));
		return *(cache[N]);
	}
}

const std::vector<std::vector<complex<fft_simd4>>>& vector_twiddles(int N1, int N2) {
	using entry_type = std::shared_ptr<std::vector<std::vector<complex<fft_simd4>>>>;
	static std::unordered_map<int, std::unordered_map<int, entry_type>> cache;
	auto iter = cache[N1].find(N2);
	if (iter != cache[N1].end()) {
		return *(iter->second);
	} else {
		std::vector<std::vector<complex<fft_simd4>>> W(N1 / SIMD_SIZE, std::vector<complex<fft_simd4>>(N2));
		for (int k2 = 0; k2 < N2; k2++) {
			int n1;
			for (n1 = 0; n1 < round_down(N1, SIMD_SIZE); n1 ++) {
				W[n1/SIMD_SIZE][k2].real()[n1 % SIMD_SIZE] = cos(-2.0 * M_PI * n1 * k2 / (N1*N2));
				W[n1/SIMD_SIZE][k2].imag()[n1 % SIMD_SIZE] = sin(-2.0 * M_PI * n1 * k2 / (N1*N2));
			}
		}
		cache[N1][N2] = std::make_shared<std::vector<std::vector<complex<fft_simd4>>>>(std::move(W));
		return *(cache[N1][N2]);
	}
}

class timer {
	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
	double time;
public:
	inline timer() {
		time = 0.0;
	}
	inline void stop() {
		std::chrono::time_point<std::chrono::high_resolution_clock> stop_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> dur = stop_time - start_time;
		time += dur.count();
	}
	inline void start() {
		start_time = std::chrono::high_resolution_clock::now();
	}
	inline void reset() {
		time = 0.0;
	}
	inline double read() {
		return time;
	}
};

bool can_cooley_tukey(int N) {
	static std::unordered_map<int, bool> cache;
	auto iter = cache.find(N);
	if (iter == cache.end()) {
		auto factors = prime_factorization(N);
		cache[N] = (factors.rbegin()->first <= FFT_NMAX);
		return cache[N];
	}
	return iter->second;
}

int cooley_tukey_radix(int N) {
	static std::unordered_map<int, int> cache;
	auto iter = cache.find(N);
	if (iter != cache.end()) {
		return iter->second;
	} else {
		int N1 = sqrt(N);
		while (N % N1 != 0) {
			N1--;
		}
		cache[N] = N1;
		return N1;
	}
}

void cooley_tukey(complex<fft_simd4>* X, complex<fft_simd4>* Y, int N) {
	if (N <= FFT_NMAX) {
		fft_complex_simd4((fft_simd4*) X, N);
		return;
	}
	static std::vector<complex<fft_simd4>> z;
	const auto& W = scalar_twiddles(N);
	const int N1 = cooley_tukey_radix(N);
	const int N2 = N / N1;
	z.resize(N1);
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			Y[N2 * n1 + n2] = X[N1 * n2 + n1];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		cooley_tukey(Y + n1 * N2, X + n1 * N2, N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		z[0] = Y[k2];
		for (int n1 = 1; n1 < N1; n1++) {
			z[n1] = Y[N2 * n1 + k2] * W[n1 * k2];
		}
		fft_complex_simd4((fft_simd4*) (z.data()), N1);
		for (int k1 = 0; k1 < N1; k1++) {
			X[N2 * k1 + k2] = z[k1];
		}
	}
}

/*void cooley_tukey_s(complex<double>* X, complex<double>* Y, int N) {
 if (N <= FFT_NMAX) {
 fft_complex((double*) X, N);
 return;
 }
 static std::vector<complex<double>> z;
 const auto& W = scalar_twiddles(N);
 const int N1 = cooley_tukey_radix(N);
 const int N2 = N / N1;
 z.resize(N1);
 for (int n1 = 0; n1 < N1; n1++) {
 for (int n2 = 0; n2 < N2; n2++) {
 Y[N2 * n1 + n2] = X[N1 * n2 + n1];
 }
 }
 for (int n1 = 0; n1 < N1; n1++) {
 cooley_tukey_s(Y + n1 * N2, X + n1 * N2, N2);
 }
 for (int k2 = 0; k2 < N2; k2++) {
 z[0] = Y[k2];
 for (int n1 = 1; n1 < N1; n1++) {
 z[n1] = Y[N2 * n1 + k2] * W[n1 * k2];
 }
 fft_complex((double*) (z.data()), N1);
 for (int k1 = 0; k1 < N1; k1++) {
 X[N2 * k1 + k2] = z[k1];
 }
 }
 }
 */
/*
 void cooley_tukey2(complex<double>* X, int N) {
 if (N <= FFT_NMAX) {
 fft_complex((double*) X, N);
 return;
 }
 const auto& W = scalar_twiddles(N);
 const int N1 = cooley_tukey_radix(N);
 const int N2 = N / N1;
 const int N1v = round_down(N1, SIMD_SIZE);
 const int N2v = round_down(N2, SIMD_SIZE);
 std::vector<std::vector<complex<fft_simd4>>>Yv(N1v/SIMD_SIZE, std::vector<complex<fft_simd4>>(N2));
 std::vector<std::vector<complex<fft_simd4>>>Zv(N2v/SIMD_SIZE, std::vector<complex<fft_simd4>>(N1));
 std::vector<std::vector<complex<fft_simd4>>>Y0v(N1v/SIMD_SIZE, std::vector<complex<fft_simd4>>(N2));
 std::vector<std::vector<complex<fft_simd4>>>Z0v(N2v/SIMD_SIZE, std::vector<complex<fft_simd4>>(N1));
 std::vector<std::vector<complex<double>>>Ys(N1 - N1v, std::vector<complex<double>>(N2));
 std::vector<std::vector<complex<double>>>Zs(N2 - N2v, std::vector<complex<double>>(N1));
 int n1;
 const auto& Wv = vector_twiddles(N1, N2);
 const auto& Ws = scalar_twiddles(N);

 for (n1 = 0; n1 < N1v; n1++) {
 for (int n2 = 0; n2 < N2; n2++) {
 Yv[n1 / SIMD_SIZE][n2].real()[n1 % SIMD_SIZE] = X[N1 * n2 + n1].real();
 Yv[n1 / SIMD_SIZE][n2].imag()[n1 % SIMD_SIZE] = X[N1 * n2 + n1].imag();
 }
 }
 for (; n1 < N1; n1++) {
 for (int n2 = 0; n2 < N2; n2++) {
 Ys[n1 - N1v][n2] = X[N1 * n2 + n1];
 }
 }
 for (n1 = 0; n1 < N1v; n1 += SIMD_SIZE) {
 cooley_tukey(Yv[n1 / SIMD_SIZE].data(), Y0v[n1 / SIMD_SIZE].data(), N2);
 }
 for (; n1 < N1; n1++) {
 cooley_tukey2(Ys[n1 - N1v].data(), N2);
 }
 for (n1 = 0; n1 < N1v; n1 += SIMD_SIZE) {
 for (int k2 = 0; k2 < N2; k2++) {
 Yv[n1 / SIMD_SIZE][k2] *= Wv[n1 / SIMD_SIZE][k2];
 }
 }
 for (; n1 < N1; n1++) {
 for (int k2 = 0; k2 < N2; k2++) {
 Ys[n1 - N1v][k2] *= Ws[n1 * k2];
 }
 }
 for (n1 = 0; n1 < N1v; n1++) {
 for (int n2 = 0; n2 < N2; n2++) {
 X[N1 * n2 + n1].real() = Yv[n1 / SIMD_SIZE][n2].real()[n1 % SIMD_SIZE];
 X[N1 * n2 + n1].imag() = Yv[n1 / SIMD_SIZE][n2].imag()[n1 % SIMD_SIZE];
 }
 }
 for (; n1 < N1; n1++) {
 for (int n2 = 0; n2 < N2; n2++) {
 X[N1 * n2 + n1] = Ys[n1 - N1v][n2];
 }
 }

 int k2;
 for (k2 = 0; k2 < N2v; k2++) {
 for (n1 = 0; n1 < N1; n1++) {
 Zv[k2 / SIMD_SIZE][n1].real()[k2 % SIMD_SIZE] = X[N1 * k2 + n1].real();
 Zv[k2 / SIMD_SIZE][n1].imag()[k2 % SIMD_SIZE] = X[N1 * k2 + n1].imag();
 }
 }
 for (; k2 < N2; k2++) {
 for (int n1 = 0; n1 < N1; n1++) {
 Zs[k2 - N2v][n1] = X[N1 * k2 + n1];
 }
 }
 for (k2 = 0; k2 < N2v; k2 += SIMD_SIZE) {
 cooley_tukey(Zv[k2 / SIMD_SIZE].data(), Z0v[k2 / SIMD_SIZE].data(), N1);
 }
 for (; k2 < N2; k2++) {
 cooley_tukey2(Zs[k2 - N2v].data(), N1);
 }
 for (k2 = 0; k2 < N2v; k2++) {
 for (int k1 = 0; k1 < N1; k1++) {
 X[N2 * k1 + k2].real() = Zv[k2 / SIMD_SIZE][k1].real()[k2 % SIMD_SIZE];
 X[N2 * k1 + k2].imag() = Zv[k2 / SIMD_SIZE][k1].imag()[k2 % SIMD_SIZE];
 }
 }
 for (; k2 < N2; k2++) {
 for (int k1 = 0; k1 < N1; k1++) {
 X[N2 * k1 + k2] = Zs[k2 - N2v][k1];
 }
 }

 }
 */
timer tm1, tm2, tm3, tm4;

void FFT(std::vector<complex<double>>& Z) {
	int N = Z.size();
	std::vector<double> Y(4 * N);
	std::vector<double> X(4 * N);
	for (int i = 0; i < 4 * N; i++) {
		X[i] = Z[i / 4].real();
		Y[i] = Z[i / 4].imag();
	}
	tm2.start();
	sfft_complex(X.data(), Y.data(), 4, 4, N);
	tm2.stop();
	for (int i = 0; i < 4 * N; i++) {
		Z[i / 4].real() = X[i];
		Z[i / 4].imag() = Y[i];
	}
}

double rand1() {
	return (rand() + 0.5) / (RAND_MAX + 1.0);
}

int main(int argc, char **argv) {
	for (int N = 2; N <= FFT_NMAX; N++) {
		while (!can_cooley_tukey(N)) {
			N++;
		}
		double avg_err = 0.0;
		for (int i = 0; i < 100; i++) {
			std::vector<complex<double>> X(N);
			std::vector<std::complex<double>> Y(N);
			for (int n = 0; n < N; n++) {
				Y[n].real(X[n].real() = rand1());
				Y[n].imag(X[n].imag() = rand1());
				//Y[n].real(X[n].real() = rand1());
//				Y[n].imag(X[n].imag() = rand1());
			}
			if (i != 0) {
				tm1.start();
				tm3.start();
			}
			fftw(Y);
			if (i != 0) {
				tm1.stop();
				tm3.stop();
				tm4.start();
			}
			FFT(X);
			if (i != 0) {
				tm4.stop();
			}
			for (int n = 0; n < N; n++) {
				double x = X[n].real() - Y[n].real();
				double y = X[n].imag() - Y[n].imag();
				double err = sqrt(x * x + y * y);
				avg_err += err;
		   //	printf("%e %e | %e %e | %e\n", X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag(), err);
			}
			//abort();
		}
		avg_err /= (255 * N);
		auto pfac = prime_factorization(N);
		std::string f;
		for (auto i = pfac.begin(); i != pfac.end(); i++) {
			f += "(" + std::to_string(i->first) + "^" + std::to_string(i->second) + ")";
		}
		printf("%i: %32s | %e %e %e %e %e\n", N, f.c_str(), avg_err, tm1.read(), tm2.read(), tm2.read() / tm1.read(), tm4.read() / tm3.read());
		tm2.reset();
		tm1.reset();
	}
	return 0;
}
