#include "fft.hpp"
#include "convolve.hpp"

int calc_padding(int N) {
	int M = 1;
	while (M < 2 * (N - 1)) {
		M *= 2;
	}
	return M;
}

std::vector<cmplx> fft_raders2(std::vector<cmplx> xin, int N, bool padded, int opts) {
	std::vector<cmplx> xout(N);
	const auto odd_index = [N]( int n1, int n2 ) {
		int i = mod_pow(3, n2, N);
		if( n1 > 0 ) {
			i = N - i;
		}
		return i;
	};
	auto prime_fac = prime_factorization(N);
	constexpr int P = 2;
	const int c = prime_fac.begin()->second;
	std::vector<std::vector<std::complex<double>>> tw(2, std::vector<std::complex<double>>(N / 4));
	std::vector<cmplx> x0(N / 2);
	std::vector<std::vector<cmplx>> x1(2, std::vector<cmplx>(N / 4));
	std::vector<std::vector<cmplx>> x2(N / 4, std::vector<cmplx>(2));
	for (int k1 = 0; k1 < 2; k1++) {
		for (int k2 = 0; k2 < N / 4; k2++) {
			int k = odd_index(k1, k2);
			tw[k1][k2] = std::polar(1.0, -2.0 * M_PI * k / N);
		}
	}
	for (int n = 0; n < N / 2; n++) {
		x0[n] = xin[2 * n];
	}
	for (int n1 = 0; n1 < 2; n1++) {
		for (int n2 = 0; n2 < N / 4; n2++) {
			int n = odd_index(n1, n2);
			x1[n1][n2] = xin[n];
		}
	}
	x0 = fft(x0, N / 2, opts);
	x1 = convolve(x1, tw, opts);

	return xout;
}

std::vector<cmplx> fft_raders(std::vector<cmplx> xin, int N, bool padded, int opts) {
	std::vector<cmplx> xout(N, cmplx( { 0.0, 0.0 }));
	auto prime_fac = prime_factorization(N);
	assert(prime_fac.size() == 1);
	int P = prime_fac.begin()->first;
	int c = prime_fac.begin()->second;
	const auto& b = raders_four_twiddle(N);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	int L = std::pow(P, c - 1);
	int M = L * (P - 1);
	std::vector<cmplx> x0(L, cmplx( { 0.0, 0.0 }));
	std::vector<cmplx> x1(L, cmplx( { 0.0, 0.0 }));
	std::vector<cmplx> x2(M, cmplx( { 0.0, 0.0 }));
	for (int n2 = 0; n2 < L; n2++) {
		for (int n1 = 0; n1 < P; n1++) {
			x0[n2] += xin[L * n1 + n2];
		}
	}
	for (int n1 = 0; n1 < L; n1++) {
		x1[n1] = xin[P * n1];
	}
	for (int q = 0; q < M; q++) {
		x2[q] = xin[gq[q]];
	}
	x0 = fft(x0, L, opts);
	x1 = fft(x1, L, opts);
	x2 = convolve(x2, b, opts);
	for (int k1 = 0; k1 < L; k1++) {
		xout[P * k1] += x0[k1];
	}
	if (L > 1) {
		for (int k1 = 0; k1 < P; k1++) {
			for (int k2 = 0; k2 < L; k2++) {
				if (k2 % P != 0) {
					xout[L * k1 + k2] += x1[k2];
				}
			}
		}
	} else {
		for (int k2 = 1; k2 < P; k2++) {
			xout[k2] = x1[0];
		}
	}
	for (int p = 0; p < M; p++) {
		xout[ginvq[p]] += x2[p];
	}
	return std::move(xout);
}

