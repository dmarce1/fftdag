#include "fft.hpp"
#include "convolve.hpp"

int calc_padding(int N) {
	int M = 1;
	while (M < 2 * (N - 1)) {
		M *= 2;
	}
	return M;
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

