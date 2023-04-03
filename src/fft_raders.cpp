#include "fft.hpp"
#include "convolve.hpp"

int calc_padding(int N) {
	int M = 2 * N - 1;
	while (1) {
		auto factors = prime_factorization(M);
		bool flag = true;
		int twocnt = 0;
		for (auto i = factors.begin(); i != factors.end(); i++) {
			if (i->first == 2) {
				twocnt = i->second;
				break;
			}
		}
		for (auto i = factors.begin(); i != factors.end(); i++) {
			if ((i->first != 2 && i->first != 3 && i->first != 5) || i->second > twocnt) {
				flag = false;
				break;
			}
		}
		if (flag) {
			break;
		}
		M++;
	}
	return M;
}

std::vector<cmplx> fft_raders_fft(std::vector<cmplx> xin, int N, bool padded, int opts) {
//	padded = true;
	std::vector<cmplx> xout(N, cmplx( { 0.0, 0.0 }));
	auto prime_fac = prime_factorization(N);
	assert(prime_fac.size() == 1);
	int P = prime_fac.begin()->first;
	int c = prime_fac.begin()->second;
	int L = std::pow(P, c - 1);
	int M = L * (P - 1);
	int K = padded ? calc_padding(M) : M;
	auto b = raders_twiddle(N, K, padded);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	std::vector<cmplx> x2(K, cmplx( { 0.0, 0.0 }));
	if (L > 1) {
		std::vector<cmplx> x0(L, cmplx( { 0.0, 0.0 }));
		std::vector<cmplx> x1(L, cmplx( { 0.0, 0.0 }));
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
		x2 = convolve_fft(x2, b);
		for (int k1 = 0; k1 < L; k1++) {
			xout[P * k1] += x0[k1];
		}
		for (int k1 = 0; k1 < P; k1++) {
			for (int k2 = 0; k2 < L; k2++) {
				if (k2 % P != 0) {
					xout[L * k1 + k2] += x1[k2];
				}
			}
		}
	} else {
		xout[0] = xin[0];
		for (int q = 0; q < M; q++) {
			x2[q] = xin[gq[q]];
		}
		x2 = fft(x2, K, opts);
		xout[0] += x2[0];
		fftw(b);
		for (int n = 0; n < K; n++) {
			x2[n] *= b[n] * (1.0 / K);
		}
		for (int n = 0; n < K; n++) {
			std::swap(x2[n].x, x2[n].y);
		}
		x2 = fft(x2, K, opts);
		for (int n = 0; n < M; n++) {
			std::swap(x2[n].x, x2[n].y);
		}
		for (int p = 1; p < P; p++) {
			xout[p] += xin[0];
		}
	}
	for (int p = 0; p < M; p++) {
		xout[ginvq[p]] += x2[p];
	}
	return std::move(xout);
}

std::vector<cmplx> fft_raders_fast(std::vector<cmplx> xin, int N, int opts) {
	std::vector<cmplx> xout(N, cmplx( { 0.0, 0.0 }));
	auto prime_fac = prime_factorization(N);
	assert(prime_fac.size() == 1);
	int P = prime_fac.begin()->first;
	int c = prime_fac.begin()->second;
	int L = std::pow(P, c - 1);
	int M = L * (P - 1);
	const auto& b = raders_twiddle(N, M, false);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
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
	x2 = convolve_fast(x2, b);
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

std::vector<cmplx> fft_raders(std::vector<cmplx> xin, int N, bool padded, int opts) {
	if (padded || N > 64) {
		return fft_raders_fft(xin, N, padded, opts);
	} else {
		int fast_cnt, fft_cnt;
		fft_cnt = math_vertex::operation_count(fft_raders_fft(xin, N, false, opts)).total();
		fast_cnt = math_vertex::operation_count(fft_raders_fast(xin, N, opts)).total();
		if (fft_cnt < fast_cnt) {
			return fft_raders_fft(xin, N, false, opts);
		} else {
			return fft_raders_fast(xin, N, opts);
		}
	}
}
