#include "fft.hpp"
#include "convolve.hpp"

std::vector<cmplx> fft_raders(std::vector<cmplx> xin, int N, bool padded, int opts) {
	assert(is_prime(N));
	std::vector<cmplx> xout(N);
	int M;
	if (padded) {
		M = 1;
		while (M < 2 * (N - 1)) {
			M *= 2;
		}
	} else {
		M = N - 1;
	}
	const auto& b = padded ? raders_four_twiddle(N, M) : raders_four_twiddle(N);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	std::vector<cmplx> a(M, cmplx( { 0.0, 0.0 }));
	std::vector<cmplx> c(M, cmplx( { 0.0, 0.0 }));
	for (int q = 0; q < N - 1; q++) {
		a[q] = xin[gq[q]];
	}
	xout[0] = cmplx( { 0.0, 0.0 });
	int o = 0;
	if (M != N - 1) {
		o = M << 16;
	}
	xout[0] = xin[0];
	for (int n = 1; n < N; n++) {
		xout[0] += xin[n];
	}
	c = convolve(a, b, opts);
	for (int p = 0; p < N - 1; p++) {
		xout[ginvq[p]] = xin[0] + c[p];
	}
	return std::move(xout);
}
