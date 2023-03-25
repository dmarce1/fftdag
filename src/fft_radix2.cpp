#include "fft.hpp"

std::vector<cmplx> fft_radix2(std::vector<cmplx> xin, int N, int opts) {
	if (N == 1) {
		return xin;
	}
	assert(N % 2 == 0);
	std::vector<cmplx> xout(N);
	std::vector<cmplx> even, odd;
	for (int n = 0; n < N / 2; n++) {
		even.push_back(xin[2 * n]);
	}
	for (int n = 0; n < N / 2; n++) {
		odd.push_back(xin[2 * n + 1]);
	}
	even = fft(even, N / 2, opts);
	odd = fft(odd, N / 2, opts);
	for (int k = 0; k < N / 2; k++) {
		double theta = -2.0 * M_PI * k / N;
		auto tw = twiddle(k, N);
		auto t = odd[k] * tw;
		xout[k] = even[k] + t;
		xout[k + N / 2] = even[k] - t;
	}
	return xout;
}
