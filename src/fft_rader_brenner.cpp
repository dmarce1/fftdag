#include "fft.hpp"


std::vector<cmplx> fft_rader_brenner(std::vector<cmplx> xin, int N, int opts) {
	std::vector<cmplx> even(N / 2);
	std::vector<cmplx> odd(N / 2);
	std::vector<cmplx> xout(N);
	cmplx x0 = (xin[0] - xin[N / 2]) - cmplx( { 0.0, 1.0 }) * (xin[N / 4] - xin[3 * N / 4]);
	cmplx x1 = (xin[0] - xin[N / 2]) + cmplx( { 0.0, 1.0 }) * (xin[N / 4] - xin[3 * N / 4]);
	for (int n = 0; n < N / 2; n++) {
		even[n] = xin[n] + xin[n + N / 2];
	}
	odd[0] = 0.0;
	odd[N / 4] = 0.0;
	for (int n = 0; n < N / 2; n++) {
		if (n % (N / 4) == 0) {
			continue;
		}
		odd[n] = (xin[n] - xin[n + N / 2]) * (1.0 / (2.0 * cos(2.0 * M_PI * n / N)));
	}
	even = fft(even, N / 2, opts);
	odd = fft(odd, N / 2, opts);
	for (int k = 0; k < N / 2; k++) {
		xout[2 * k] = even[k];
		xout[2 * k + 1] = odd[k] + odd[(k + 1) % (N / 2)];
		if (k % 2 == 0) {
			xout[2 * k + 1] += x0;
		} else {
			xout[2 * k + 1] += x1;
		}
	}
	return xout;
}
