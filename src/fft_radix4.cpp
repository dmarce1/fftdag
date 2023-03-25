#include "fft.hpp"

std::vector<cmplx> fft_radix4(std::vector<cmplx> xin, int N, int opts) {
	cmplx I( { 0.0, 1.0 });
	if (N == 1) {
		return xin;
	} else if (N == 2) {
		std::vector<cmplx> xout(2);
		auto x0 = xin[0] + xin[1];
		auto x1 = xin[0] - xin[1];
		xout[0] = x0;
		xout[1] = x1;
		return xout;
	}
	assert(N % 2 == 0);
	std::vector<cmplx> xout(N);
	std::vector<cmplx> even, odd1, odd3;
	for (int n = 0; n < N / 2; n++) {
		even.push_back(xin[2 * n]);
	}
	for (int n = 0; n < N / 4; n++) {
		odd1.push_back(xin[4 * n + 1]);
		odd3.push_back(xin[(4 * n - 1 + N) % N]);
	}
	even = fft(even, N / 2, opts);
	odd1 = fft(odd1, N / 4, opts);
	odd3 = fft(odd3, N / 4, opts);
	for (int k = 0; k < N / 4; k++) {
		auto t1 = twiddle(+k, N) * odd1[k];
		auto t3 = twiddle(-k, N) * odd3[k];
		auto zs = t1 + t3;
		auto zd = t1 - t3;
		auto ur0 = even[k + 0 / 4] + 1 * zs;
		auto ur1 = even[k + N / 4] - I * zd;
		auto ur2 = even[k + 0 / 4] - 1 * zs;
		auto ur3 = even[k + N / 4] + I * zd;
		xout[k + 0 * N / 4] = ur0;
		xout[k + 1 * N / 4] = ur1;
		xout[k + 2 * N / 4] = ur2;
		xout[k + 3 * N / 4] = ur3;
	}
	return xout;
}
