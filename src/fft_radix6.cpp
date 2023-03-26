#include "fft.hpp"

std::vector<cmplx> fft_radix6(std::vector<cmplx> xin, int N, int opts) {
	cmplx I( { 0.0, 1.0 });
	assert(N % 3 == 0);
	std::vector<cmplx> xout(N);
	std::vector<cmplx> even, odd1, odd2, odd3, odd4;
	for (int n = 0; n < N / 3; n++) {
		even.push_back(xin[3 * n]);
	}
	for (int n = 0; n < N / 6; n++) {
		odd1.push_back(xin[6 * n + 1]);
		odd2.push_back(xin[6 * n + 2]);
		odd3.push_back(xin[6 * n + 4]);
		odd4.push_back(xin[6 * n + 5]);
	}
	even = fft(even, N / 3, opts);
	odd1 = fft(odd1, N / 6, opts);
	odd2 = fft(odd2, N / 6, opts);
	odd3 = fft(odd3, N / 6, opts);
	odd4 = fft(odd4, N / 6, opts);
	for (int k = 0; k < N / 6; k++) {
		auto t0 = twiddle(1 * k, N) * odd1[k];
		auto t1 = twiddle(2 * k, N) * odd2[k];
		auto t2 = twiddle(4 * k, N) * odd3[k];
		auto t3 = twiddle(5 * k, N) * odd4[k];
		auto e0 = even[k + 0 * N / 6];
		auto e1 = even[k + 1 * N / 6];
		cmplx mu0 = twiddle(0, 6);
		cmplx mu1 = twiddle(1, 6);
		cmplx mu2 = twiddle(2, 6);
		cmplx mu3 = twiddle(3, 6);
		cmplx mu4 = twiddle(4, 6);
		cmplx mu5 = twiddle(5, 6);
		/*
		 * 0 0 0 0 0 0
		 * 0 1 2 3 4 5
		 * 0 2 4 0 2 4
		 * 0 3 0 3 0 3
		 * 0 4 2 0 4 2
		 * 0 5 4 3 2 1
		 */
		auto s1 = t0 + t2;
		auto s2 = t1 + t3;
		xout[k + 0 * N / 6] = e0 + t0 + t1 + t2 + t3;
		xout[k + 1 * N / 6] = e1 + mu1 * t0 + mu2 * t1 + mu4 * t2 + mu5 * t3;
		xout[k + 2 * N / 6].x = e0.x + mu2.x * (s1.x + s2.x) + mu2.y * (s2.y - s1.y);
		xout[k + 2 * N / 6].y = e0.y + mu2.x * (s1.y + s2.y) + mu2.y * (s1.x - s2.x);
		xout[k + 3 * N / 6] = e1 + t1 + t2 + mu3 * (t0 + t3);
		xout[k + 4 * N / 6].x = e0.x + mu2.x * (s1.x + s2.x) + mu2.y * (s1.y - s2.y);
		xout[k + 4 * N / 6].y = e0.y + mu2.x * (s1.y + s2.y) - mu2.y * (s1.x - s2.x);
		xout[k + 5 * N / 6] = e1 + mu5 * t0 + mu4 * t1 + mu2 * t2 + mu1 * t3;
	}
	return xout;
}
