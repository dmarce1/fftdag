#include "fft.hpp"

std::vector<cmplx> fft_prime_factor(int N1, int N2, std::vector<cmplx> xin, int opts) {
	assert(are_coprime(N1, N2));
	int N = N1 * N2;
	std::vector<cmplx> xout(N);
	std::vector<std::vector<cmplx>> sub1(N1, std::vector<cmplx>(N2));
	std::vector<std::vector<cmplx>> sub2(N2, std::vector<cmplx>(N1));
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			sub1[n1][n2] = xin[(n1 * N2 + n2 * N1) % N];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		sub1[n1] = fft(sub1[n1], N2, opts);
	}
	for (int n1 = 0; n1 < N1; n1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			sub2[k2][n1] = sub1[n1][k2];
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		sub2[k2] = fft(sub2[k2], N1, opts);
	}
	for (int k = 0; k < N; k++) {
		const int k1 = k % N1;
		const int k2 = k % N2;
		xout[k] = sub2[k2][k1];
	}
	return xout;
}
