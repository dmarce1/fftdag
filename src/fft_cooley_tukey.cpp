#include "fft.hpp"

std::vector<cmplx> fft_cooley_tukey(int N1, int N2, std::vector<cmplx> xin, int opts) {
	int N = N1 * N2;
	std::vector<cmplx> xout(N);
	std::vector<std::vector<cmplx>> X1(N1, std::vector<cmplx>(N2));
	std::vector<std::vector<cmplx>> X2(N2, std::vector<cmplx>(N1));
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			X1[n1][n2] = xin[(n1 + n2 * N1) % N];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		X1[n1] = fft(X1[n1], N2, opts);
	}
	for (int n1 = 0; n1 < N1; n1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			auto tw = twiddle(n1 * k2, N);
			X2[k2][n1] = X1[n1][k2] * tw;
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		X2[k2] = fft(X2[k2], N1, opts);
	}
	for (int k1 = 0; k1 < N1; k1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			xout[k2 + k1 * N2] = X2[k2][k1];
		}
	}
	return xout;

}
std::vector<cmplx> fft_cooley_tukey2(int N1, int N2, std::vector<cmplx> xin, int opts) {
	int N = N1 * N2;
	std::vector<cmplx> xout(N);
	std::vector<std::vector<cmplx>> X1(N1, std::vector<cmplx>(N2));
	std::vector<std::vector<cmplx>> X2(N2, std::vector<cmplx>(N1));
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			X1[n1][n2] = xin[(n1 + n2 * N1) % N];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		X1[n1] = fft(X1[n1], N2, opts);
	}
	for (int n1 = 0; n1 < N1; n1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			auto tw = twiddle(n1 * k2, N);
			X2[k2][n1] = X1[n1][k2] * tw;
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		X2[k2] = fft(X2[k2], N1, opts);
	}
	for (int k1 = 0; k1 < N1; k1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			xout[k2 + k1 * N2] = X2[k2][k1];
		}
	}
	return xout;

}
