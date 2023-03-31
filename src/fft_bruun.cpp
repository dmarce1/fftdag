#include "fft.hpp"
#include "polynomial.hpp"

using complex_polynomial = polynomial<std::complex<double>>;

std::vector<std::complex<double>> fft_bruun(const complex_polynomial& X, const complex_polynomial& D, const std::vector<complex_polynomial>& W, int N, bool root = false) {
	printf("D = %s\n", D.to_str().c_str());
	if (N == 1) {
		return std::vector<std::complex<double>>(1, X[0]);
	}
	std::vector<std::complex<double>> Y;
	complex_polynomial D1;
	complex_polynomial D2;
	std::vector<complex_polynomial> W1;
	std::vector<complex_polynomial> W2;
	std::vector<int> I1;
	std::vector<int> I2;
	if (N > 2) {
		D1 = D;
		D2 = D;
		for (int n = 0; n < N / 2; n++) {
			int n0 = 2 * n + 0;
			int n1 = 2 * n + 1;
			D1 = D1 / (W[n1]);
			D2 = D2 / (W[n0]);
			W1.push_back(W[n0]);
			W2.push_back(W[n1]);
			I1.push_back(n0);
			I2.push_back(n1);
		}
		complex_polynomial X1 = X % D1;
		complex_polynomial X2 = X % D2;
		auto Y1 = fft_bruun(X1, D1, W1, N / 2);
		auto Y2 = fft_bruun(X2, D2, W2, N / 2);
		Y.resize(N);
		for (int n = 0; n < I1.size(); n++) {
			Y[I1[n]] = Y1[n];
		}
		for (int n = 0; n < I2.size(); n++) {
			Y[I2[n]] = Y2[n];
		}
	} else {
		D1 = D / W[0];
		D2 = D / W[1];
		complex_polynomial X1 = X % D1;
		complex_polynomial X2 = X % D2;
		auto Y1 = fft_bruun(X1, D1, W1, 1);
		auto Y2 = fft_bruun(X2, D2, W2, 1);
		Y.resize(N);
		Y[0] = Y1[0];
		Y[1] = Y2[0];
	}

	printf("\n");
	return Y;
}

std::vector<std::complex<double>> fft_bruun(std::vector<std::complex<double>> xin, int N) {
	std::vector<std::complex<double>> xout;
	complex_polynomial X;
	complex_polynomial D;
	std::vector<complex_polynomial> W;
	for (int n = 0; n < N; n++) {
		complex_polynomial w;
		w[1] = 1.0;
		w[0] = std::polar(1.0, -2.0 * M_PI * mod(n, N) / N);
		W.push_back(w);
		X[n] = xin[n];// * std::polar(1.0, -M_PI * n);
	}
	D[N] = 1.0;
	D[0] = -1.0;
	auto Y = fft_bruun(X, D, W, N, true);
	X = X % D;
	for (int n = 0; n < N; n++) {
		xout.push_back(Y[n]);
	}
	return xout;
}
