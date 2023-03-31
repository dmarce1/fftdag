#include "fft.hpp"
#include "polynomial.hpp"

using complex_polynomial = polynomial<std::complex<long double>>;

std::vector<std::complex<long double>> fft_bruun(const complex_polynomial& X0, const complex_polynomial& F0, const std::vector<complex_polynomial>& W0, int N, bool root = false) {
	int R = least_prime_factor(N);
	if (N == 1) {
		printf("R = 1  N = %i D = %s\n", N, F0.to_str().c_str());
		return std::vector<std::complex<long double>>(1, (X0 % F0)[0]);
	} else {
		printf("R = %i N = %i D = %s\n", R, N, F0.to_str().c_str());
	}
	std::vector<std::vector<std::complex<long double>>>Y1(R);
	std::vector<complex_polynomial> F1(R, F0);
	std::vector<complex_polynomial> X1(R, F0);
	std::vector<complex_polynomial> Wp(R);
	std::vector<std::vector<complex_polynomial>> W1(R);
	std::vector<std::complex<long double>> Y0(N);
	std::vector<std::vector<int>> I(R);
	if (N == R || N % 4 != 0 || (F0[0].real() < 0.0 && close2(F0[0].imag(), 0.0))) {
		for (int n = 0; n < N / R; n++) {
			for (int r = 0; r < R; r++) {
				int nn = R * n + r;
				I[r].push_back(nn);
				W1[r].push_back(W0[nn]);
			}
		}
	} else {
		for (int n = 0; n < N / R; n++) {
			for (int r = 0; r < R; r++) {
				int r0 = n % 2 == 0 ? r : R - r - 1;
				int nn = R * n + r0;
				I[r].push_back(nn);
				W1[r].push_back(W0[nn]);
			}
		}
	}
	for (int r = 0; r < R; r++) {
		auto Wp = W0[I[r][0]];
		for (int i = 1; i < I[r].size(); i++) {
			Wp *= W0[I[r][i]];
		}
		F1[r] *= Wp;
		F1[r] /= F0;
	}
	for (int r = 0; r < R; r++) {
		auto X1 = X0 % F1[r];
		auto Y1 = fft_bruun(X1, F1[r], W1[r], I[r].size());
		for (int n = 0; n < I[r].size(); n++) {
			Y0[I[r][n]] = Y1[n];
		}
	}
	return Y0;
}

std::vector<std::complex<double>> fft_bruun(std::vector<std::complex<double>> xin, int N) {
	std::vector<std::complex<double>> xout;
	complex_polynomial X;
	complex_polynomial D;
	std::vector<complex_polynomial> W;
	for (int n = 0; n < N; n++) {
		complex_polynomial w;
		w[1] = 1.0;
		w[0] = std::polar(1.0l, -2.0l * M_PIl * mod(n, N) / N);
		W.push_back(w);
		X[n].real(xin[n].real());
		X[n].imag(xin[n].imag());
		X[n] *= std::polar(1.0l, -M_PIl * n);
	}
	D[N] = 1.0;
	D[0] = -1.0;
	auto Y = fft_bruun(X, D, W, N, true);
	X = X % D;
	for (int n = 0; n < N; n++) {
		std::complex<double> y(Y[n].real(), Y[n].imag());
		xout.push_back(y);
	}
	return xout;
}
