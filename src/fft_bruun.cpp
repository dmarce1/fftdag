#include "fft.hpp"
#include "polynomial.hpp"

using complex_polynomial = polynomial<std::complex<double>>;

std::vector<std::complex<double>> fft_bruun(const complex_polynomial& X, const complex_polynomial& D, const std::vector<complex_polynomial>& W, int N) {
//	printf("D = %s\n", D.to_str().c_str());
	if (N == 1) {
		return std::vector<std::complex<double>>(1, X[0]);
	} else {
		std::vector<std::complex<double>> Y;
		complex_polynomial D1(N / 2);
		complex_polynomial D2(N / 2);
		std::vector<complex_polynomial> W1;
		std::vector<complex_polynomial> W2;
		if (N > 2) {
			D1[N / 2] = D2[N / 2] = 1.0;
			if (D[0].real() < 0.0) {
				D1[0] = -1.0;
				D2[0] = 1.0;
			} else {
				D1[0] = D2[0] = 1.0;
				D1[N / 4] = sqrt(2.0 - D[N / 2]);
				D2[N / 4] = -D1[N / 4];
			}
		} else {
			D1 = D / W[1];
			D2 = D / W[0];
		}
		for (int n = 0; n < N; n++) {
			auto w = W[n];
			if (close2(std::abs((D1 % w)[0]), 0.0)) {
				W1.push_back(w);
				printf("e ");
			} else {
				W2.push_back(w);
				printf("o ");
			}
		}
		printf("\n");
	//	printf("D1 = %s  D2 = %s\n", D1.to_str().c_str(), D2.to_str().c_str());
		complex_polynomial X1 = X % D1;
		complex_polynomial X2 = X % D2;
		auto Y1 = fft_bruun(X1, D1, W1, N / 2);
		auto Y2 = fft_bruun(X2, D2, W2, N / 2);
		for (int n = 0; n < N / 2; n++) {
			Y.push_back(Y1[n]);
			Y.push_back(Y2[n]);
		}
		return Y;
	}
}

std::vector<std::complex<double>> fft_bruun(std::vector<std::complex<double>> xin, int N) {
	std::vector<std::complex<double>> xout;
	complex_polynomial X(N);
	complex_polynomial D(N);
	std::vector<complex_polynomial> W;
	for (int n = 0; n < N; n++) {
		complex_polynomial w(1);
		w[1] = 1.0;
		w[0] = std::polar(1.0, -2.0 * M_PI * mod(n - N / 2, N) / N);
		W.push_back(w);
		X[n] = xin[n];
	}
	D[N] = 1.0;
	D[0] = -1.0;
	auto Y = fft_bruun(X, D, W, N);
	for (int n = 0; n < N; n++) {
		xout.push_back(Y[n]);
	}
	return xout;
}
