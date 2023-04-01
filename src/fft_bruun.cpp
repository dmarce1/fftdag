#include "fft.hpp"
#include "polynomial.hpp"

using con_poly = polynomial<std::complex<double>>;
using var_poly = polynomial<cmplx>;

std::vector<cmplx> fft_bruun(const var_poly& X0, const con_poly& F0, const std::vector<con_poly>& W0, int N) {
	std::vector<cmplx> Y(N);
/*	int R = least_prime_factor(N);
	if (N == 1) {
		return std::vector<cmplx>(1, (X0 % F0)[0]);
	}
	std::vector<con_poly> F1(R, F0);
	std::vector<con_poly> Wp(R);
	std::vector<std::vector<con_poly>> W1(R);
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
			Y[I[r][n]] = Y1[n];
		}
	}*/
	return Y;
}

std::vector<cmplx> fft_bruun(std::vector<cmplx> xin, int N, int opts) {
	std::vector<cmplx> xout;
	con_poly D;
	var_poly X;
	std::vector<con_poly> W;
	for (int n = 0; n < N; n++) {
		con_poly w;
		w[1] = 1.0;
		w[0] = std::polar(1.0, -2.0 * M_PI * mod(n, N) / N);
		W.push_back(w);
		auto z = std::polar(1.0, -M_PI * n);
		X[n] = xin[n] * cmplx({z.real(), z.imag()});
	}
	D[N] = 1.0;
	D[0] = -1.0;
	return fft_bruun(X, D, W, N);
}
