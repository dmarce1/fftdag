#include "dag.hpp"


#include <cmath>
#include <complex>

std::vector<dag_node> fft_prime_power(int R, std::vector<dag_node> xin, int N) {
	if (N == 1) {
		return xin;
	}
	const int N1 = R;
	const int N2 = N / R;
	std::vector<dag_node> xout(2 * N);
	std::vector<std::vector<dag_node>> sub(N1, std::vector<dag_node>(2 * N2));

	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			sub[n1][2 * n2] = xin[2 * (n2 * N1 + n1)];
			sub[n1][2 * n2 + 1] = xin[2 * (n2 * N1 + n1) + 1];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		sub[n1] = fft_prime_power(N1, sub[n1], N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int k1 = 0; k1 < N1; k1++) {
			std::complex<dag_node> x, t, w;
			x.real(sub[0][2 * k2]);
			x.imag(sub[0][2 * k2 + 1]);
			for (int n1 = 1; n1 < N1; n1++) {
				const double theta = -2.0 * M_PI * (n1 * (N2 * k1 + k2)) / N;
				t.real(sub[n1][2 * k2]);
				t.imag(sub[n1][2 * k2 + 1]);
				w.real(cos(theta));
				w.imag(sin(theta));
				x = x + w * t;
			}
			xout[2 * (k2 + k1 * N2)] = x.real();
			xout[2 * (k2 + k1 * N2) + 1] = x.imag();
		}
	}
	return xout;
}
