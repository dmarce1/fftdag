#include "dag.hpp"

#include <cmath>
#include <complex>

struct cmplx {
	dag_node x;
	dag_node y;
	cmplx operator+(cmplx other) {
		x += other.x;
		y += other.y;
		return *this;
	}
	cmplx operator*(cmplx other) {
		auto tmpx = x;
		auto tmpy = y;
		x = tmpx * other.x - tmpy * other.y;
		y = tmpy * other.x + tmpx * other.y;
		return *this;
	}
};

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
	assert(dag_node::sanity_check());
	for (int k2 = 0; k2 < N2; k2++) {
		for (int k1 = 0; k1 < N1; k1++) {
			cmplx x;
			cmplx t;
			cmplx w;
			x.x = (sub[0][2 * k2]);
			x.y = (sub[0][2 * k2 + 1]);
			for (int n1 = 1; n1 < N1; n1++) {
				const double theta = -2.0 * M_PI * (n1 * (N2 * k1 + k2)) / N;
				t.x = (sub[n1][2 * k2]);
				t.y = (sub[n1][2 * k2 + 1]);
				w.x = (cos(theta));
				w.y = (sin(theta));
				x = x + w * t;
			}
			xout[2 * (k2 + k1 * N2)] = x.x;
			xout[2 * (k2 + k1 * N2) + 1] = x.y;
		}
	}
	return xout;
}
