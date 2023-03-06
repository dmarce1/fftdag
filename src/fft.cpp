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


std::vector<dag_node> fft_radix2(std::vector<dag_node> xin, int N) {
	if (N == 1) {
		return xin;
	}
	std::vector<dag_node> xout(2 * N);
	std::vector<dag_node> even, odd;
	for (int n = 0; n < N / 2; n++) {
		even.push_back(xin[4 * n]);
		even.push_back(xin[4 * n + 1]);
	}
	for (int n = 0; n < N / 2; n++) {
		odd.push_back(xin[4 * n + 2]);
		odd.push_back(xin[4 * n + 3]);
	}
	even = fft_radix2(even, N / 2);
	odd = fft_radix2(odd, N / 2);
	for (int k = 0; k < N / 2; k++) {
		double theta = -2.0 * M_PI * k / N;
		auto twr = dag_node(cos(theta));
		auto twi = dag_node(sin(theta));
		auto tr = odd[2 * k] * twr - odd[2 * k + 1] * twi;
		auto ti = odd[2 * k] * twi + odd[2 * k + 1] * twr;
		xout[2 * k] = even[2 * k] + tr;
		xout[2 * (k + N / 2)] = even[2 * k] - tr;
		xout[2 * k + 1] = even[2 * k + 1] + ti;
		xout[2 * (k + N / 2) + 1] = even[2 * k + 1] - ti;
	}
	return xout;
}



std::vector<dag_node> fft_radix4(std::vector<dag_node> xin, int N) {
	if (N == 1) {
		return xin;
	} else if (N == 2) {
		std::vector<dag_node> xout(4);
		auto x0 = xin[0] + xin[2];
		auto x1 = xin[0] - xin[2];
		auto y0 = xin[1] + xin[3];
		auto y1 = xin[1] - xin[3];
		xout[0] = x0;
		xout[1] = y0;
		xout[2] = x1;
		xout[3] = y1;
		return xout;
	}
	std::vector<dag_node> xout(2 * N);
	std::vector<dag_node> even, odd1, odd3;
	for (int n = 0; n < N / 2; n++) {
		even.push_back(xin[4 * n]);
		even.push_back(xin[4 * n + 1]);
	}
	for (int n = 0; n < N / 4; n++) {
		odd1.push_back(xin[8 * n + 2]);
		odd1.push_back(xin[8 * n + 3]);
		odd3.push_back(xin[8 * n + 6]);
		odd3.push_back(xin[8 * n + 7]);
	}
	if (N == 4) {
		even = fft_radix2(even, 2);
		odd1 = fft_radix2(odd1, 1);
		odd3 = fft_radix2(odd3, 1);
	} else if (N == 8) {
		even = fft_radix4(even, 4);
		odd1 = fft_radix2(odd1, 2);
		odd3 = fft_radix2(odd3, 2);
	} else {
		even = fft_radix4(even, N / 2);
		odd1 = fft_radix4(odd1, N / 4);
		odd3 = fft_radix4(odd3, N / 4);
	}
	const auto tw_mult = [N](int k, dag_node r, dag_node i) {
		std::pair<dag_node, dag_node> rc;
		double theta = -2.0 * M_PI * k / N;

		auto twr = dag_node(cos(theta));
		auto twi = dag_node(sin(theta));
		rc.first = r * twr - i * twi;
		rc.second = i * twr + r * twi;
		return rc;
	};
	for (int k = 0; k < N / 4; k++) {
		auto odds1 = tw_mult(k, odd1[2 * k], odd1[2 * k + 1]);
		auto odds3 = tw_mult(3 * k, odd3[2 * k], odd3[2 * k + 1]);
		auto zsr = odds1.first + odds3.first;
		auto zsi = odds1.second + odds3.second;
		auto zdr = odds1.first - odds3.first;
		auto zdi = odds1.second - odds3.second;
		auto ur0 = even[2 * k + 0] + zsr;
		auto ui0 = even[2 * k + 1] + zsi;
		auto ur1 = even[2 * (k + N / 4) + 0] + zdi;
		auto ui1 = even[2 * (k + N / 4) + 1] - zdr;
		auto ur2 = even[2 * k + 0] - zsr;
		auto ui2 = even[2 * k + 1] - zsi;
		auto ur3 = even[2 * (k + N / 4) + 0] - zdi;
		auto ui3 = even[2 * (k + N / 4) + 1] + zdr;
		xout[2 * k] = ur0;
		xout[2 * k + 1] = ui0;
		xout[2 * (k + N / 4)] = ur1;
		xout[2 * (k + N / 4) + 1] = ui1;
		xout[2 * (k + N / 2)] = ur2;
		xout[2 * (k + N / 2) + 1] = ui2;
		xout[2 * (k + 3 * N / 4)] = ur3;
		xout[2 * (k + 3 * N / 4) + 1] = ui3;
	}
	return xout;
}

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
