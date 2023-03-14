#include "fft.hpp"

#include <cmath>

struct cmplx {
	math_vertex x;
	math_vertex y;
	cmplx operator+(cmplx other) {
		x += other.x;
		y += other.y;
		return *this;
	}
	cmplx operator*(cmplx other) {
		auto tmpx = x;
		auto tmpy = y;
		if (true) {
			x = tmpx * other.x - tmpy * other.y;
			y = tmpy * other.x + tmpx * other.y;
		} else {
			auto tmpx = x;
			auto tmpy = y;
			auto a = (tmpx + tmpy);
			auto b = (other.x - other.y);
			auto c = tmpy * other.x;
			auto d = tmpx * other.y;
			x = a * b + d - c;
			y = c + d;
			y = tmpy * other.x + tmpx * other.y;
		}
		return *this;
	}
};

std::vector<math_vertex> fft_singleton(std::vector<math_vertex> xin, int N) {
	const int M = (N - 1) / 2;
	std::vector<math_vertex> xout(2 * N);
	std::vector<math_vertex> tp(M + 1);
	std::vector<math_vertex> tm(M + 1);
	std::vector<math_vertex> up(M + 1);
	std::vector<math_vertex> um(M + 1);
	std::vector<math_vertex> ap(M + 1);
	std::vector<math_vertex> am(M + 1);
	std::vector<math_vertex> bp(M + 1);
	std::vector<math_vertex> bm(M + 1);
	std::vector<std::vector<bool>> pc(M + 1, std::vector<bool>(M + 1, false));
	std::vector<bool> xc(M + 1, false);
	for (int j = 1; j <= M; j++) {
		tp[j] = xin[2 * j] + xin[2 * (N - j)];
		tm[j] = xin[2 * j] - xin[2 * (N - j)];
		up[j] = xin[2 * j + 1] + xin[2 * (N - j) + 1];
		um[j] = xin[2 * j + 1] - xin[2 * (N - j) + 1];
	}
	for (int i = 1; i <= M; i++) {
		am[i] = 0.0;
		bm[i] = 0.0;
		ap[i] = 0.0;
		bp[i] = 0.0;
	}
	xout[0] = math_vertex(0.0);
	xout[1] = math_vertex(0.0);
	std::vector<bool> first(M + 1, true);
	for (int i0 = 1; i0 <= M; i0++) {
		for (int j0 = 1; j0 <= M; j0++) {
			for (int i1 = i0 + 1; i1 <= M; i1++) {
				for (int j1 = j0 + 1; j1 <= M; j1++) {
					{
						double c00 = cos(2.0 * M_PI * i0 * j0 / N);
						double c10 = cos(2.0 * M_PI * i1 * j0 / N);
						double c01 = cos(2.0 * M_PI * i0 * j1 / N);
						double c11 = cos(2.0 * M_PI * i1 * j1 / N);
						if (close2(c00, c11) && close2(c10, c01)) {
							double c0 = c00;
							double c1 = c01;
							math_vertex a = math_vertex(0.5 * (c0 + c1));
							math_vertex b = math_vertex(0.5 * (c0 - c1));
							math_vertex s0 = (tp[j0] + tp[j1]);
							math_vertex s1 = (tp[j0] - tp[j1]);
							math_vertex r0 = (up[j0] + up[j1]);
							math_vertex r1 = (up[j0] - up[j1]);
							bool f = first[i0] && first[i1];
							auto as0 = (f ? xin[0] : math_vertex(0.0)) + a * s0;
							auto ar0 = (f ? xin[1] : math_vertex(0.0)) + a * r0;
							first[i0] = false;
							first[i1] = false;
							ap[i0] += as0 + b * s1;
							ap[i1] += as0 - b * s1;
							bp[i0] += ar0 + b * r1;
							bp[i1] += ar0 - b * r1;
							if (!xc[j0] && !xc[j1]) {
								xout[0] += s0;
								xout[1] += r0;
								xc[j0] = true;
								xc[j1] = true;
							}
							pc[i0][j0] = true;
							pc[i0][j1] = true;
							pc[i1][j0] = true;
							pc[i1][j1] = true;
						}
					}
				}
			}
		}
	}
	for (int i = 1; i <= M; i++) {
		for (int j = 1; j <= M; j++) {
			auto co = cos(2.0 * M_PI * j * i / N);
			if (!pc[i][j]) {
				ap[i] += tp[j] * co;
				bp[i] += up[j] * co;
			}
		}
	}
	for (int i = 1; i <= M; i++) {
		if (!xc[i]) {
			xout[0] += tp[i];
			xout[1] += up[i];
		}
	}
	xout[0] += xin[0];
	xout[1] += xin[1];
	for (int i = 1; i <= M; i++) {
		for (int j = 1; j <= M; j++) {
			auto si = sin(2.0 * M_PI * j * i / N);
			am[i] -= um[j] * si;
			bm[i] -= tm[j] * si;
		}
	}
	for (int i = 1; i <= M; i++) {
		if (first[i]) {
			ap[i] += xin[0];
			bp[i] += xin[1];
		}
	}
	for (int i = 1; i <= M; i++) {
		xout[2 * i] = ap[i] - am[i];
		xout[2 * i + 1] = bp[i] + bm[i];
		xout[2 * (N - i)] = ap[i] + am[i];
		xout[2 * (N - i) + 1] = bp[i] - bm[i];
	}
	return xout;
}

std::vector<math_vertex> fft_radix2(std::vector<math_vertex> xin, int N) {
	if (N == 1) {
		return xin;
	}
	std::vector<math_vertex> xout(2 * N);
	std::vector<math_vertex> even, odd;
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
		auto twr = math_vertex(cos(theta));
		auto twi = math_vertex(sin(theta));
		auto tr = odd[2 * k] * twr - odd[2 * k + 1] * twi;
		auto ti = odd[2 * k] * twi + odd[2 * k + 1] * twr;
		xout[2 * k] = even[2 * k] + tr;
		xout[2 * (k + N / 2)] = even[2 * k] - tr;
		xout[2 * k + 1] = even[2 * k + 1] + ti;
		xout[2 * (k + N / 2) + 1] = even[2 * k + 1] - ti;
	}
	return xout;
}

std::vector<math_vertex> fft_prime_power(int R, std::vector<math_vertex> xin, int N) {
	if (N == 1) {
		return xin;
	}
	const int N1 = R;
	const int N2 = N / R;
	std::vector<math_vertex> xout(2 * N);
	std::vector<std::vector<math_vertex>> sub(N1, std::vector<math_vertex>(2 * N2));
	int begin = -(N1 / 2);
	int end = begin + N1;
	if (false) {
		begin = -(N1 / 2);
		end = begin + N1;
	} else {
		begin = 0;
		end = N1;
	}

	const auto I = [N]( int n ) {
		while( n < 0 ) {
			n += N;
		}
		return n % N;
	};
	const auto I1 = [N1]( int n ) {
		while( n < 0 ) {
			n += N1;
		}
		return n % N1;
	};
	const auto I2 = [N2]( int n ) {
		while( n < 0 ) {
			n += N2;
		}
		return n % N2;
	};

	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = begin; n1 < end; n1++) {
			sub[I1(n1)][2 * n2] = xin[2 * I(n2 * N1 + n1)];
			sub[I1(n1)][2 * n2 + 1] = xin[2 * I(n2 * N1 + n1) + 1];
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
			x.x = 0.0;
			x.y = 0.0;
			for (int n1 = begin; n1 < end; n1++) {
				const double theta = -2.0 * M_PI * (n1 * (N2 * k1 + k2)) / N;
				t.x = (sub[I1(n1)][2 * k2]);
				t.y = (sub[I1(n1)][2 * k2 + 1]);
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

std::vector<math_vertex> fft_radix4(std::vector<math_vertex> xin, int N) {
	if (N == 1) {
		return xin;
	} else if (N == 2) {
		std::vector<math_vertex> xout(4);
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
	const auto I = [N](int n) {
		while( n < N) {
			n+= N;
		}
		return n % N;
	};
	std::vector<math_vertex> xout(2 * N);
	std::vector<math_vertex> even, odd1, odd3;
	for (int n = 0; n < N / 2; n++) {
		even.push_back(xin[4 * n]);
		even.push_back(xin[4 * n + 1]);
	}
	for (int n = 0; n < N / 4; n++) {
		odd1.push_back(xin[2 * (4 * n + 1)]);
		odd1.push_back(xin[2 * (4 * n + 1) + 1]);
		odd3.push_back(xin[2 * (4 * n + 3)]);
		odd3.push_back(xin[2 * (4 * n + 3) + 1]);
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
	const auto tw_mult = [N](int k, math_vertex r, math_vertex i) {
		std::pair<math_vertex, math_vertex> rc;
		double theta = -2.0 * M_PI * k / N;

		auto twr = math_vertex(cos(theta));
		auto twi = math_vertex(sin(theta));
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
