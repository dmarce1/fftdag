#include "fft.hpp"

#include <cmath>
#include <complex>
#include <fftw3.h>

static int mod_pow(int a, int b, int m) {
	int rc = 1;
	int apow = a;
	while (b) {
		if (b & 1) {
			rc = ((rc % m) * (apow % m)) % m;
		}
		b >>= 1;
		apow = ((apow % m) * (apow % m)) % m;
	}
	return rc;
}

static int mod_inv(int a, int m) {
	return mod_pow(a, m - 2, m);
}

static int generator(int N) {
	static thread_local std::unordered_map<int, int> values;
	auto i = values.find(N);
	if (i == values.end()) {
		for (int g = 2;; g++) {
			std::set<int> I;
			bool fail = false;
			for (int m = 0; m < N - 2; m++) {
				int n = mod_pow(g, m, N);
				if (I.find(n) == I.end()) {
					I.insert(n);
				} else {
					fail = true;
					break;
				}
			}
			if (!fail) {
				values[N] = g;
				i = values.find(N);
				break;
			}
		}
	}
	return i->second;
}

std::vector<int> raders_ginvq(int N) {
	const int g = generator(N);
	std::vector<int> ginvq;
	for (int q = 0; q < N - 1; q++) {
		ginvq.push_back(mod_inv(mod_pow(g, q, N), N));
	}
	return ginvq;
}

const std::vector<int> raders_gq(int N) {
	const int g = generator(N);
	std::vector<int> gq;
	for (int q = 0; q < N - 1; q++) {
		gq.push_back(mod_pow(g, q, N));
	}
	return gq;
}

void fftw(std::vector<std::complex<double>>& x) {
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, fftw_complex*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);
		plans[N] = fftw_plan_dft_1d(N, in[N], out[N], FFTW_FORWARD, FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n][0] = x[n].real();
		i[n][1] = x[n].imag();
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		x[n].real(o[n][0]);
		x[n].imag(o[n][1]);
	}

}

const std::vector<std::complex<double>> twiddles(int N) {
	std::vector<std::complex<double>> tw(N);
	for (int k = 0; k < N; k++) {
		tw[k] = std::polar(1.0, -2 * M_PI * k / N);
	}
	return tw;
}

const std::vector<std::complex<double>> raders_four_twiddle(int N) {
	std::vector<std::complex<double>> b(N - 1);
	const auto tws = twiddles(N);
	const auto ginvq = raders_ginvq(N);
	for (int q = 0; q < N - 1; q++) {
		b[q] = tws[ginvq[q]];
	}
	fftw(b);
	return b;
}

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

std::vector<math_vertex> fft(std::vector<math_vertex> xin, int N, int opts) {
	if (opts & FFT_REAL) {
		auto tmp = std::move(xin);
		xin.resize(2 * N);
		for (int n = 0; n < N; n++) {
			xin[2 * n] = tmp[n];
			xin[2 * n + 1] = math_vertex(0.0);
		}
	}
	if (opts & FFT_SYM) {
		for (int n = 1; n < N - n; n++) {
			xin[2 * (N - n)] = xin[2 * n];
			xin[2 * (N - n) + 1] = xin[2 * n + 1];
		}
	}
	std::vector<math_vertex> xout;
	auto pfac = prime_factorization(N);
	if (pfac.size() == 1) {
		int R = pfac.begin()->first;
		if (R == 2) {
			xout = fft_radix4(xin, N);
		} else if (R < 29) {
			xout = fft_prime_power(pfac.begin()->first, xin, N);
		} else {
			xout = fft_raders(xin, N);
		}
	} else {
		int N1 = 1;
		int N2 = 1;
		bool flag = false;
		for (auto i = pfac.begin(); i != pfac.end(); i++) {
			if (flag) {
				N2 *= std::pow(i->first, i->second);
			} else {
				N1 *= std::pow(i->first, i->second);
			}
			flag = !flag;
		}
		xout = fft_prime_factor(N1, N2, xin);
	}
	if (opts & (FFT_REAL | FFT_SYM)) {
		auto tmp = std::move(xout);
		xout.resize(N);
		xout[0] = tmp[0];
		if (N % 2 == 0) {
			xout[N / 2] = tmp[2 * (N / 2)];
		}
		for (int n = 1; n < N - n; n++) {
			xout[n] = tmp[2 * n];
			xout[N - n] = tmp[2 * n + 1];
		}
		if (opts & FFT_SYM) {
			xout.resize(N / 2 + 1);
		}
	}
	return std::move(xout);
}

std::vector<math_vertex> fft_raders(std::vector<math_vertex> xin, int N) {
	std::vector<math_vertex> xout(2 * N);
	const auto& b = raders_four_twiddle(N);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	std::vector<math_vertex> a(2 * (N - 1));
	for (int q = 0; q < N - 1; q++) {
		a[2 * q] = xin[2 * gq[q]];
		a[2 * q + 1] = xin[2 * gq[q] + 1];
	}
	a = fft(a, N - 1, false);
	for (int q = 0; q < N - 1; q++) {
		cmplx A, B, C;
		A.x = a[2 * q];
		A.y = a[2 * q + 1];
		B.x = b[q].real();
		B.y = b[q].imag();
		C = A * B;
		C.y = -C.y;
		a[2 * q] = C.x;
		a[2 * q + 1] = C.y;
	}
	a = fft(a, N - 1, false);
	const auto Nm1inv = 1.0 / (N - 1.0);
	for (int q = 0; q < N - 1; q++) {
		a[2 * q + 1] = a[2 * q + 1] * Nm1inv;
		a[2 * q] = a[2 * q] * Nm1inv;
		a[2 * q + 1] = -a[2 * q + 1];
	}
	xout[0] = 0.0;
	xout[1] = 0.0;
	for (int n = 0; n < N; n++) {
		xout[0] += xin[2 * n];
		xout[1] += xin[2 * n + 1];
	}
	for (int p = 0; p < N - 1; p++) {
		xout[2 * ginvq[p]] = xin[0] + a[2 * p];
		xout[2 * ginvq[p] + 1] = xin[1] + a[2 * p + 1];
	}
	return std::move(xout);
}

std::vector<math_vertex> fft_prime_factor(int N1, int N2, std::vector<math_vertex> xin) {
	int N = N1 * N2;
	std::vector<math_vertex> xout(2 * N);
	std::vector<std::vector<math_vertex>> sub1(N1, std::vector<math_vertex>(2 * N2));
	std::vector<std::vector<math_vertex>> sub2(N2, std::vector<math_vertex>(2 * N1));
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			sub1[n1][2 * n2] = xin[2 * ((n1 * N2 + n2 * N1) % N)];
			sub1[n1][2 * n2 + 1] = xin[2 * ((n1 * N2 + n2 * N1) % N) + 1];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		sub1[n1] = fft(sub1[n1], N2);
	}
	for (int n1 = 0; n1 < N1; n1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			sub2[k2][2 * n1] = sub1[n1][2 * k2];
			sub2[k2][2 * n1 + 1] = sub1[n1][2 * k2 + 1];
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		sub2[k2] = fft(sub2[k2], N1, false);
	}
	for (int k = 0; k < N; k++) {
		const int k1 = k % N1;
		const int k2 = k % N2;
		xout[2 * k] = sub2[k2][2 * k1];
		xout[2 * k + 1] = sub2[k2][2 * k1 + 1];
	}
	return xout;
}

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

std::vector<math_vertex> fft_prime_power(int R, std::vector<math_vertex> xin, int N) {
	if (N == R) {
		return fft_singleton(xin, N);
	}
	const int N1 = R;
	const int N2 = N / R;
	std::vector<math_vertex> xout(2 * N);
	std::vector<std::vector<math_vertex>> sub(N1, std::vector<math_vertex>(2 * N2));
	int begin = -(N1 / 2);
	int end = begin + N1;
	if (true) {
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
	{
		for (int k1 = begin; k1 < end; k1++) {
			xout[2 * I(k1 * N2)] = sub[I1(0)][0];
			xout[2 * I(k1 * N2) + 1] = sub[I1(0)][1];
		}
		for (int n1 = 1; n1 < end; n1++) {
			cmplx x;
			cmplx t_0, t_1;
			cmplx w;
			for (int k1 = begin; k1 < end; k1++) {
				const double theta_0 = -2.0 * M_PI * n1 * k1 / N1;
				t_0.x = sub[I1(n1)][0];
				t_0.y = sub[I1(n1)][1];
				t_1.x = sub[I1(-n1)][0];
				t_1.y = sub[I1(-n1)][1];
				w.x = cos(theta_0);
				w.y = sin(theta_0);
				xout[2 * I(k1 * N2)] += w.x * (t_0.x + t_1.x);
				xout[2 * I(k1 * N2) + 1] += w.x * (t_0.y + t_1.y);
				xout[2 * I(k1 * N2)] += w.y * (t_1.y - t_0.y);
				xout[2 * I(k1 * N2) + 1] += w.y * (t_0.x - t_1.x);
			}
		}
	}
	for (int k2 = 1; k2 < N2; k2++) {
		for (int k1 = begin; k1 < end; k1++) {
			xout[2 * I(k2 + k1 * N2)] = sub[I1(0)][2 * k2];
			xout[2 * I(k2 + k1 * N2) + 1] = sub[I1(0)][2 * k2 + 1];
		}
		for (int k1 = begin + 1; k1 < end; k1++) {
			cmplx x;
			cmplx t_0, t_1;
			cmplx w;
			x.x = x.y = 0.0;
			for (int n1 = 1; n1 < end; n1++) {
				const double theta_0 = -2.0 * M_PI * (n1 * (N2 * k1 + k2)) / N;
				t_0.x = sub[I1(n1)][2 * k2];
				t_0.y = sub[I1(n1)][2 * k2 + 1];
				t_1.x = sub[I1(-n1)][2 * k2];
				t_1.y = sub[I1(-n1)][2 * k2 + 1];
				w.x = cos(theta_0);
				w.y = sin(theta_0);
				x.x = x.x + w.x * (t_0.x + t_1.x) + w.y * (t_1.y - t_0.y);
				x.y = x.y + w.x * (t_0.y + t_1.y) + w.y * (t_0.x - t_1.x);
			}
			xout[2 * I(k2 + k1 * N2)] += x.x;
			xout[2 * I(k2 + k1 * N2) + 1] += x.y;
			xout[2 * I(k2 + begin * N2)] -= x.x;
			xout[2 * I(k2 + begin * N2) + 1] -= x.y;
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
	even = fft_radix4(even, N / 2);
	odd1 = fft_radix4(odd1, N / 4);
	odd3 = fft_radix4(odd3, N / 4);
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
