#include "fft.hpp"

#include <cmath>
#include <complex>
#include <fftw3.h>

std::vector<cmplx> fft_prime_power(int R, std::vector<cmplx> xin, int N, int opts = 0);
std::vector<cmplx> fft_radix4(std::vector<cmplx> xin, int N, int opts = 0);
std::vector<cmplx> fft_singleton(std::vector<cmplx> xin, int N, int opts = 0);
std::vector<cmplx> fft_prime_factor(int N1, int N2, std::vector<cmplx> xin, int opts = 0);
std::vector<cmplx> fft_raders(std::vector<cmplx> xin, int N, bool padded, int opts = 0);
std::vector<cmplx> fft(std::vector<cmplx>& xin, int N, int opts = 0);

std::vector<math_vertex> fft(std::vector<math_vertex>& xin, int N, int opts) {
	std::vector<cmplx> x(N, cmplx( { 0.0, 0.0 }));
	std::vector<math_vertex> xout;
	if (opts & FFT_INV) {
		if (opts & FFT_REAL) {
			x[0].x = xin[0];
			if (N % 2 == 0) {
				x[N / 2].x = xin[N / 2];
			}
			for (int n = 1; n < N - n; n++) {
				x[N - n].x = xin[n];
				x[N - n].y = xin[N - n];
				x[n] = x[N - n].conj();
			}
		} else if (opts & FFT_DCT2) {
			x[0].x = xin[0];
			for (int n = 1; n < N / 4; n++) {
				x[N - n].x = xin[n];
				x[N - (n + N / 2)].x = -xin[n];
			}
			for (int n = 1; n < N - n; n++) {
				x[n].x = x[N - n].x;
			}
		} else {
			x[0].x = xin[0];
			x[0].y = xin[1];
			for (int n = 1; n < N; n++) {
				x[n].x = xin[2 * (N - n)];
				x[n].y = xin[2 * (N - n) + 1];
			}
		}
	} else {
		if (opts & FFT_REAL) {
			for (int n = 0; n < N; n++) {
				x[n].x = xin[n];
			}
		} else if (opts & FFT_DCT2) {
			for (int n = 1; n < N - n; n += 2) {
				x[n].x = xin[(n - 1) / 2];
				x[N - n].x = xin[(n - 1) / 2];
			}
		} else {
			for (int n = 0; n < N; n++) {
				x[n].x = xin[2 * n];
				x[n].y = xin[2 * n + 1];
			}
		}
	}

	x = fft(x, N);

	if (opts & FFT_INV) {
		if (opts & FFT_REAL) {
			xout.resize(N);
			for (int n = 0; n < N; n++) {
				xout[n] = x[n].x;
			}
		} else if (opts & FFT_DCT2) {
			xout.resize(N / 4);
			for (int n = 0; n < xout.size(); n++) {
				xout[n] = x[2 * n + 1].x;
			}
		} else {
			xout.resize(2 * N);
			for (int n = 0; n < N; n++) {
				xout[2 * n] = x[n].x;
				xout[2 * n + 1] = x[n].y;
			}
		}
	} else {
		if (opts & FFT_REAL) {
			xout.resize(N);
			xout[0] = x[0].x;
			if (N % 2 == 0) {
				xout[N / 2] = x[N / 2].x;
			}
			for (int n = 1; n < N - n; n++) {
				xout[n] = x[n].x;
				xout[N - n] = x[n].y;
			}
		} else if (opts & FFT_DCT2) {
			xout.resize(N / 4);
			for (int n = 0; n < xout.size(); n++) {
				xout[n] = x[n].x;
			}
		} else {
			xout.resize(2 * N);
			for (int n = 0; n < N; n++) {
				xout[2 * n] = x[n].x;
				xout[2 * n + 1] = x[n].y;
			}
		}
	}
	return std::move(xout);
}

std::vector<cmplx> fft(std::vector<cmplx>& xin, int N, int opts) {
	if (N == 1) {
		return xin;
	}
	std::vector<cmplx> xout;
	auto pfac = prime_factorization(N);
	if (pfac.size() == 1) {
		int R = pfac.begin()->first;
		int a = pfac.begin()->second;
		if (a == 1) {
			if (R == 2) {
				xout = fft_radix4(xin, N);
			} else if (R < 29) {
				xout = fft_singleton(xin, N);
			} else {
				xout = fft_raders(xin, N, false);
			}
		} else {
			if (R == 2) {
				xout = fft_modsplit(xin, N);
			} else if (R < 29) {
				xout = fft_prime_power(pfac.begin()->first, xin, N);
			}
		}
	} else {
		{
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
	}
	return std::move(xout);
}

std::vector<cmplx> fft_radix4(std::vector<cmplx> xin, int N, int opts) {
	static cmplx I( { 0.0, 1.0 });
	if (N == 1) {
		return xin;
	} else if (N == 2) {
		std::vector<cmplx> xout(2);
		auto x0 = xin[0] + xin[1];
		auto x1 = xin[0] - xin[1];
		xout[0] = x0;
		xout[1] = x1;
		return xout;
	}
	std::vector<cmplx> xout(N);
	std::vector<cmplx> even, odd1, odd3;
	for (int n = 0; n < N / 2; n++) {
		even.push_back(xin[2 * n]);
	}
	for (int n = 0; n < N / 4; n++) {
		odd1.push_back(xin[4 * n + 1]);
		odd3.push_back(xin[(4 * n - 1 + N) % N]);
	}
	even = fft(even, N / 2, opts);
	odd1 = fft(odd1, N / 4, opts);
	odd3 = fft(odd3, N / 4, opts);
	for (int k = 0; k < N / 4; k++) {
		auto t1 = twiddle(+k, N) * odd1[k];
		auto t3 = twiddle(-k, N) * odd3[k];
		auto zs = t1 + t3;
		auto zd = t1 - t3;
		auto ur0 = even[k + 0 / 4] + 1 * zs;
		auto ur1 = even[k + N / 4] - I * zd;
		auto ur2 = even[k + 0 / 4] - 1 * zs;
		auto ur3 = even[k + N / 4] + I * zd;
		xout[k + 0 * N / 4] = ur0;
		xout[k + 1 * N / 4] = ur1;
		xout[k + 2 * N / 4] = ur2;
		xout[k + 3 * N / 4] = ur3;
	}
	return xout;
}

std::vector<cmplx> fft_prime_power(int R, std::vector<cmplx> xin, int N, int opts) {
	static cmplx I( { 0.0, 1.0 });
	if (N == R) {
		return fft_singleton(xin, N, opts);
	}
	const int N1 = R;
	const int N2 = N / R;
	std::vector<cmplx> xout(N);
	std::vector<std::vector<cmplx>> sub(N1, std::vector<cmplx>(N2));
	int begin = -(N1 / 2);
	int end = begin + N1;
	const auto I0 = [N]( int n ) {
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
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = begin; n1 < end; n1++) {
			sub[I1(n1)][n2] = xin[I0(n2 * N1 + n1)];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		sub[n1] = fft(sub[n1], N2, opts);
	}
	{
		for (int k1 = begin; k1 < end; k1++) {
			xout[I0(k1 * N2)] = sub[I1(0)][0];
		}
		for (int n1 = 1; n1 < end; n1++) {
			for (int k1 = begin; k1 < end; k1++) {
				const auto w = twiddle(n1 * k1, N1);
				const auto t_0 = sub[I1(n1)][0];
				const auto t_1 = sub[I1(-n1)][0];
				xout[I0(k1 * N2)].x += w.x * (t_0.x + t_1.x);
				xout[I0(k1 * N2)].y += w.x * (t_0.y + t_1.y);
				xout[I0(k1 * N2)].x += w.y * (t_1.y - t_0.y);
				xout[I0(k1 * N2)].y += w.y * (t_0.x - t_1.x);
			}
		}
	}
	for (int k2 = 1; k2 < N2; k2++) {
		for (int k1 = begin; k1 < end; k1++) {
			xout[I0(k2 + k1 * N2)] = sub[I1(0)][k2];
		}
		for (int k1 = begin + 1; k1 < end; k1++) {
			cmplx x = cmplx( { 0.0, 0.0 });
			for (int n1 = 1; n1 < end; n1++) {
				const auto w = twiddle(n1 * (N2 * k1 + k2), N);
				const auto t_0 = sub[I1(n1)][k2];
				const auto t_1 = sub[I1(-n1)][k2];
				x.x = x.x + w.x * (t_0.x + t_1.x) + w.y * (t_1.y - t_0.y);
				;
				x.y = x.y + w.x * (t_0.y + t_1.y) + w.y * (t_0.x - t_1.x);
			}
			xout[I0(k2 + k1 * N2)].x += x.x;
			xout[I0(k2 + k1 * N2)].y += x.y;
			xout[I0(k2 + begin * N2)].x -= x.x;
			xout[I0(k2 + begin * N2)].y -= x.y;
		}
	}
	return xout;
}

std::vector<cmplx> fft_singleton(std::vector<cmplx> xin, int N, int opts) {
	const int M = (N - 1) / 2;
	std::vector<cmplx> xout(N);
	std::vector<cmplx> tp(M + 1);
	std::vector<cmplx> tm(M + 1);
	std::vector<math_vertex> up(M + 1);
	std::vector<math_vertex> um(M + 1);
	std::vector<math_vertex> ap(M + 1);
	std::vector<math_vertex> am(M + 1);
	std::vector<math_vertex> bp(M + 1);
	std::vector<math_vertex> bm(M + 1);
	std::vector<std::vector<bool>> pc(M + 1, std::vector<bool>(M + 1, false));
	std::vector<bool> xc(M + 1, false);
	for (int j = 1; j <= M; j++) {
		tp[j] = xin[j] + xin[N - j];
		tm[j] = xin[j] - xin[N - j];
	}
	for (int i = 1; i <= M; i++) {
		am[i] = 0.0;
		bm[i] = 0.0;
		ap[i] = 0.0;
		bp[i] = 0.0;
	}
	xout[0] = cmplx( { 0.0, 0.0 });
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
							cmplx s0 = (tp[j0] + tp[j1]);
							cmplx s1 = (tp[j0] - tp[j1]);
							bool f = first[i0] && first[i1];
							auto as0 = (f ? xin[0] : cmplx( { 0.0, 0.0 })) + a * s0;
							first[i0] = false;
							first[i1] = false;
							ap[i0] += as0.x + b * s1.x;
							ap[i1] += as0.x - b * s1.x;
							bp[i0] += as0.y + b * s1.y;
							bp[i1] += as0.y - b * s1.y;
							if (!xc[j0] && !xc[j1]) {
								xout[0] += s0;
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
				ap[i] += tp[j].x * co;
				bp[i] += tp[j].y * co;
			}
		}
	}
	for (int i = 1; i <= M; i++) {
		if (!xc[i]) {
			xout[0] += tp[i];
		}
	}
	xout[0] += xin[0];
	for (int i = 1; i <= M; i++) {
		for (int j = 1; j <= M; j++) {
			auto si = sin(2.0 * M_PI * j * i / N);
			am[i] -= tm[j].y * si;
			bm[i] -= tm[j].x * si;
		}
	}
	for (int i = 1; i <= M; i++) {
		if (first[i]) {
			ap[i] += xin[0].x;
			bp[i] += xin[0].y;
		}
	}
	for (int i = 1; i <= M; i++) {
		xout[i].x = ap[i] - am[i];
		xout[i].y = bp[i] + bm[i];
		xout[N - i].x = ap[i] + am[i];
		xout[N - i].y = bp[i] - bm[i];
	}
	return std::move(xout);
}

std::vector<cmplx> fft_prime_factor(int N1, int N2, std::vector<cmplx> xin, int opts) {
	int N = N1 * N2;
	std::vector<cmplx> xout(N);
	std::vector<std::vector<cmplx>> sub1(N1, std::vector<cmplx>(N2));
	std::vector<std::vector<cmplx>> sub2(N2, std::vector<cmplx>(N1));
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			sub1[n1][n2] = xin[(n1 * N2 + n2 * N1) % N];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		sub1[n1] = fft(sub1[n1], N2, opts);
	}
	for (int n1 = 0; n1 < N1; n1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			sub2[k2][n1] = sub1[n1][k2];
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		sub2[k2] = fft(sub2[k2], N1, opts);
	}
	for (int k = 0; k < N; k++) {
		const int k1 = k % N1;
		const int k2 = k % N2;
		xout[k] = sub2[k2][k1];
	}
	return xout;
}

std::vector<cmplx> fft_raders(std::vector<cmplx> xin, int N, bool padded, int opts) {
	std::vector<cmplx> xout(N);
	int M;
	if (padded) {
		M = 1;
		while (M < 2 * (N - 1)) {
			M *= 2;
		}
	} else {
		M = N - 1;
	}
	const auto& b = padded ? raders_four_twiddle(N, M) : raders_four_twiddle(N);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	std::vector<cmplx> a(M, cmplx( { 0.0, 0.0 }));
	for (int q = 0; q < N - 1; q++) {
		a[q] = xin[gq[q]];
	}
	a = fft(a, M, 0);
	for (int q = 0; q < M; q++) {
		a[q] = (a[q] * cmplx(b[q])).conj();
	}
	a = fft(a, M, 0);
	for (int q = 0; q < M; q++) {
		a[q] = a[q].conj();
	}
	xout[0] = cmplx( { 0.0, 0.0 });
	for (int n = 0; n < N; n++) {
		xout[0] += xin[n];
	}
	for (int p = 0; p < N - 1; p++) {
		xout[ginvq[p]] = xin[0] + a[p];
	}
	return std::move(xout);
}
