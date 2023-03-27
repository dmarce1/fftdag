#include "fft.hpp"
#include "util.hpp"

std::vector<cmplx> convolve_karatsuba(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> convolve_fast(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> convolve_fft(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<std::vector<cmplx>> convolve_aperiodic(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h);
std::vector<cmplx> operator*(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> operator*(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> convolve_toomcook(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> operator+(std::vector<cmplx> x1, std::vector<cmplx> x2) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] = x1[n] + x2[n];
	}
	return x1;
}
std::vector<cmplx> operator-(std::vector<cmplx> x1, std::vector<cmplx> x2) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] = x1[n] - x2[n];
	}
	return x1;
}
std::vector<cmplx> operator-(std::vector<cmplx> x1) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] = -x1[n];
	}
	return x1;
}
std::vector<std::complex<double>> operator+(std::vector<std::complex<double>> x1, std::vector<std::complex<double>> x2) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] = x1[n] + x2[n];
	}
	return x1;
}
std::vector<std::complex<double>> operator-(std::vector<std::complex<double>> x1, std::vector<std::complex<double>> x2) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] = x1[n] - x2[n];
	}
	return x1;
}
std::vector<std::complex<double>> operator*(double a, std::vector<std::complex<double>> x1) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] *= a;
	}
	return x1;
}
std::vector<cmplx> operator*(double a, std::vector<cmplx> x1) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] *= a;
	}
	return x1;
}
std::vector<cmplx> operator*(std::vector<cmplx> x1, double a) {
	return a * x1;
}

std::vector<std::vector<cmplx>> convolve_radix2(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h) {
	int N = x.size();
	if( N == 1 ) {
		return std::vector<std::vector<cmplx>>(1, convolve_fast(x[0], h[0]));
	}
	std::vector<std::vector<cmplx>> y(2 * N - 1);
	std::vector<std::vector<cmplx>> xlo(N / 2), xhi(N / 2), xmd(N / 2);
	std::vector<std::vector<std::complex<double>>> hlo(N / 2), hhi(N / 2), hmd(N / 2);
	for (int n = 0; n < N / 2; n++) {
		xlo[n] = x[n];
		xhi[n] = x[n + N / 2];
		xmd[n] = xlo[n] + xhi[n];
		hlo[n] = h[n];
		hhi[n] = h[n + N / 2];
		hmd[n] = hlo[n] + hhi[n];
	}
	auto ylo = convolve_aperiodic(xlo, hlo);
	auto ymd = convolve_aperiodic(xmd, hmd);
	auto yhi = convolve_aperiodic(xhi, hhi);
	for (int n = 0; n < N - 1; n++) {
		ymd[n] = ymd[n] - ylo[n] - yhi[n];
	}
	for (int n = 0; n < 2 * N - 1; n++) {
		y[n] = std::vector<cmplx>(x[0].size(), cmplx( {0.0,0.0}));
	}
	for (int n = 0; n < N - 1; n++) {
		y[n] = y[n] + ylo[n];
		y[n + N / 2] = y[n + N / 2] + ymd[n];
		y[n + N] = y[n + N] + yhi[n];
	}
	for (int n = 0; n < 2*N-1; n++) {
		for( auto& yy : y[n]) {
			yy.set_goal();
		}
	}
	return y;
}

std::vector<std::vector<cmplx>> convolve_radix3(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h) {
	int N = x.size();
	if( N == 1 ) {
		return std::vector<std::vector<cmplx>>(1, convolve_fast(x[0], h[0]));
	}
	std::vector<std::vector<cmplx>> y(2 * N - 1);
	std::vector<std::vector<cmplx>> a2(N / 3), a3(N / 3), a4(N / 3), a5(N / 3), a6(N / 3);
	std::vector<std::vector<std::complex<double>>> b0(N / 3), b1(N / 3), b2(N / 3), b3(N / 3), b4(N / 3);
	std::vector<std::vector<cmplx>> u1(2 * N / 3 - 1), u2(2 * N / 3 - 1), u3(2 * N / 3 - 1), u4(2 * N / 3 - 1), u5(2 * N / 3 - 1);
	for (int n = 0; n < N / 3; n++) {
		auto x0 = x[n];
		auto x1 = x[n + N / 3];
		auto x2 = x[n + 2 * N / 3];
		auto h0 = h[n];
		auto h1 = h[n + N / 3];
		auto h2 = h[n + 2 * N / 3];
		auto a0 = x1 + x2;
		auto a1 = x2 - x1;
		a2[n] = x0;
		a3[n] = x0 + a0;
		a4[n] = x0 + a1;
		a5[n] = 2.0 * a0 + a1 + a3[n];
		a6[n] = x2;
		b0[n] = 0.5 * h0;
		b1[n] = 0.5 * (h0 + h1 + h2);
		b2[n] = (1.0 / 6.0) * (h0 - h1 + h2);
		b3[n] = (1.0 / 6.0) * (h0 + 2.0 * h1 + 4.0 * h2);
		b4[n] = h2;
	}
	auto m0 = convolve_aperiodic(a2, b0);
	auto m1 = convolve_aperiodic(a3, b1);
	auto m2 = convolve_aperiodic(a4, b2);
	auto m3 = convolve_aperiodic(a5, b3);
	auto m4 = convolve_aperiodic(a6, b4);
	for (int n = 0; n < 2 * N / 3 - 1; n++) {
		auto u0 = 2.0 * m4[n];
		u1[n] = 2.0 * m1[n];
		u2[n] = 2.0 * m0[n];
		u3[n] = 2.0 * m2[n];
		u4[n] = u0 - m0[n] - m3[n];
		u5[n] = m1[n] + m2[n];
	}

	for (int n = 0; n < 2 * N - 1; n++) {
		y[n] = std::vector<cmplx>(x[0].size(), cmplx( {0.0,0.0}));
	}
	for (int n = 0; n < 2 * N / 3 - 1; n++) {
		y[n + 0 * N / 3] = y[n + 0 * N / 3] + u2[n];
		y[n + 1 * N / 3] = y[n + 1 * N / 3] + u1[n] - u3[n] + u4[n];
		y[n + 2 * N / 3] = y[n + 2 * N / 3] + -u2[n] + u3[n] + u5[n] - m4[n];
		y[n + 3 * N / 3] = y[n + 3 * N / 3] + -u4[n] - u5[n];
		y[n + 4 * N / 3] = y[n + 4 * N / 3] + m4[n];

	}
	for (int n = 0; n < 2*N-1; n++) {
		for( auto& yy : y[n]) {
			yy.set_goal();
		}
	}
	return y;
}

std::vector<std::vector<cmplx>> convolve_fast_2(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h) {
	std::vector<std::vector<cmplx>> y(2, std::vector<cmplx>(x[0].size()));
	std::vector<cmplx> a0 = x[0] + x[1];
	std::vector<cmplx> a1 = x[0] - x[1];
	std::vector<std::complex<double>> b0 = 0.5 * (h[0] + h[1]);
	std::vector<std::complex<double>> b1 = 0.5 * (h[0] - h[1]);
	std::vector<cmplx> m0 = a0 * b0;
	std::vector<cmplx> m1 = a1 * b1;
	y[0] = m0 + m1;
	y[1] = m0 - m1;
	return std::move(y);
}

std::vector<std::vector<cmplx>> convolve_fast_3(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h) {
	std::vector<std::vector<cmplx>> y(3);
	auto a0 = x[0] + x[1] + x[2];
	auto a1 = x[0] - x[2];
	auto a2 = x[1] - x[2];
	auto a3 = a1 + a2;
	auto b0 = (1.0 / 3.0) * (h[0] + h[1] + h[2]);
	auto b1 = h[0] - h[2];
	auto b2 = h[1] - h[2];
	auto b3 = (1.0 / 3.0) * (b1 + b2);
	auto m0 = a0 * b0;
	auto m1 = a1 * b1;
	auto m2 = a2 * b2;
	auto m3 = a3 * b3;
	auto u0 = m1 - m3;
	auto u1 = m2 - m3;
	y[0] = m0 + u0;
	y[1] = m0 - u0 - u1;
	y[2] = m0 + u1;
	return std::move(y);
}

std::vector<std::vector<cmplx>> convolve_fast_4(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h) {
	std::vector<std::vector<cmplx>> y(4);
	auto a0 = x[0] + x[2];
	auto a1 = x[1] + x[3];
	auto a2 = a0 + a1;
	auto a3 = a0 - a1;
	auto a4 = x[0] - x[2];
	auto a5 = x[1] - x[3];
	auto a6 = a4 + a5;
	auto b0 = h[0] + h[2];
	auto b1 = h[1] + h[3];
	auto b2 = 0.25 * (b0 + b1);
	auto b3 = 0.25 * (b0 - b1);
	auto b4 = 0.5 * ((h[0] - h[2]) - (h[1] - h[3]));
	auto b5 = 0.5 * ((h[0] - h[2]) + (h[1] - h[3]));
	auto b6 = 0.5 * (h[0] - h[2]);
	auto m0 = a2 * b2;
	auto m1 = a3 * b3;
	auto m2 = a4 * b4;
	auto m3 = a5 * b5;
	auto m4 = a6 * b6;
	auto u0 = m0 + m1;
	auto u1 = m0 - m1;
	auto u2 = m4 - m3;
	auto u3 = m4 - m2;
	y[0] = u0 + u2;
	y[1] = u1 + u3;
	y[2] = u0 - u2;
	y[3] = u1 - u3;
	return std::move(y);
}

std::vector<std::vector<cmplx>> convolve_fast_5(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h) {
	std::vector<std::vector<cmplx>> y(5);
	auto a0 = x[0] - x[4];
	auto a1 = x[1] - x[4];
	auto a2 = a0 + a1;
	auto a3 = x[2] - x[4];
	auto a4 = x[3] - x[4];
	auto a5 = a3 + a4;
	auto a6 = a0 - a3;
	auto a7 = a1 - a4;
	auto a8 = a2 - a5;
	auto a9 = x[0] + x[1] + x[2] + x[3] + x[4];
	auto b0 = h[0] - h[2] + h[3] - h[4];
	auto b1 = h[1] - h[2] + h[3] - h[4];
	auto b2 = 0.2 * (-2.0*h[0] - 2.0*h[1] + 3.0*h[2] - 2.0*h[3] + 3.0*h[4]);
	auto b3 = h[1] - h[0] - h[2] + h[3];
	auto b4 = h[1] - h[0] - h[2] + h[4];
	auto b5 = 0.2 * (3.0*h[0] - 2.0*h[1] + 3.0*h[2] - 2.0*h[3] - 2.0*h[4]);
	auto b6 = h[3] - h[2];
	auto b7 = h[1] - h[2];
	auto b8 = 0.2*(4.0 *h[2] - h[3] - h[4] - h[0] - h[1]);
	auto b9 = 0.2*(h[0] + h[1] + h[2] + h[3] + h[4]);
	auto m0 = a0 * b0;
	auto m1 = a1 * b1;
	auto m2 = a2 * b2;
	auto m3 = a3 * b3;
	auto m4 = a4 * b4;
	auto m5 = a5 * b5;
	auto m6 = a6 * b6;
	auto m7 = a7 * b7;
	auto m8 = a8 * b8;
	auto m9 = a9 * b9;
	auto u0 = m0 + m2;
	auto u1 = m1 + m2;
	auto u2 = m3 + m5;
	auto u3 = m4 + m5;
	auto u4 = m6 + m8;
	auto u5 = m7 + m8;
	y[0] = u0 - u4 + m9;
	y[1] = m9 - u0 - u1 - u2 - u3;
	y[2] = u3 + u5 + m9;
	y[3] = u2 + u4 + m9;
	y[4] = u1 - u5 + m9;
	return y;
}

std::vector<std::vector<cmplx>> convolve_fast_7(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h) {
	std::vector<std::vector<cmplx>> y(7);
	auto a0 = x[1] - x[6];
	auto a1 = x[2] - x[6];
	auto a2 = x[4] - x[6];
	auto a3 = x[5] - x[6];
	auto a4 = a0 + a1;
	auto a5 = a1 - a0;
	auto a6 = a2 + a3;
	auto a7 = a3 - a2;
	auto a8 = x[0] - x[6];
	auto a9 = a8 + a4;
	auto a10 = a8 + a5;
	auto a11 = a9 + 2.0 * a4 + a5;
	auto a12 = a1;
	auto a13 = x[3] - x[6];
	auto a14 = a13 + a6;
	auto a15 = a13 + a7;
	auto a16 = a14 + 2.0 * a6 + a7;
	auto a17 = a3;
	auto a18 = a13 - a8;
	auto a19 = a14 - a9;
	auto a20 = a15 - a10;
	auto a21 = a16 - a11;
	auto a22 = a17 - a12;
	auto a23 = a19 + 2.0 * (x[2] + x[1] + x[0]) + x[6];
	auto b0 = (1.0/2.0) * (-1.0 * h[6] - 2.0 * h[5] + 3.0 * h[4] - 1.0 * h[3] - 2.0 * h[2] + 1.0 * h[1] + 2.0 * h[0]);
	auto b1 = (1.0/14.0) * (10.0 * h[6] + 3.0 * h[5] - 11.0 * h[4] + 10.0 * h[3] + 3.0 * h[2] - 11.0 * h[1] - 4.0 * h[0]);
	auto b2 = (1.0/6.0) * (-2.0 * h[6] + 3.0 * h[5] - 1.0 * h[4] - 2.0 * h[3] + 3.0 * h[2] - 1.0 * h[1] + 0.0 * h[0]);
	auto b3 = (1.0/6.0) * (-1.0 * h[6] + 0.0 * h[5] + 1.0 * h[4] - 1.0 * h[3] + 0.0 * h[2] + 1.0 * h[1] + 0.0 * h[0]);
	auto b4 = (2.0 * h[6] - 1.0 * h[5] - 2.0 * h[4] + 3.0 * h[3] - 1.0 * h[2] - 2.0 * h[1] + 1.0 * h[0]);
	auto b5 = (1.0/2.0) * (-2.0 * h[6] + 1.0 * h[5] + 2.0 * h[4] - 1.0 * h[3] - 2.0 * h[2] + 3.0 * h[1] - 1.0 * h[0]);
	auto b6 = (1.0/14.0) * (3.0 * h[6] - 11.0 * h[5] - 4.0 * h[4] + 10.0 * h[3] + 3.0 * h[2] - 11.0 * h[1] + 10.0 * h[0]);

	auto b7 = (1.0/6.0) * (3.0 * h[6] - 1.0 * h[5] - 0.0 * h[4] - 2.0 * h[3] + 3.0 * h[2] - 1.0 * h[1] - 2.0 * h[0]);
	auto b8 = (1.0/6.0) * (h[5] - h[3] + h[1] - h[0]);
	auto b9 = (-1.0 * h[6] - 2.0 * h[5] + 1.0 * h[4] + 2.0 * h[3] - 1.0 * h[2] - 2.0 * h[1] + 3.0 * h[0]);
	auto b10 = (1.0/2.0) * (2.0 * h[4] - h[3] - 2.0 * h[2]+h[1]);
	auto b11 = (1.0/14.0)*(-2.0 * h[6] - 2.0 * h[5] - 2.0 * h[4] + 12.0 * h[3] + 5.0 * h[2] - 9.0 * h[1] - 2.0 * h[0]);
	auto b12 = (1.0/6.0) * (-2.0 * h[3] + 3.0 * h[2] - h[1]);
	auto b13 = (1.0/6.0) * (h[1] - h[3]);
	auto b14 = 2.0 * h[3] - h[2] - 2.0 * h[1] + h[0];
	auto b15 = (1.0/7.0) * (h[0] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6]);
	auto m0 = a8 * b0;
	auto m1 = a9 * b1;
	auto m2 = a10 * b2;
	auto m3 = a11 * b3;
	auto m4 = a12 * b4;
	auto m5 = a13 * b5;
	auto m6 = a14 * b6;
	auto m7 = a15 * b7;
	auto m8 = a16 * b8;
	auto m9 = a17 * b9;
	auto m10 = a18 * b10;
	auto m11 = a19 * b11;
	auto m12 = a20 * b12;
	auto m13 = a21 * b13;
	auto m14 = a22 * b14;
	auto m15 = a23 * b15;
	auto u0 = m0 + m10;
	auto u1 = m1 + m11;
	auto u2 = m2 + m12;
	auto u3 = m3 + m13;
	auto u4 = m4 + m14;
	auto u5 = m5 - m10;
	auto u6 = m6 - m11;
	auto u7 = m7 - m12;
	auto u8 = m8 - m13;
	auto u9 = m9 - m14;
	auto u10 = u1 + u3;
	auto u11 = u10 + u2;
	auto u12 = u0 + u11;
	auto u13 = u10 + u3;
	auto u14 = u13 - u2;
	auto u15 = u13 + 2.0 * u3 + u4 + u2;
	auto u16 = -u12 - 2.0 * (u13 + u3) - u4;
	auto u17 = u6 + u8;
	auto u18 = u17 + u7;
	auto u19 = u5 + u18;
	auto u20 = u17 + u8;
	auto u21 = u20 - u7;
	auto u22 = u20 + 2.0 * u8 + u9 + u7;
	auto u23 = -u19 - 2.0 * (u20 + u8) - u9;
	y[0] = u12 + m15;
	y[1] = u16 + u23 + m15;
	y[2] = u22 + m15;
	y[3] = u21 + m15;
	y[4] = u19 + m15;
	y[5] = u15 + m15;
	y[6] = u14 + m15;
	return y;
}

std::vector<std::vector<cmplx>> convolve_fast(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h) {
	int N1 = x.size();
	if( N1 == 2 ) {
		return convolve_fast_2(x,h);
	}
	if( N1 == 3 ) {
		return convolve_fast_3(x,h);
	}
	if( N1 == 4 ) {
		return convolve_fast_4(x,h);
	}
	if( N1 == 5 ) {
		return convolve_fast_5(x,h);
	}
	if( N1 == 7 ) {
		return convolve_fast_7(x,h);
	}
	assert(false);
	return std::vector<std::vector<cmplx>>();
}

bool can_fast_cyclic(int N) {
	switch (N) {
	case 2:
	case 3:
	case 4:
	case 5:
	case 7:
		return true;
	default:
		return false;
	};
}

std::vector<cmplx> convolve_fast(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	if (N == 1) {
		return std::vector<cmplx>(1, x[0] * h[0]);
	}
	auto pfacs = prime_factorization(N);
	std::vector<cmplx> y(N);
	if (can_fast_cyclic(N)) {
		std::vector<std::vector<cmplx>> X(N);
		std::vector<std::vector<std::complex<double>>>H(N);
		for (int n = 0; n < N; n++) {
			X[n].push_back(x[n]);
			H[n].push_back(h[n]);
		}
		std::vector<std::vector<cmplx>> Y;
		Y = convolve_fast(X, H);
		for (int n = 0; n < N; n++) {
			y[n] = Y[n][0];
		}
	} else if (pfacs.size() == 1) {
		std::vector<std::vector<cmplx>> X(N);
		std::vector<std::vector<std::complex<double>>>H(N);
		for (int n = 0; n < N; n++) {
			X[n].push_back(x[n]);
			H[n].push_back(h[n]);
		}
		auto ya = convolve_aperiodic(X, H);
		for (int n = 0; n < N - 1; n++) {
			y[n] = ya[n][0] + ya[n + N][0];
		}
		y[N - 1] = ya[N - 1][0];
	} else {
		int N2 = 1;
		auto p = prime_factorization(N);
		std::vector<std::pair<int, int>> pfacs(p.begin(), p.end());
		int N1 = std::pow(pfacs[0].first, pfacs[0].second);
		for (int n = 1; n < pfacs.size(); n++) {
			N2 *= std::pow(pfacs[n].first, pfacs[n].second);
		}
		//	fprintf( stderr, "%i = %i x %i\n", N, N1, N2);
		std::vector<std::vector<cmplx>> X(N1, std::vector<cmplx>(N2));
		std::vector<std::vector<std::complex<double>>>H(N1, std::vector<std::complex<double>>(N2));
		for (int n1 = 0; n1 < N1; n1++) {
			X[n1].resize(N2);
			H[n1].resize(N2);
			for (int n2 = 0; n2 < N2; n2++) {
				const int n = (n1 * N2 + n2 * N1) % N;
				X[n1][n2] = x[n];
				H[n1][n2] = h[n];
			}
		}
		auto Y = convolve_fast(X, H);
		for (int n1 = 0; n1 < N1; n1++) {
			for (int n2 = 0; n2 < N2; n2++) {
				const int n = (n1 * N2 + n2 * N1) % N;
				y[n] = Y[n1][n2];
			}
		}
	}
	for (int n = 0; n < N; n++) {
		y[n].set_goal();
	}
	return y;
}

std::vector<cmplx> convolve_N2(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	std::vector<cmplx> y(N, cmplx( { 0.0, 0.0 }));
	for (int n = 0; n < N; n++) {
		for (int m = 0; m < N; m++) {
			y[n] += x[m] * h[(n - m + N) % N];
		}
	}
	return y;
}

std::vector<cmplx> convolve(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	if (N == 1) {
		return std::vector<cmplx>(1, x[0] * h[0]);
	}
//	return convolve_fast(x, h);
	int fast_cnt;
	int fft_cnt;
	std::vector<cmplx> y;
	fast_cnt = math_vertex::operation_count(x * h).total();
	fft_cnt = math_vertex::operation_count(convolve_fft(x, h)).total();
	fprintf(stderr, "%i : %i %i\n", N, fft_cnt, fast_cnt);
	if (fast_cnt < fft_cnt) {
		return convolve_fast(x, h);
	} else {
		return convolve_fft(x, h);
	}
}

std::vector<cmplx> operator*(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	if (N == 1) {
		return std::vector<cmplx>(1, x[0] * h[0]);
	}
	int fast_cnt;
	int naive_cnt;
	std::vector<cmplx> y;
	fast_cnt = math_vertex::operation_count(convolve_fast(x, h)).total();
	naive_cnt = math_vertex::operation_count(convolve_N2(x, h)).total();
	if (naive_cnt < fast_cnt) {
		return convolve_N2(x, h);
	} else {
		return convolve_fast(x, h);
	}
}

std::vector<std::vector<cmplx>> convolve_aperiodic(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h)
{
	std::vector<std::vector<cmplx>> y;
	int N = x.size();
	if (N == 1) {
		return std::vector<std::vector<cmplx>>(1, convolve_fast(x[0], h[0]));
	}
	if (N % 2 != 0 && N % 3 != 0) {
		x.push_back(std::vector<cmplx>(x[0].size(), cmplx( {0.0, 0.0})));
		h.push_back(std::vector<std::complex<double>>(x[0].size(), {0.0, 0.0}));
	}
	int M = x.size();
	if (M % 3 == 0) {
		y = convolve_radix3(x, h);
	} else {
		y = convolve_radix2(x, h);
	}
	if (M != N) {
		y.pop_back();
		y.pop_back();
	}
	return y;
}

std::vector<cmplx> convolve_fft(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	auto X = fft(x, N, 0);
	for (int n = 0; n < N; n++) {
		h[n] *= (1.0 / N);
	}
	auto H = h;
	fftw(H);
	std::vector<cmplx> Y(N);
	for (int n = 0; n < N; n++) {
		Y[n] = X[n] * H[n];
		std::swap(Y[n].x, Y[n].y);
	}
	auto y = fft(Y, N, 0);
	for (int n = 0; n < N; n++) {
		std::swap(y[n].x, y[n].y);
	}
	return std::move(y);
}

std::vector<cmplx> convolve_karatsuba2(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	if (N == 1) {
		return std::vector<cmplx>(1, x[0] * h[0]);
	}
	std::vector<cmplx> y(N, cmplx( { 0.0, 0.0 }));
	std::vector<cmplx> xlo(N / 2), xhi(N / 2), xmd(N / 2);
	std::vector<std::complex<double>> hlo(N / 2), hhi(N / 2), hmd(N / 2);
	for (int n = 0; n < N / 2; n++) {
		xlo[n] = x[2 * n];
		xhi[n] = x[2 * n + 1];
		hlo[n] = h[2 * n];
		hhi[n] = h[2 * n + 1];
		xmd[n] = x[2 * n] + x[2 * n + 1];
		hmd[n] = h[2 * n] + h[(2 * n + 1 + N) % N];
	}
	auto ylo = convolve_karatsuba2(xlo, hlo);
	auto ymd = convolve_karatsuba2(xmd, hmd);
	auto yhi = convolve_karatsuba2(xhi, hhi);
	for (int n = 0; n < N / 2; n++) {
		ymd[n] -= ylo[n] + yhi[n];
	}
	for (int n = 0; n < N / 2; n++) {
		y[2 * n] = ylo[n] + yhi[(n - 1 + (N / 2)) % (N / 2)];
		y[2 * n + 1] = ymd[n];
	}
	for (int n = 0; n < N; n++) {
		y[n].set_goal();
	}
	return y;
}

