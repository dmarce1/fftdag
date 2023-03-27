#include "fft.hpp"
#include "util.hpp"

std::vector<cmplx> convolve_karatsuba(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> convolve_fast(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> convolve_fft(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> convolve_aperiodic(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> operator*(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> operator*(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> convolve_toomcook(std::vector<cmplx> x, std::vector<std::complex<double>> h);
bool can_agarwal_cooley(int N);
std::vector<cmplx> convolve_agarwal_cooley(std::vector<cmplx> x, std::vector<std::complex<double>> h);
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

std::vector<cmplx> convolve_karatsuba(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	std::vector<cmplx> y(2 * N - 1, cmplx( { 0.0, 0.0 }));
	std::vector<cmplx> xlo(N / 2), xhi(N / 2), xmd(N / 2);
	std::vector<std::complex<double>> hlo(N / 2), hhi(N / 2), hmd(N / 2);
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
		ymd[n] -= ylo[n] + yhi[n];
	}
	for (int n = 0; n < N - 1; n++) {
		y[n] += ylo[n];
		y[n + N / 2] += ymd[n];
		y[n + N] += yhi[n];
	}
	for (int n = 0; n < 2 * N - 1; n++) {
		y[n].set_goal();
	}
	return y;
}

std::vector<cmplx> convolve_toomcook(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	std::vector<cmplx> y(2 * N - 1, cmplx( { 0.0, 0.0 }));
	std::vector<cmplx> a2(N / 3), a3(N / 3), a4(N / 3), a5(N / 3), a6(N / 3);
	std::vector<std::complex<double>> b0(N / 3), b1(N / 3), b2(N / 3), b3(N / 3), b4(N / 3);
	std::vector<cmplx> u1(2 * N / 3 - 1), u2(2 * N / 3 - 1), u3(2 * N / 3 - 1), u4(2 * N / 3 - 1), u5(2 * N / 3 - 1);
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
		b1[n] = (h0 + h1 + h2) * 0.5;
		b2[n] = (h0 - h1 + h2) * (1.0 / 6.0);
		b3[n] = (h0 + 2.0 * h1 + 4.0 * h2) * (1.0 / 6.0);
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
	for (int n = 0; n < 2 * N / 3 - 1; n++) {
		y[n + 0 * N / 3] += u2[n];
		y[n + 1 * N / 3] += u1[n] - u3[n] + u4[n];
		y[n + 2 * N / 3] += -u2[n] + u3[n] + u5[n] - m4[n];
		y[n + 3 * N / 3] += -u4[n] - u5[n];
		y[n + 4 * N / 3] += m4[n];

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
	assert(false);
	return std::vector<std::vector<cmplx>>();
}

bool can_fast_cyclic(int N) {
	switch (N) {
	case 2:
	case 3:
	case 4:
	case 5:
		return true;
	default:
		return false;
	};
}

bool can_agarwal_cooley(int N) {
	int N2 = 1;
	auto p = prime_factorization(N);
	std::vector<std::pair<int, int>> pfacs(p.begin(), p.end());
	if (pfacs.size() <= 1) {
		return false;
	}
	int N1 = std::pow(pfacs[0].first, pfacs[0].second);
	for (int n = 1; n < pfacs.size(); n++) {
		N2 *= std::pow(pfacs[n].first, pfacs[n].second);
	}
	if (can_fast_cyclic(N1)) {
		if (can_fast_cyclic(N2)) {
			return true;
		} else {
			return can_agarwal_cooley(N2);
		}
	}
	return false;
}

std::vector<cmplx> convolve_agarwal_cooley(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	int N2 = 1;
	auto p = prime_factorization(N);
	std::vector<std::pair<int, int>> pfacs(p.begin(), p.end());
	int N1 = std::pow(pfacs[0].first, pfacs[0].second);
	for (int n = 1; n < pfacs.size(); n++) {
		N2 *= std::pow(pfacs[n].first, pfacs[n].second);
	}
	std::vector<cmplx> y(N);
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
	return std::move(y);
}

std::vector<cmplx> convolve_fast(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	if (N == 1) {
		return std::vector<cmplx>(1, x[0] * h[0]);
	}
	std::vector<std::vector<cmplx>> X(N);
	std::vector<std::vector<std::complex<double>>>H(N);
	std::vector<cmplx> y(N);
	{
		for (int n = 0; n < N; n++) {
			X[n].push_back(x[n]);
			H[n].push_back(h[n]);
		}
		std::vector<std::vector<cmplx>> Y;
		if (can_fast_cyclic(N)) {
			Y = convolve_fast(X, H);
			for (int n = 0; n < N; n++) {
				y[n] = Y[n][0];
			}
			return y;
		}
	}
	if (can_agarwal_cooley(N)) {
		return convolve_agarwal_cooley(x, h);
	} else {
		/*std::vector<cmplx> X(N);
		 std::vector<std::complex<double>> H(N);
		 int M = 2 * N - 1;
		 while (!can_agarwal_cooley(M)) {
		 M++;
		 }
		 X.resize(M);
		 H.resize(M);
		 for (int n = 0; n < N; n++) {
		 X[n] = x[n];
		 }
		 for (int n = N; n < M; n++) {
		 X[n] = cmplx( { 0.0, 0.0 });
		 }
		 H[0] = h[0];
		 for (int n = 1; n < N - 1; n++) {
		 H[M + n + 1 - N] = H[n] = h[n];
		 }
		 fprintf(stderr, "c %i %i\n", N, M);
		 y = convolve_fast(X, H);
		 y.resize(N);
		 return y;*/
		std::vector<cmplx> yc(N);
		auto ya = convolve_aperiodic(x, h);
		for (int n = 0; n < N; n++) {
			yc.push_back(ya[n]);
		}
		for (int n = 0; n < N - 1; n++) {
			yc[n] = ya[n] + ya[n + N];
		}
		yc[N - 1] = ya[N - 1];
		return yc;
	}
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
	int fast_cnt;
	int fft_cnt;
	std::vector<cmplx> y;
	fast_cnt = math_vertex::operation_count(x * h).total();
	fft_cnt = math_vertex::operation_count(convolve_fft(x, h)).total();
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
//	fprintf(stderr, "%i : %i %i %i\n", N, naive_cnt, fft_cnt, fast_cnt);
	if (naive_cnt < fast_cnt) {
		return convolve_N2(x, h);
	} else {
		return convolve_fast(x, h);
	}
}

std::vector<cmplx> convolve_aperiodic(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	std::vector<cmplx> y;
	int N = x.size();
	if (N == 1) {
		return std::vector<cmplx>(1, x[0] * h[0]);
	}
	if (N % 2 != 0 && N % 3 != 0) {
		x.push_back(cmplx( { 0.0, 0.0 }));
		h.push_back(std::complex<double>( { 0.0, 0.0 }));
	}
	int M = x.size();
	if (M % 2 == 0) {
		y = convolve_karatsuba(x, h);
	} else {
		y = convolve_toomcook(x, h);
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

