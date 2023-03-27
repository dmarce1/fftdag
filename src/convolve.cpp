#include "fft.hpp"
#include "util.hpp"

std::vector<cmplx> convolve_karatsuba(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> convolve_fast(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> convolve_fft(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> operator&(std::vector<cmplx> x, std::vector<std::complex<double>> h);
std::vector<cmplx> operator*(std::vector<cmplx> x, std::vector<std::complex<double>> h);
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
	auto a8 = a1 - a5;
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
			for (int n = 0; n < N; n++) {
				Y = convolve_fast(X, H);
				y[n] = Y[n][0];
			}
			return y;
		}
	}
	if (can_agarwal_cooley(N)) {
//		fprintf( stderr, "agarwal_cooley: %i\n", N);
		return convolve_agarwal_cooley(x, h);
	}
	std::vector<cmplx> yc(N);
	auto ya = x & h;
	for (int n = 0; n < N; n++) {
		yc.push_back(ya[n]);
	}
	for (int n = 0; n < N - 1; n++) {
		yc[n] = ya[n] + ya[n + N];
	}
	yc[N - 1] = ya[N - 1];
	return yc;
}

std::vector<cmplx> operator*(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	if (N == 1) {
		return std::vector<cmplx>(1, x[0] * h[0]);
	}
	int fast_cnt;
	int fft_cnt;
	std::vector<cmplx> y;
	fast_cnt = math_vertex::operation_count(convolve_fast(x, h)).total();
	fft_cnt = math_vertex::operation_count(convolve_fft(x, h)).total();
//	fprintf(stderr, "%i %i %i\n", N, fft_cnt, fast_cnt);
	if (fast_cnt < fft_cnt) {
		return convolve_fast(x, h);
	} else {
		return convolve_fft(x, h);
	}
}

std::vector<cmplx> operator&(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	std::vector<cmplx> y;
	int N = x.size();
	if (N == 1) {
		return std::vector<cmplx>(1, x[0] * h[0]);
	}
	if (N % 2 != 0) {
		x.push_back(cmplx( { 0.0, 0.0 }));
		h.push_back(std::complex<double>( { 0.0, 0.0 }));
	}
	y = convolve_karatsuba(x, h);
	if (N % 2 != 0) {
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
	auto ylo = xlo & hlo;
	auto ymd = xmd & hmd;
	auto yhi = xhi & hhi;
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

/*



 bool can_fast_convolve(int N) {
 if (N == 2) {
 return true;
 } else if (N == 3) {
 return true;
 } else if (N == 4) {
 return true;
 } else {
 auto pfacs = prime_factorization(N);
 if (pfacs.size() > 1) {
 for (int n = 0; n < pfacs.size(); n++) {
 auto k = std::pow(pfacs[n].first, pfacs[n].second);
 if (!can_fast_convolve(k)) {
 return false;
 }
 }
 return true;
 }
 }
 return false;
 }

 std::vector<cmplx> fast_convolve(std::vector<cmplx> x, std::vector<cmplx> h, int N) {
 assert(can_fast_convolve(N));
 if (N == 2) {
 return convolve2(x, h);
 } else if (N == 3) {
 return convolve3(x, h);
 } else if (N == 4) {
 return convolve3(x, h);
 } else {
 auto pfacs = prime_factorization(N);
 auto N1 = std::pow(pfacs[0].first, pfacs[0].second);
 int N2 = 1;
 for (int n = 1; n < pfacs.size(); n++) {
 N2 *= std::pow(pfacs[n].first, pfacs[n].second);
 }
 std::vector<std::vector<cmplx>> y1(N1, std::vector<cmplx>(N2));
 std::vector<std::vector<cmplx>> y2(N2, std::vector<cmplx>(N1));
 std::vector<std::vector<cmplx>> x1(N1, std::vector<cmplx>(N2));
 std::vector<std::vector<cmplx>> x2(N2, std::vector<cmplx>(N1));
 std::vector<std::vector<cmplx>> h1(N1, std::vector<cmplx>(N2));
 std::vector<std::vector<cmplx>> h2(N2, std::vector<cmplx>(N1));
 for (int n1 = 0; n1 < N1; n1++) {
 for (int n2 = 0; n2 < N2; n2++) {
 int n = n1 % N2 + n2 % N1;
 x1[n1][n2] = x[n];
 h1[n1][n2] = h[n];
 }
 }
 for (int n1 = 0; n1 < N1; n1++) {
 y[n1] = fast_convolve(x1[n1], h1[n1], N1);
 }
 for (int n1 = 0; n1 < N1; n1++) {
 for (int n2 = 0; n2 < N2; n2++) {
 }
 }
 }
 }*/
