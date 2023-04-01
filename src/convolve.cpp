#include "fft.hpp"
#include "util.hpp"
#include "convolve.hpp"

#define FFT_CONVOLVE 0
#define FAST_CONVOLVE 1
#define SLOW_CONVOLVE 2

static std::array<std::array<double, 2>, 3> toom2a = { { { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } } };

static std::array<std::array<double, 3>, 3> toom2b = { { { 1.0, 0.0, 0.0 }, { -1.0, 1.0, -1.0 }, { 0.0, 0.0, 1.0 } } };

static std::array<std::array<double, 3>, 5> toom3a = { { { 1.0, 0.0, 0.0 }, { 1.0, 1 / 2.0, 1 / 4.0 }, { 1.0, -1.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } } };

static std::array<std::array<double, 5>, 5> toom3b = { { { 1.0, 0.0, 0.0, 0.0, 0.0 }, { -2.0, 8 / 3.0, -1 / 6.0, -1 / 2.0, 1 / 2.0 }, { -1.0, 0.0, 1 / 2.0, 1 / 2.0, -1.0 }, { 2.0, -8 / 3.0, -1 / 3.0, 1.0, -1 / 2.0 }, { 0.0, 0.0, 0.0, 0.0, 1.0 } } };

static std::array<std::array<double, 4>, 7> toom4a = { { { 1.0, 0.0, 0.0, 0.0 }, { 1.0, 2.0, 4.0, 8.0 }, { 1.0, -2.0, 4.0, -8.0 }, { 1.0, -1 / 4.0, 1 / 16.0, -1 / 64.0 }, { 1.0, -1.0, 1.0, -1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 0.0, 1.0 } } };

static std::array<std::array<double, 7>, 7> toom4b = { { { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 4.0, -1 / 108.0, -1 / 84.0, -4096 / 945.0, 2 / 9.0, 2 / 15.0, 1.0 }, { -5 / 4.0, -1 / 24.0, -1 / 24.0, 0.0, 2 / 3.0, 2 / 3.0, 4.0 }, { -5.0, -1 / 108.0, 1
		/ 28.0, 1024 / 189.0, -17 / 18.0, 1 / 2.0, -5 / 4.0 }, { 1 / 4.0, 1 / 24.0, 1 / 24.0, 0.0, -1 / 6.0, -1 / 6.0, -5.0 }, { 1.0, 1 / 54.0, -1 / 42.0, -1024 / 945.0, 2 / 9.0, -2 / 15.0, 1 / 4.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } } };

static std::array<std::array<double, 5>, 9> toom5a = { { { 1.0, 0.0, 0.0, 0.0, 0.0 }, { 1.0, 3 / 2.0, 9 / 4.0, 27 / 8.0, 81 / 16.0 }, { 1.0, 2 / 3.0, 4 / 9.0, 8 / 27.0, 16 / 81.0 }, { 1.0, -1.0, 1.0, -1.0, 1.0 }, { 1.0, 1.0, 1.0, 1.0, 1.0 }, { 1.0, -2
		/ 3.0, 4 / 9.0, -8 / 27.0, 16 / 81.0 }, { 1.0, -3 / 2.0, 9 / 4.0, -27 / 8.0, 81 / 16.0 }, { 1.0, 1 / 2.0, 1 / 4.0, 1 / 8.0, 1 / 16.0 }, { 0.0, 0.0, 0.0, 0.0, 1.0 } } };

static std::array<std::array<double, 9>, 9> toom5b = { { { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { -2.0, -32 / 975.0, -6561 / 1300.0, 6 / 25.0, 18 / 25.0, -6561 / 9100.0, -16 / 975.0, 48 / 7.0, 1 / 2.0 }, { -133 / 36.0, 128 / 2925.0, 6561
		/ 2600.0, -18 / 25.0, -18 / 25.0, 6561 / 2600.0, 128 / 2925.0, 0.0, -1.0 }, { 133 / 18.0, 88 / 585.0, 1458 / 65.0, -1 / 6.0, -169 / 50.0, -729 / 650.0, 92 / 2925.0, -76 / 3.0, -133 / 72.0 }, { 133 / 36.0, -32 / 225.0, -729 / 200.0, 97 / 50.0, 97
		/ 50.0, -729 / 200.0, -32 / 225.0, 0.0, 133 / 36.0 }, { -133 / 18.0, -632 / 2925.0, -31347 / 1300.0, -79 / 75.0, 23 / 5.0, 729 / 260.0, 4 / 117.0, 76 / 3.0, 133 / 72.0 }, { -1.0, 32 / 325.0, 729 / 650.0, -18 / 25.0, -18 / 25.0, 729 / 650.0, 32
		/ 325.0, 0.0, -133 / 36.0 }, { 2.0, 32 / 325.0, 2187 / 325.0, 12 / 25.0, -36 / 25.0, -2187 / 2275.0, -16 / 325.0, -48 / 7.0, -1 / 2.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } } };

static std::array<std::array<double, 6>, 11> toom6a = { { { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 1.0, 1 / 3.0, 1 / 9.0, 1 / 27.0, 1 / 81.0, 1 / 243.0 }, { 1.0, -1 / 3.0, 1 / 9.0, -1 / 27.0, 1 / 81.0, -1 / 243.0 }, { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }, { 1.0,
		-1.0, 1.0, -1.0, 1.0, -1.0 }, { 1.0, -1 / 2.0, 1 / 4.0, -1 / 8.0, 1 / 16.0, -1 / 32.0 }, { 1.0, -3 / 2.0, 9 / 4.0, -27 / 8.0, 81 / 16.0, -243 / 32.0 }, { 1.0, 1 / 4.0, 1 / 16.0, 1 / 64.0, 1 / 256.0, 1 / 1024.0 }, { 1.0, -2 / 3.0, 4 / 9.0, -8 / 27.0,
		16 / 81.0, -32 / 243.0 }, { 1.0, -2.0, 4.0, -8.0, 16.0, -32.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } } };

static std::array<std::array<double, 11>, 11> toom6b = { { { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 2 / 3.0, -2187 / 1540.0, -6561 / 980.0, 1 / 900.0, 3 / 20.0, 256 / 45.0, -256 / 40425.0, 262144 / 72765.0, -2187 / 1100.0, 1 / 3780.0,
		-1 / 36.0 }, { -85 / 4.0, -729 / 140.0, 2187 / 140.0, 1 / 540.0, -1 / 20.0, -1024 / 135.0, 0.0, 524288 / 31185.0, 729 / 440.0, 1 / 22680.0, -1 / 54.0 }, { -95 / 3.0, 89667 / 6160.0, 373977 / 3920.0, -47 / 2160.0, -251 / 80.0, -14272 / 135.0, 1088
		/ 8085.0, -2031616 / 218295.0, 2187 / 55.0, -16 / 2835.0, 85 / 144.0 }, { 623 / 6.0, 78003 / 880.0, -41553 / 560.0, -41 / 720.0, -129 / 80.0, 1408 / 45.0, 128 / 1155.0, -524288 / 3465.0, 729 / 220.0, -1 / 180.0, 95 / 108.0 }, { 252.0, 9477 / 80.0,
		-264627 / 560.0, 631 / 10800.0, 275 / 16.0, 71296 / 135.0, -128 / 175.0, -131072 / 567.0, -21141 / 100.0, 49 / 1620.0, -623 / 216.0 }, { 905 / 12.0, -2187 / 880.0, -150903 / 560.0, 731 / 2160.0, 1649 / 80.0, 50944 / 135.0, -256 / 231.0, -524288
		/ 31185.0, -80919 / 440.0, 167 / 3240.0, -7.0 }, { -185.0, -8019 / 70.0, 74358 / 245.0, 19 / 45.0, -93 / 10.0, -14656 / 45.0, 64 / 245.0, 14876672 / 72765.0, 13851 / 110.0, -11 / 1890.0, -905 / 432.0 }, { -159.0, -124659 / 1540.0, 6561 / 20.0, 13
		/ 60.0, -369 / 20.0, -6016 / 15.0, 384 / 385.0, 524288 / 3465.0, 19683 / 110.0, -29 / 630.0, 185 / 36.0 },
		{ -36.0, -6561 / 385.0, 19683 / 245.0, 1 / 25.0, -27 / 5.0, -512 / 5.0, 4608 / 13475.0, 262144 / 8085.0, 13122 / 275.0, -2 / 105.0, 53 / 12.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } }

};

inline double toomA(int N, int n, int m) {
	if (N == 2) {
		return toom2a[n][m];
	} else if (N == 3) {
		return toom3a[n][m];
	} else if (N == 4) {
		return toom4a[n][m];
	} else if (N == 5) {
		return toom5a[n][m];
	} else if (N == 6) {
		return toom6a[n][m];
	}
	return 0.0 / 0.0;
}

inline double toomB(int N, int n, int m) {
	if (N == 2) {
		return toom2b[n][m];
	} else if (N == 3) {
		return toom3b[n][m];
	} else if (N == 4) {
		return toom4b[n][m];
	} else if (N == 5) {
		return toom5b[n][m];
	} else if (N == 6) {
		return toom6b[n][m];
	}
	return 0.0 / 0.0;
}

std::vector<cmplx> operator*(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	return convolve_fast(x, h);
}

template<class T>
std::vector<T> operator+(std::vector<T> x1, std::vector<T> x2) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] = x1[n] + x2[n];
	}
	return x1;
}

template<class T>
std::vector<T> operator-(std::vector<T> x1, std::vector<T> x2) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] = x1[n] - x2[n];
	}
	return x1;
}

template<class T>
std::vector<T> operator-(std::vector<T> x1) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] = -x1[n];
	}
	return x1;
}

template<class T>
std::vector<T> operator*(double a, std::vector<T> x1) {
	for (int n = 0; n < x1.size(); n++) {
		x1[n] = x1[n] * a;
	}
	return x1;
}

template<class T>
std::vector<T> operator*(std::vector<T> x1, double a) {
	return a * x1;
}

bool can_fast_cyclic(int N) {
	switch (N) {
	case 2:
	case 3:
	case 4:
	case 5:
	case 7:
	case 8:
	case 9:
		return true;
	default:
		return false;
	};
}

template<class T, class V>
std::vector<T> convolve_fast_2(std::vector<T> x, std::vector<V> h) {
	std::vector<T> y = x;
	T a0 = x[0] + x[1];
	T a1 = x[0] - x[1];
	V b0 = 0.5 * (h[0] + h[1]);
	V b1 = 0.5 * (h[0] - h[1]);
	T m0 = a0 * b0;
	T m1 = a1 * b1;
	y[0] = m0 + m1;
	y[1] = m0 - m1;
	return std::move(y);
}

template<class T, class V>
std::vector<T> convolve_fast_3(std::vector<T> x, std::vector<V> h) {
	std::vector<T> y = x;
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

template<class T, class V>
std::vector<T> convolve_fast_4(std::vector<T> x, std::vector<V> h) {
	std::vector<T> y = x;
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

template<class T, class V>
std::vector<T> convolve_fast_5(std::vector<T> x, std::vector<V> h) {
	std::vector<T> y = x;
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
	auto b2 = 0.2 * (-2.0 * h[0] - 2.0 * h[1] + 3.0 * h[2] - 2.0 * h[3] + 3.0 * h[4]);
	auto b3 = h[1] - h[0] - h[2] + h[3];
	auto b4 = h[1] - h[0] - h[2] + h[4];
	auto b5 = 0.2 * (3.0 * h[0] - 2.0 * h[1] + 3.0 * h[2] - 2.0 * h[3] - 2.0 * h[4]);
	auto b6 = h[3] - h[2];
	auto b7 = h[1] - h[2];
	auto b8 = 0.2 * (4.0 * h[2] - h[3] - h[4] - h[0] - h[1]);
	auto b9 = 0.2 * (h[0] + h[1] + h[2] + h[3] + h[4]);
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

template<class T, class V>
std::vector<T> convolve_fast_7(std::vector<T> x, std::vector<V> h) {
	std::vector<T> y = x;
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
	auto b0 = (1.0 / 2.0) * (-1.0 * h[6] - 2.0 * h[5] + 3.0 * h[4] - 1.0 * h[3] - 2.0 * h[2] + 1.0 * h[1] + 2.0 * h[0]);
	auto b1 = (1.0 / 14.0) * (10.0 * h[6] + 3.0 * h[5] - 11.0 * h[4] + 10.0 * h[3] + 3.0 * h[2] - 11.0 * h[1] - 4.0 * h[0]);
	auto b2 = (1.0 / 6.0) * (-2.0 * h[6] + 3.0 * h[5] - 1.0 * h[4] - 2.0 * h[3] + 3.0 * h[2] - 1.0 * h[1] + 0.0 * h[0]);
	auto b3 = (1.0 / 6.0) * (-1.0 * h[6] + 0.0 * h[5] + 1.0 * h[4] - 1.0 * h[3] + 0.0 * h[2] + 1.0 * h[1] + 0.0 * h[0]);
	auto b4 = (2.0 * h[6] - 1.0 * h[5] - 2.0 * h[4] + 3.0 * h[3] - 1.0 * h[2] - 2.0 * h[1] + 1.0 * h[0]);
	auto b5 = (1.0 / 2.0) * (-2.0 * h[6] + 1.0 * h[5] + 2.0 * h[4] - 1.0 * h[3] - 2.0 * h[2] + 3.0 * h[1] - 1.0 * h[0]);
	auto b6 = (1.0 / 14.0) * (3.0 * h[6] - 11.0 * h[5] - 4.0 * h[4] + 10.0 * h[3] + 3.0 * h[2] - 11.0 * h[1] + 10.0 * h[0]);
	auto b7 = (1.0 / 6.0) * (3.0 * h[6] - 1.0 * h[5] - 0.0 * h[4] - 2.0 * h[3] + 3.0 * h[2] - 1.0 * h[1] - 2.0 * h[0]);
	auto b8 = (1.0 / 6.0) * (h[5] - h[3] + h[1] - h[0]);
	auto b9 = (-1.0 * h[6] - 2.0 * h[5] + 1.0 * h[4] + 2.0 * h[3] - 1.0 * h[2] - 2.0 * h[1] + 3.0 * h[0]);
	auto b10 = (1.0 / 2.0) * (2.0 * h[4] - h[3] - 2.0 * h[2] + h[1]);
	auto b11 = (1.0 / 14.0) * (-2.0 * h[6] - 2.0 * h[5] - 2.0 * h[4] + 12.0 * h[3] + 5.0 * h[2] - 9.0 * h[1] - 2.0 * h[0]);
	auto b12 = (1.0 / 6.0) * (-2.0 * h[3] + 3.0 * h[2] - h[1]);
	auto b13 = (1.0 / 6.0) * (h[1] - h[3]);
	auto b14 = 2.0 * h[3] - h[2] - 2.0 * h[1] + h[0];
	auto b15 = (1.0 / 7.0) * (h[0] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6]);
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

template<class T, class V>
std::vector<T> convolve_fast_8(std::vector<T> x, std::vector<V> h) {
	std::vector<T> y = x;
	auto a0 = x[0] + x[4];
	auto a1 = x[1] + x[5];
	auto a2 = x[2] + x[6];
	auto a3 = x[3] + x[7];
	auto a4 = a0 + a2;
	auto a5 = a1 + a3;
	auto a6 = x[0] - x[4];
	auto a7 = x[1] - x[5];
	auto a8 = x[2] - x[6];
	auto a9 = x[3] - x[7];
	auto a10 = a0 - a2;
	auto a11 = a1 - a3;
	auto a12 = a4 + a5;
	auto a13 = a4 - a5;
	auto a14 = a7 + a9;
	auto a15 = a6 + a8;
	auto a16 = a15 - a14;
	auto a17 = a8 - a9;
	auto a18 = a6 - a7;
	auto a19 = a10 + a11;
	auto b0 = h[0] + h[4];
	auto b1 = h[1] + h[5];
	auto b2 = h[2] + h[6];
	auto b3 = h[3] + h[7];
	auto b4 = b0 + b2;
	auto b5 = b1 + b3;
	auto b6 = 0.5 * (((h[2] - h[6])) - (h[0] - h[4]) - ((h[1] - h[5]) - (h[3] - h[7])));
	auto b7 = 0.5 * (((h[2] - h[6])) - (h[0] - h[4]) + ((h[1] - h[5]) + (h[3] - h[7])));
	auto b8 = 0.5 * (((h[0] - h[4]) + (h[2] - h[6])) + ((h[1] - h[5]) + (h[3] - h[7])));
	auto b9 = 0.5 * (((h[0] - h[4]) + (h[2] - h[6])) + ((h[1] - h[5]) - (h[3] - h[7])));
	auto b10 = 0.25 * ((b1 - b3) - (b0 - b2));
	auto b11 = 0.25 * ((b0 - b2) + (b1 - b3));
	auto b12 = 0.125 * (b4 + b5);
	auto b13 = 0.125 * (b4 - b5);
	auto b14 = 0.5 * ((h[0] - h[4]) - (h[3] - h[7]));
	auto b15 = 0.5 * ((h[0] - h[4]) + (h[1] - h[5]));
	auto b16 = 0.5 * (h[0] - h[4]);
	auto b17 = 0.5 * ((h[0] - h[4]) + (h[2] - h[6]));
	auto b18 = 0.5 * ((h[2] - h[6]) - (h[0] - h[4]));
	auto b19 = 0.25 * (b0 - b2);
	auto m0 = a6 * b6;
	auto m1 = a7 * b7;
	auto m2 = a8 * b8;
	auto m3 = a9 * b9;
	auto m4 = a10 * b10;
	auto m5 = a11 * b11;
	auto m6 = a12 * b12;
	auto m7 = a13 * b13;
	auto m8 = a14 * b14;
	auto m9 = a15 * b15;
	auto m10 = a16 * b16;
	auto m11 = a17 * b17;
	auto m12 = a18 * b18;
	auto m13 = a19 * b19;
	auto u0 = m8 + m10;
	auto u1 = m9 - m10;
	auto u2 = m3 + m11;
	auto u3 = m11 - m2;
	auto u4 = m1 + m12;
	auto u5 = m0 - m12;
	auto u6 = m13 - m5;
	auto u7 = m13 + m4;
	auto u8 = m7 + m6;
	auto u9 = m6 - m7;
	auto u10 = u0 - u2;
	auto u11 = u6 + u8;
	auto u12 = u1 + u3;
	auto u13 = u7 + u9;
	auto u14 = u0 + u4;
	auto u15 = u8 - u6;
	auto u16 = u1 + u5;
	auto u17 = u9 - u7;
	y[0] = u10 + u11;
	y[1] = u12 + u13;
	y[2] = u14 + u15;
	y[3] = u16 + u17;
	y[4] = u11 - u10;
	y[5] = u13 - u12;
	y[6] = u15 - u14;
	y[7] = u17 - u16;
	return y;
}

template<class T, class V>
std::vector<T> convolve_fast_9(std::vector<T> x, std::vector<V> h) {
	std::vector<T> y = x;
	const T a0 = x[0] - x[6];
	const T a1 = x[1] - x[7];
	const T a2 = x[2] - x[8];
	const T a3 = x[3] - x[6];
	const T a4 = x[4] - x[7];
	const T a5 = x[5] - x[8];
	const T a6 = x[0] + x[3] + x[6];
	const T a7 = x[1] + x[4] + x[7];
	const T a8 = x[2] + x[5] + x[8];
	const T a9 = a0 + a2;
	const T a10 = a3 + a5;
	const T a11 = a6 + a7 + a8;
	const T a12 = a10 + a4;
	const T a13 = a9 + a1;
	const T a14 = a13 - a12;
	const T a15 = a10 - a4;
	const T a16 = a9 - a1;
	const T a17 = a16 - a15;
	const T a18 = a3;
	const T a19 = a0 - a3;
	const T a20 = a0;
	const T a21 = a5;
	const T a22 = a2 - a5;
	const T a23 = a2;
	const T a24 = a0 - a4 - a22;
	const T a25 = a19 + a5 - a1;
	const T a26 = a24 - a25;
	const T a27 = a6 - a8;
	const T a28 = a7 - a8;
	const T a29 = a27 + a28;
	const V b0 = 2.0 * h[6] - h[0] - h[3];
	const V b1 = 2.0 * h[7] - h[1] - h[4];
	const V b2 = 2.0 * h[8] - h[2] - h[5];
	const V b3 = h[0] - 2.0 * h[3] + h[6];
	const V b4 = h[1] - 2.0 * h[4] + h[7];
	const V b5 = h[2] - 2.0 * h[5] + h[8];
	const V b6 = h[0] + h[3] + h[6];
	const V b7 = h[1] + h[4] + h[7];
	const V b8 = h[2] + h[5] + h[8];
	const V b9 = (1.0 / 9.0) * (b6 + b7 + b8);
	const V b10 = (1.0 / 18.0) * (b0 + 3.0 * b1 + 2.0 * b2 - 2.0 * b3 - 3.0 * b4 - b5);
	const V b11 = (1.0 / 18.0) * (b0 - b2 + b3 + 3.0 * b4 + 2.0 * b5);
	const V b12 = b10 + b11;
	const V b13 = (1.0 / 6.0) * (b1 - b0 - b4 + b5);
	const V b14 = (1.0 / 6.0) * (b0 - b2 - b3 + b4);
	const V b15 = b13 + b14;
	const V b16 = (1.0 / 3.0) * (2.0 * b0 + b1 - b2 - 2.0 * b3 + b5);
	const V b17 = (1.0 / 3.0) * (2.0 * b0 - b2 + b4);
	const V b18 = b17 - b16;
	const V b19 = (1.0 / 3.0) * (b0 - b1 - 2.0 * b2 + b4);
	const V b20 = (1.0 / 3.0) * (b3 - b1 - 2.0 * b5);
	const V b21 = b20 - b19;
	const V b22 = (1.0 / 9.0) * (b0 - b2 - 2.0 * b3 + 2.0 * b5);
	const V b23 = (1.0 / 9.0) * (b2 - b0 - b3 + b5);
	const V b24 = b23 - b22;
	const V b25 = (1.0 / 3.0) * (b6 - b8);
	const V b26 = (1.0 / 3.0) * (b7 - b8);
	const V b27 = (1.0 / 3.0) * (b25 + b26);
	const T m0 = a11 * b9;
	const T m1 = a12 * b10;
	const T m2 = a13 * b11;
	const T m3 = a14 * b12;
	const T m4 = a15 * b13;
	const T m5 = a16 * b14;
	const T m6 = a17 * b15;
	const T m7 = a18 * b16;
	const T m8 = a19 * b17;
	const T m9 = a20 * b18;
	const T m10 = a21 * b19;
	const T m11 = a22 * b20;
	const T m12 = a23 * b21;
	const T m13 = a24 * b22;
	const T m14 = a25 * b23;
	const T m15 = a26 * b24;
	const T m16 = a27 * b25;
	const T m17 = a28 * b26;
	const T m18 = a29 * b27;
	const T u0 = m1 + m2;
	const T u1 = m4 + m5;
	const T u2 = m14 + m15;
	const T u3 = u0 + u1;
	const T u4 = m1 + m3;
	const T u5 = m4 + m6;
	const T u6 = m13 + m15;
	const T u7 = m7 - u3;
	const T u8 = u4 + u5;
	const T u9 = m10 - u6;
	const T u10 = m8 + u2 + u7;
	const T u11 = u8 + m11 + u9;
	const T u12 = u4 - u5 + u2;
	const T u13 = u7 + u8 + m9 + u6;
	const T u14 = u3 + m12 + u9 + u2;
	const T u15 = u0 - u1 + u6;
	const T u16 = m16 - m18;
	const T u17 = m17 - m18;
	const T u18 = m0 + u16;
	const T u19 = m0 - u16 - u17;
	const T u20 = m0 + u17;
	y[0] = u13 - u10 + u18;
	y[1] = u14 - u11 + u19;
	y[2] = u15 - u12 + u20;
	y[3] = u18 - u13;
	y[4] = u19 - u14;
	y[5] = u20 - u15;
	y[6] = u10 + u18;
	y[7] = u11 + u19;
	y[8] = u12 + u20;
	return y;
}

template<class T, class V>
std::vector<T> convolve_toom(int R, std::vector<T> x, std::vector<V> h) {
	assert(R < 7);
	int N = x.size();
	assert(N % R == 0);
	if (N == 1) {
		return std::vector<T>(1, x[0] * h[0]);
	}
	std::vector<T> y(N);
	std::vector<std::vector<T>> a(2 * R - 1, std::vector<T>(N / R));
	std::vector<std::vector<V>> b(2 * R - 1, std::vector<V>(N / R));
	std::vector<std::vector<T>> c(2 * R - 1, std::vector<T>(N / R));
	std::vector<std::vector<T>> m(2 * R - 1, std::vector<T>(N / R));
	for (int r = 0; r < 2 * R - 1; r++) {
		for (int n = 0; n < N / R; n++) {
			a[r][n] = toomA(R, r, 0) * x[R * n];
			b[r][n] = toomA(R, r, 0) * h[R * n];
			for (int i = 1; i < R; i++) {
				a[r][n] = a[r][n] + toomA(R, r, i) * x[R * n + i];
				b[r][n] = b[r][n] + toomA(R, r, i) * h[R * n + i];
			}
		}
	}
	for (int r = 0; r < 2 * R - 1; r++) {
		m[r] = convolve_dispatch(a[r], b[r]);
	}
	for (int n = 0; n < N / R; n++) {
		for (int r = 0; r < 2 * R - 1; r++) {
			c[r][n] = toomB(R, r, 0) * m[0][n];
			for (int i = 1; i < 2 * R - 1; i++) {
				c[r][n] = c[r][n] + toomB(R, r, i) * m[i][n];
			}
		}
	}
	for (int n = 0; n < N / R; n++) {
		y[R * n + R - 1] = c[R - 1][n];
		for (int r = 0; r < R - 1; r++) {
			y[R * n + r] = c[r][n] + c[R + r][mod(n - 1, N / R)];
		}
	}
	return y;
}

bool can_agarwal(int N) {
	auto pfacs = prime_factorization(N);
	if (pfacs.size() > 1) {
		for (auto fac : pfacs) {
			if (!can_fast_cyclic(std::pow(fac.first, fac.second))) {
				if (fac.first == 2) {
					continue;
				}
				if (fac.first == 3) {
					continue;
				}
				if (fac.first == 5) {
					continue;
				}
				return false;
			}
		}
		return true;
	}
	return false;
}

template<class T, class V>
std::vector<T> convolve_tiny(std::vector<T> x, std::vector<V> h);

std::vector<cmplx> convolve_dispatch(std::vector<cmplx> X, std::vector<std::complex<double>> H) {
	int N = X.size();
	std::vector<cmplx> y;
	if (can_fast_cyclic(N)) {
		y = convolve_tiny(X, H);
	} else if (N % 2 == 0) {
		y = convolve_toom(2, X, H);
	} else if (N % 3 == 0) {
		y = convolve_toom(3, X, H);
	} else if (N % 5 == 0) {
		y = convolve_toom(5, X, H);
	} else {
		y = convolve_fft(X, H);
	}
	for (auto z : y) {
		z.set_goal();
	}
	return y;
}

std::vector<std::vector<cmplx>> convolve_dispatch(std::vector<std::vector<cmplx>> X, std::vector<std::vector<std::complex<double>>>H) {
	int N = X.size();
	std::vector<std::vector<cmplx>> y;
	if (can_fast_cyclic(N)) {
		y = convolve_tiny(X, H);
	} else if (N % 2 == 0) {
		y = convolve_toom(2, X, H);
	} else if (N % 3 == 0) {
		y = convolve_toom(3, X, H);
	} else if (N % 5 == 0) {
		y = convolve_toom(5, X, H);
	} else {
		assert(false);
	}
	for (auto z : y) {
		for( auto q : z ) {
			q.set_goal();
		}
	}
	return y;
}

std::vector<std::vector<cmplx>> convolve(std::vector<std::vector<cmplx>> x, std::vector<std::vector<std::complex<double>>>h, int opts) {
	return convolve_dispatch(x, h);
}

template<class T, class V>
std::vector<T> convolve_agarwal(int N1, int N2, std::vector<T> x, std::vector<V> h) {
	assert(can_agarwal(N1 * N2));
	int N = N1 * N2;
	std::vector<T> y(N);
	std::vector<std::vector<T>> X(N1, std::vector<T>(N2));
	std::vector<std::vector<V>> H(N1, std::vector<V>(N2));
	std::vector<std::vector<T>> Y;
	for (int n1 = 0; n1 < N1; n1++) {
		X[n1].resize(N2);
		H[n1].resize(N2);
		for (int n2 = 0; n2 < N2; n2++) {
			const int n = (n1 * N2 + n2 * N1) % N;
			X[n1][n2] = x[n];
			H[n1][n2] = h[n];
		}
	}
	Y = convolve_dispatch(X, H);
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			const int n = (n1 * N2 + n2 * N1) % N;
			y[n] = Y[n1][n2];
		}
	}
	return y;
}

template<class T, class V>
std::vector<T> convolve_tiny(std::vector<T> x, std::vector<V> h) {
	int N1 = x.size();
	if (N1 == 2) {
		return convolve_fast_2(x, h);
	}
	if (N1 == 3) {
		return convolve_fast_3(x, h);
	}
	if (N1 == 4) {
		return convolve_fast_4(x, h);
	}
	if (N1 == 5) {
		return convolve_fast_5(x, h);
	}
	if (N1 == 7) {
		return convolve_fast_7(x, h);
	}
	if (N1 == 8) {
		return convolve_fast_8(x, h);
	}
	if (N1 == 9) {
		return convolve_fast_9(x, h);
	}
	assert(false);
	return std::vector<T>();
}

int operation_count(std::vector<cmplx> x) {
	return math_vertex::operation_count(x).total();
}

int operation_count(std::vector<std::vector<cmplx>> x) {
	int cnt = 0;
	for (auto y : x) {
		cnt += math_vertex::operation_count(y).total();
	}
	return cnt;
}

std::vector<cmplx> convolve_fast(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	if (N == 1) {
		return std::vector<cmplx>(1, x[0] * h[0]);
	}
	auto pfacs = prime_factorization(N);
	std::vector<cmplx> y(N);
	best_x X;
	X.N = N;
	X.sig = fft_input_signature(x);
	X.opts = 0;
	if (can_fast_cyclic(N)) {
		y = convolve_tiny(x, h);
	} else if (can_agarwal(N)) {
		static std::map<best_x, int> cache;
		if (cache.find(X) == cache.end()) {
			auto p = prime_factorization(N);
			int best_N1;
			int best_cnt = 1000000000;
			for (auto fac : p) {
				int N1 = std::pow(fac.first, fac.second);
				if (!can_fast_cyclic(N) && !(N % 2 == 0 || N % 3 == 0 || N % 5 == 0)) {
					continue;
				}
				int N2 = N / N1;
				int cnt = math_vertex::operation_count(convolve_agarwal(N1, N2, x, h)).total();
				if (cnt < best_cnt) {
					best_cnt = cnt;
					best_N1 = N1;
				}
			}
			cache[X] = best_N1;
		}
		int N1 = cache[X];
		int N2 = N / N1;
		y = convolve_agarwal(N1, N2, x, h);
	} else {
		if (N % 5 == 0) {
			y = convolve_toom(5, x, h);
		} else if (N % 3 == 0) {
			y = convolve_toom(3, x, h);
		} else if (N % 2 == 0) {
			y = convolve_toom(2, x, h);
		} else {
			y = convolve_fft(x, h);
		}
	}
	for (auto z : y) {
		z.set_goal();
	}
	return y;
}

std::vector<cmplx> convolve_slow(std::vector<cmplx> x, std::vector<std::complex<double>> h) {
	int N = x.size();
	std::vector<cmplx> y(N, cmplx( { 0.0, 0.0 }));
	for (int n = 0; n < N; n++) {
		for (int m = 0; m < N; m++) {
			y[n] += x[m] * h[(n - m + N) % N];
		}
	}
	return y;
}

std::vector<cmplx> convolve(std::vector<cmplx> x, std::vector<std::complex<double>> h, int opts) {
	int N = x.size();
	if (N == 1) {
		return std::vector<cmplx>(1, x[0] * h[0]);
	}
	static std::map<best_x, int> cache;
	best_x X;
	X.opts = opts;
	X.N = N;
	X.sig = fft_input_signature(x);
	if (cache.find(X) == cache.end()) {
		int slow_cnt;
		int fast_cnt;
		int fft_cnt;
		std::vector<cmplx> y;
		fast_cnt = math_vertex::operation_count(x * h).total();
		fft_cnt = math_vertex::operation_count(convolve_fft(x, h)).total();
		if (fast_cnt < fft_cnt && fast_cnt > 0) {
			cache[X] = FAST_CONVOLVE;
		} else {
			cache[X] = FFT_CONVOLVE;
		}
	}
	switch (cache[X]) {
	case FFT_CONVOLVE:
		return convolve_fft(x, h);
	case FAST_CONVOLVE:
		return x * h;
	default:
		assert(false);
		return std::vector<cmplx>();
	}
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

