#include "fft.hpp"
#include "polynomial.hpp"
#include "util.hpp"

using poly = polynomial<std::complex<double>>;

poly cyclotomic(int n) {
	poly P;
	for (int k = 1; k <= n; k++) {
		if (are_coprime(n, k)) {
			poly R;
			R[1] = 1.0;
			R[0] = -std::polar(1.0, 2.0 * M_PI * k / n);
			if (k == 1) {
				P = R;
			} else {
				P *= R;
			}
		}
	}
	return P;
}
template<class T>
polynomial<T> inverse(polynomial<T> f, polynomial<T> m) {
	if (f.degree() == 0) {
		f[0] = 1.0 / f[0];
		return f;
	}
	polynomial<T> t2, t1, t0;
	polynomial<T> r2, r1, r0;
	t1[0] = 1.0;
	r0 = m;
	r1 = f;
	do {
		polynomial<T> quotient = (r0 / r1);
		//	printf("%32s | %32s | %32s  | %32s  | %32s \n", r0.to_str().c_str(), r1.to_str().c_str(), (quotient).to_str().c_str(),  r2.to_str().c_str(),  t2.to_str().c_str(), ((f*t2)%m).to_str().c_str());
		t2 = (t0 - quotient * t1);
		r2 = (r0 - quotient * r1);
		r2 = (r2 + m) % m;
		r0 = r1;
		r1 = r2;
		t0 = t1;
		t1 = t2 % m;
		if (((t1 * f) % m).degree() == 0) {
			break;
		}
	} while (!r1.zero() && !close2(std::abs(r1[0]), 0.0));
	auto c = 1.0 / ((t1 * f) % m)[0];
	for (int n = 0; n <= t1.degree(); n++) {
		t1[n] *= c;
	}
	return t1;
}

static double rand1() {
	return 2.0 * (rand() + 0.5) / (RAND_MAX + 1.0) - 1.0;
}

/*muls += (xq.degree() + 1) * (hq.degree() + 1);
 adds += (xq.degree() + 1) * (hq.degree() + 1) - hq.degree() - xq.degree() - 1;
 adds += (q[i].degree()) * (q[i].degree());
 adds += (hxq.degree() + 1) * (s.degree() + 1) - hxq.degree() - s.degree() - 1;
 */

std::vector<std::complex<double>> winograd_convolve(const std::vector<std::complex<double>>& X, const std::vector<std::complex<double>>& H);

template<class T, class U>
std::vector<polynomial<T>> polynomial_transform(const std::vector<polynomial<T>>& xn, polynomial<U> m, int P, int r = 1) {
	using polyT = polynomial<T>;
	using polyU = polynomial<U>;
	int N = xn.size();
	std::vector<polyT> xk(N);
	if (N <= 2 || N % 2 != 0) {
		using polyT = polynomial<T>;
		using polyU = polynomial<U>;
		for (int k = 0; k < N; k++) {
			for (int n = 0; n < N; n++) {
				polyU w;
				w[r * n * k] = U(1);
				xk[k] += (xn[n] * w) % m;
			}
		}
	} else {
		std::vector<polyT> xe(N / 2);
		std::vector<polyT> xo(N / 2);
		for (int n = 0; n < N / 2; n++) {
			xe[n] = xn[2 * n];
			xo[n] = xn[2 * n + 1];
		}
		xe = polynomial_transform(xe, m, N / 2, 2 * r);
		xo = polynomial_transform(xo, m, N / 2, 2 * r);
		for (int k = 0; k < N / 2; k++) {
			poly w1;
			poly w2;
			w1[k * r] = T(1);
			w1 = w1 % m;
			w2[(k + N / 2) * r] = T(1);
			w2 = w2 % m;
			xk[k] = (xe[k] + w1 * xo[k]) % m;
			xk[k + N / 2] = (xe[k] + w2 * xo[k]) % m;
		}
	}
	return xk;
}

template<class T, class U>
std::vector<polynomial<T>> inverse_polynomial_transform(const std::vector<polynomial<T>>& xk, polynomial<U> m, int P, int n = 1) {
	using polyT = polynomial<T>;
	using polyU = polynomial<U>;
	std::vector<polyT> Xn(P);
	for (int r = 0; r < P; r++) {
		for (int k = 0; k < P; k++) {
			polyU root;
			root[(P - 1) * n * r * k] = 1.0 / P;
			Xn[r] = (Xn[r] + xk[k] * root) % m;
		}
	}
	return Xn;
}

template<class T, class U>
polynomial<T> convolve_1d(const polynomial<T>& x, const polynomial<U>& h, int N) {
	poly m, y;
	std::vector<poly> q;
	m[0] = U(-1);
	m[N] = U(1);
	for (int d = 1; d <= N; d++) {
		if (N % d == 0) {
			q.push_back(cyclotomic(d));
		}
	}
	int K = q.size();
	for (int i = 0; i < K; i++) {
		auto t = (m / q[i]) % m;
		auto s = (t * inverse(t % q[i], q[i])) % m;
		auto hq = h % q[i];
		auto xq = x % q[i];
		auto hx = xq * hq;
		auto hxq = hx % q[i];
		auto shx = (s * hxq) % m;
		y += shx;
	}
	return y;
}

template<class T, class U>
std::vector<polynomial<T>> convolve_2d_2t(const std::vector<polynomial<T>>& x0, const std::vector<polynomial<U>>& h0) {
	using polyT = polynomial<T>;
	using polyU = polynomial<U>;
	int N = x0.size();
	int t = ilogb(N);
	if (t == 0) {
		std::vector<polynomial<T>> rc;
		polynomial<T> product;
		product = x0[0] * h0[0];
		rc.push_back(product);
		return rc;
	}
	std::vector<polyT> x1(N);
	std::vector<polyU> h1(N);
	std::vector<polyU> h2(N);
	std::vector<polyT> x2(N);
	std::vector<polyU> h3(N / 2);
	std::vector<polyT> x3(N / 2);
	std::vector<polyU> h4(N / 2);
	std::vector<polyT> x4(N / 2);
	std::vector<polyU> h5(N / 2);
	std::vector<polyT> x5(N / 2);
	std::vector<polyT> y4k(N / 2);
	std::vector<polyT> y2(N / 2);
	std::vector<polyT> y0(N);
	std::vector<polyT> y5k(N / 2);
	std::vector<polyT> y1k(N);
	polyU m, p1, p2;
	polyT tx;
	polyU th;
	m[1 << t] = U(1);
	m[0] = U(-1);
	p1[1 << (t - 1)] = U(1);
	p1[0] = U(+1);
	p2[1 << (t - 1)] = U(1);
	p2[0] = U(-1);
	for (int n = 0; n < N; n++) {
		x1[n] = x0[n] % p1;
		h1[n] = h0[n] % p1;
		x2[n] = x0[n] % p2;
		h2[n] = h0[n] % p2;
	}
	auto x1k = polynomial_transform(x1, m, N);
	auto h1k = polynomial_transform(h1, m, N);
	for (int k = 0; k < N; k++) {
		y1k[k] = x1k[k] * h1k[k];
	}
	auto y1 = inverse_polynomial_transform(y1k, m, N);
	for (int n = 0; n < N / 2; n++) {
		for (int m = 0; m < N; m++) {
			h3[n][m] = h2[m][n];
			x3[n][m] = x2[m][n];
		}
	}
	for (int n = 0; n < N / 2; n++) {
		x4[n] = x3[n] % p1;
		x5[n] = x3[n] % p2;
		h4[n] = h3[n] % p1;
		h5[n] = h3[n] % p2;
	}
	auto x4k = polynomial_transform(x4, m, N / 2, 2);
	auto h4k = polynomial_transform(h4, m, N / 2, 2);
	for (int k = 0; k < N / 2; k++) {
		y4k[k] = x4k[k] * h4k[k];
	}
	auto y4 = inverse_polynomial_transform(y4k, m, N / 2, 2);
	auto y5 = convolve_2d_2t(x5, h5);
	for (int n = 0; n < N / 2; n++) {
		polyU t1 = (m / p1) % m;
		polyU t2 = (m / p2) % m;
		polyU s1 = (t1 * inverse(t1 % p1, p1)) % m;
		polyU s2 = (t2 * inverse(t2 % p2, p2)) % m;
		y2[n] = (y4[n] * s1 + y5[n] * s2) % m;
	}
	for (int n = 0; n < N / 2; n++) {
		for (int l = 0; l < N; l++) {
			y0[l][n] = y2[n][l];
		}
	}
	for (int n = 0; n < N; n++) {
		polyU t1 = (m / p1) % m;
		polyU t2 = (m / p2) % m;
		polyU s1 = (t1 * inverse(t1 % p1, p1)) % m;
		polyU s2 = (t2 * inverse(t2 % p2, p2)) % m;
		y0[n] = (y1[n] * s1 + y0[n] * s2) % m;
	}
	return y0;

}

template<class T, class U>
std::vector<polynomial<T>> convolve_2d_p(const std::vector<polynomial<T>>& x0, const std::vector<polynomial<U>>& h0) {
	using polyT = polynomial<T>;
	using polyU = polynomial<U>;
	int P = x0.size();
	poly m, p1, p2;
	std::vector<polyT> x1(P);
	std::vector<polyU> h1(P);
	std::vector<polyU> h2(P);
	std::vector<polyT> x2(P);
	std::vector<polyT> y0(P);
	std::vector<polyT> y2(P);
	polyT tx;
	polyU th;
	m[P] = U(1);
	m[0] = U(-1);
	p2[1] = U(1);
	p2[0] = U(-1);
	p1 = m / p2;
	for (int n1 = 0; n1 < P; n1++) {
		x1[n1] = x0[n1] % p1;
		h1[n1] = h0[n1] % p1;
		x2[n1] = x0[n1] % p2;
		h2[n1] = h0[n1] % p2;
	}
	auto X1k = polynomial_transform(x1, p1, P);
	auto H1k = polynomial_transform(h1, p1, P);
	std::vector<polyT> Y1k(P);
	for (int n = 0; n < P; n++) {
		Y1k[n] = (X1k[n] * H1k[n]) % p1;
	}
	auto y1 = inverse_polynomial_transform(Y1k, p1, P);
	for (int n = 0; n < P; n++) {
		tx[n] = x2[n][0];
		th[n] = h2[n][0];
	}
	auto ty = convolve_1d(tx, th, P);
	for (int n = 0; n < P; n++) {
		y2[n][0] = ty[n];
	}
	for (int n1 = 0; n1 < P; n1++) {
		poly t1 = (m / p1) % m;
		poly t2 = (m / p2) % m;
		poly s1 = (t1 * inverse(t1 % p1, p1)) % m;
		poly s2 = (t2 * inverse(t2 % p2, p2)) % m;
		y0[n1] = (y1[n1] * s1 + y2[n1] * s2) % m;
	}
	return y0;
}

template<class T, class U>
std::vector<polynomial<T>> convolve_2d_p2(const std::vector<polynomial<T>>& x0, const std::vector<polynomial<U>>& h0) {
	using polyT = polynomial<T>;
	using polyU = polynomial<U>;
	int P2 = x0.size();
	std::vector<std::vector<T>> yr(P2, std::vector<T>(P2));
	int P = lround(sqrt(P2));
	poly m, r1, r2;
	std::vector<polyT> h1(P2);
	std::vector<polyT> h2(P2);
	std::vector<polyT> h3(P);
	std::vector<polyT> h4(P);
	std::vector<polyT> h5(P);
	std::vector<polyU> x1(P2);
	std::vector<polyU> x2(P2);
	std::vector<polyU> x3(P);
	std::vector<polyU> x4(P);
	std::vector<polyU> x5(P);
	std::vector<polyU> y0(P2);
	std::vector<polyU> y2(P);
	std::vector<polyU> y3(P2);
	std::vector<polyU> y5(P);
	std::vector<polyU> y1k(P2);
	std::vector<polyU> y4k(P2);
	polyT tx;
	polyU th;
	m[P2] = U(1);
	m[0] = U(-1);
	r2[P] = U(1);
	r2[0] = U(-1);
	r1 = (m / r2);
	for (int p = 0; p < P2; p++) {
		x1[p] = x0[p] % r1;
		h1[p] = h0[p] % r1;
		x2[p] = x0[p] % r2;
		h2[p] = h0[p] % r2;
	}
	for (int p = 0; p < P; p++) {
		for (int q = 0; q < P2; q++) {
			h3[p][q] = h2[q][p];
			x3[p][q] = x2[q][p];
		}
	}
	for (int p = 0; p < P; p++) {
		h4[p] = h3[p] % r1;
		x4[p] = x3[p] % r1;
		h5[p] = h3[p] % r2;
		x5[p] = x3[p] % r2;
	}
	auto x1k = polynomial_transform(x1, r1, P2);
	auto h1k = polynomial_transform(h1, r1, P2);
	for (int p = 0; p < P2; p++) {
		y1k[p] = x1k[p] * h1k[p];
	}
	auto y1 = inverse_polynomial_transform(y1k, r1, P2);
	auto x4k = polynomial_transform(x4, r1, P, P);
	auto h4k = polynomial_transform(h4, r1, P, P);
	for (int p = 0; p < P; p++) {
		y4k[p] = x4k[p] * h4k[p];
	}
	auto y4 = inverse_polynomial_transform(y4k, r1, P, P);
	for (int p = 0; p < P; p++) {
		for (int q = 0; q < P; q++) {
			y5[p] += (x5[q] * h5[mod(p - q, P)]) % r2;
		}
	}
	y5 = convolve_2d_p(x5, h5);
	for (int p = 0; p < P; p++) {
		polyU t1 = (m / r1) % m;
		polyU t2 = (m / r2) % m;
		polyU s1 = (t1 * inverse(t1 % r1, r1)) % m;
		polyU s2 = (t2 * inverse(t2 % r2, r2)) % m;
		y2[p] = (y4[p] * s1 + y5[p] * s2) % m;
	}
	for (int p = 0; p < P; p++) {
		for (int q = 0; q < P2; q++) {
			y0[q][p] = y2[p][q];
		}
	}
	for (int p = 0; p < P2; p++) {
		polyU t1 = (m / r1) % m;
		polyU t2 = (m / r2) % m;
		polyU s1 = (t1 * inverse(t1 % r1, r1)) % m;
		polyU s2 = (t2 * inverse(t2 % r2, r2)) % m;
		y0[p] = (y1[p] * s1 + y0[p] * s2) % m;
	}
	return y0;
}

void test_poly() {
	constexpr int N1 = 8;
	constexpr int N2 = 8;
	polynomial<std::complex<double>> tx, th, ty0, ty1;
	for (int n = 0; n < N1; n++) {
		tx[n].real(rand1());
		th[n].real(rand1());
		tx[n].imag(0);
		th[n].imag(0);
	}
	for (int n = 0; n < N1; n++) {
		for (int m = 0; m < N1; m++) {
			ty0[n] += tx[m] * th[mod(n - m, N1)];
		}
	}
	ty1 = convolve_1d(tx, th, N1);
	double avg_err = 0.0;
	for (int n1 = 0; n1 < N1; n1++) {
		double err = std::abs(ty0[n1] - ty1[n1]);
		avg_err += err / (N1 * N2);
		printf("%i | %e %e | %e %e | %e\n", n1, ty0[n1].real(), ty0[n1].imag(), ty1[n1].real(), ty1[n1].imag(), err);
	}
	printf("%e\n", avg_err);

	std::vector<polynomial<std::complex<double>>> X(N1,polynomial<std::complex<double>>());
	auto H = X;
	auto Y0 = X;
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			X[n1][n2].real(n1 * n2 + n1);
			H[n1][n2].real(sqrt(n1 * n2) - n2);
			X[n1][n2].imag(0);
			H[n1][n2].imag(0);
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			Y0[n1][n2] = 0.0;
			for (int m1 = 0; m1 < N1; m1++) {
				for (int m2 = 0; m2 < N2; m2++) {
					Y0[n1][n2] += X[m1][m2] * H[mod(n1 - m1, N1)][mod(n2 - m2, N2)];
				}
			}
		}
	}
	auto Y1 = convolve_2d_2t(X, H);
	avg_err = 0.0;
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			double err = std::abs(Y0[n1][n2] - Y1[n1][n2]);
			avg_err += err / (N1 * N2);
			printf("%i %i | %e %e | %e %e | %e\n", n1, n2, Y0[n1][n2].real(), Y0[n1][n2].imag(), Y1[n1][n2].real(), Y1[n1][n2].imag(), err);
		}
	}
	printf("%e\n", avg_err);
}

/*
 std::vector<cmplx> winograd_convolve(std::vector<cmplx> X, std::vector<std::complex<double>> H) {
 int N = X.size();
 polynomial<cmplx> y, x;
 polynomial<std::complex<double>> m, h;
 std::vector<polynomial<std::complex<double>>> q;
 std::vector<cmplx> yr(N);
 for (int n = 0; n < N; n++) {
 x[n] = X[n];
 h[n] = H[n];
 }
 m[0] = -1.0;
 m[N] = 1.0;
 for (int d = 1; d <= N; d++) {
 if (N % d == 0) {
 q.push_back(cyclotomic(d));
 }
 }
 int K = q.size();
 for (int i = 0; i < K; i++) {
 polynomial<std::complex<double>> t = (m / q[i]) % m;
 polynomial<std::complex<double>> s = (t * inverse(t % q[i], q[i])) % m;
 polynomial<std::complex<double>> hq = h % q[i];
 polynomial<cmplx> xq = x % q[i];
 polynomial<cmplx> hx = xq * hq;
 polynomial<cmplx> hxq = hx % q[i];
 polynomial<cmplx> shx = hxq * s;
 y += shx;
 }
 for (int n = 0; n < N; n++) {
 yr[n] = y[n];
 }
 return yr;
 }
 */
