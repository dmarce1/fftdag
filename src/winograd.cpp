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

std::vector<std::complex<double>> winograd_convolve(const std::vector<std::complex<double>>& X, const std::vector<std::complex<double>>& H);

std::vector<poly> polynomial_transform(const std::vector<poly>& xn, poly m, int P) {
	std::vector<poly> Xk(P);
	for (int k = 0; k < P; k++) {
		for (int r = 0; r < P; r++) {
			poly root;
			root[r * k] = 1.0;
			Xk[k] = (Xk[k] + xn[r] * root) % m;
		}
	}
	return Xk;
}

std::vector<poly> inverse_polynomial_transform(const std::vector<poly>& xk, poly m, int P) {
	std::vector<poly> Xn(P);
	for (int r = 0; r < P; r++) {
		for (int k = 0; k < P; k++) {
			poly root;
			root[(P - 1) * r * k] = 1.0 / P;
			Xn[r] = (Xn[r] + xk[k] * root) % m;
		}
	}
	return Xn;
}

std::vector<std::vector<std::complex<double>>>convolve2d(const std::vector<std::vector<std::complex<double>>>& X, const std::vector<std::vector<std::complex<double>>>& H) {
	int N1 = X.size();
	int N2 = X[0].size();
	poly m, p1, p2;
	std::vector<poly> x(N1);
	std::vector<poly> h(N1);
	std::vector<poly> h2(N1);
	std::vector<poly> x2(N1);
	std::vector<poly> h1(N1);
	std::vector<poly> x1(N1);
	std::vector<poly> y2(N1);
	std::vector<poly> y(N1);
	std::vector<std::vector<std::complex<double>>> yr(N1, std::vector<std::complex<double>>(N2));
	m[N2] = 1.0;
	m[0] = -1.0;
	p2[1] = 1.0;
	p2[0] = -1.0;
	p1 = m / p2;
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			x[n1][n2] = X[n1][n2];
			h[n1][n2] = H[n1][n2];
		}
	}
	for( int n1 = 0; n1 < N1; n1++) {
		x1[n1] = x[n1] % p1;
		h1[n1] = h[n1] % p1;
		x2[n1] = x[n1] % p2;
		h2[n1] = h[n1] % p2;
	}
	auto X1k = polynomial_transform(x1, p1, N2);
	auto H1k = polynomial_transform(h1, p1, N2);
	std::vector<poly> Y1k(N1);
	for( int n = 0; n < N1; n++) {
		Y1k[n] = (X1k[n] * H1k[n]) % p1;
	}
	auto y1 = inverse_polynomial_transform(Y1k, p1, N2);
	for( int n = 0; n < N1; n++) {
		for( int m = 0; m < N1; m++) {
			y2[n] += (x2[m] * h2[mod(n - m, N1)]) % p2;
		}
	}
	for( int n1 = 0; n1 < N1; n1++) {
		poly t1 = (m / p1) % m;
		poly t2 = (m / p2) % m;
		poly s1 = (t1 * inverse(t1 % p1, p1)) % m;
		poly s2 = (t2 * inverse(t2 % p2, p2)) % m;
		y[n1] = (y1[n1] * s1 + y2[n1] * s2) % m;
	}
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			yr[n1][n2] = y[n1][n2];
		}
	}
	return yr;
}

void test_poly() {
	constexpr int N1 = 5;
	constexpr int N2 = 5;
	std::vector<std::vector<std::complex<double>>>X(N1,std::vector<std::complex<double>>(N2));
	auto H = X;
	auto Y0 = X;
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			X[n1][n2].real(rand1());
			H[n1][n2].real(rand1());
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
	auto Y1 = convolve2d(X, H);
	double avg_err = 0.0;
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			double err = std::abs(Y0[n1][n2] - Y1[n1][n2]);
			avg_err += err / (N1 * N2);
			printf("%i %i | %e %e | %e %e | %e\n", n1, n2, Y0[n1][n2].real(), Y0[n1][n2].imag(), Y1[n1][n2].real(), Y1[n1][n2].imag(), err);
		}
	}
	printf("%e\n", avg_err);
}

std::vector<std::complex<double>> winograd_convolve(const std::vector<std::complex<double>>& X, const std::vector<std::complex<double>>& H) {
	int N = X.size();
	poly m, y, x, h;
	std::vector<poly> q;
	std::vector<std::complex<double>> yr(N);
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
	int muls = 0;
	int adds = 0;
	for (int i = 0; i < K; i++) {
		auto t = (m / q[i]) % m;
		auto s = (t * inverse(t % q[i], q[i])) % m;
		auto hq = h % q[i];
		auto xq = x % q[i];
		adds += (q[i].degree()) * (q[i].degree());
		auto hx = xq * hq;
		auto hxq = hx % q[i];
		auto shx = (s * hxq) % m;
		muls += (xq.degree() + 1) * (hq.degree() + 1);
		adds += (xq.degree() + 1) * (hq.degree() + 1) - hq.degree() - xq.degree() - 1;
		adds += (q[i].degree()) * (q[i].degree());
		adds += (hxq.degree() + 1) * (s.degree() + 1) - hxq.degree() - s.degree() - 1;
		y += shx;
		printf("x =      %s\n", x.to_str().c_str());
		printf("h =      %s\n", h.to_str().c_str());
		printf("q =      %s\n", q[i].to_str().c_str());
		printf("s =      %s\n", s.to_str().c_str());
		printf("xq =      %s\n", xq.to_str().c_str());
		printf("hq =     %s\n", hq.to_str().c_str());
		printf("hx =     %s\n", hx.to_str().c_str());
		printf("hxq =    %s\n", hxq.to_str().c_str());
		printf("shx =    %s\n", shx.to_str().c_str());
		printf("\n\n");
	}
	for (int n = 0; n < N; n++) {
		yr[n] = y[n];
	}
	printf("%i adds %i muls\n", adds, muls);
	return yr;
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
