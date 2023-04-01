/*#include "fft.hpp"
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
		muls += (xq.degree() + 1) * (hq.degree() + 1);
		adds += (xq.degree() + 1) * (hq.degree() + 1) - hq.degree() - xq.degree() - 1;
		auto hxq = hx % q[i];
		adds += (q[i].degree()) * (q[i].degree());
		auto shx = s * hxq;
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
