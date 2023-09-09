#include "fft.hpp"

#include <cmath>
#include <complex>
#include <fftw3.h>

#define RADIX2 0
#define RADIX4 1
#define SINGLETON 3
#define RADERS_FFT 4
#define RADERS_FAST 10
#define RADERS_PADDED 5
#define SUZUKI 6
#define PRIME_FACTOR 7
#define TANGENT 8
#define RADER_BRENNER 11
#define COOLEYTUKEY 9
#define NFFT 12

std::vector<cmplx> fft_suzuki(int R, std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_radix4(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_singleton(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_prime_factor(int N1, int N2, std::vector<cmplx> xin, int opts);
std::vector<cmplx> fft_radix2(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_split_prime(int N1, int N2, std::vector<cmplx> x, int opts);
std::vector<cmplx> fft_raders_fast2(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_rader_brenner(std::vector<cmplx> xin, int N, int opts);

static int next_group = 1;

std::vector<math_vertex> fft(std::vector<math_vertex> xin, int N, int opts) {
	int group = next_group++;
	for (auto v : xin) {
		v.set_group_id(group);
	}
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
		} else {
			x[0].x = xin[0];
			x[0].y = xin[1];
			for (int n = 1; n < N; n++) {
				x[n].x = xin[2 * (N - n)];
				x[n].y = xin[2 * (N - n) + 1];
			}
		}
	} else {
		if (opts & FFT_SKEW) {
			for (int n = 0; n < N / 2; n++) {
				x[n].x = math_vertex(0.5) * xin[n];
				x[n + N / 2].x = math_vertex(-0.5) * xin[n];
			}
		} else if (opts & FFT_REAL) {
			for (int n = 0; n < N; n++) {
				x[n].x = xin[n];
			}
		} else if (opts & FFT_DCT1) {
			x[0].x = xin[0];
			x[N / 2].x = xin[(N + 1) / 2];
			for (int n = 1; n < N - n; n++) {
				x[n].x = xin[n];
				x[N - n].x = xin[n];
			}
		} else if (opts & FFT_DCT2) {
			for (int n = 1; n < N - n; n += 2) {
				x[n].x = xin[(n - 1) / 2];
				x[N - n].x = xin[(n - 1) / 2];
			}
		} else if (opts & FFT_DCT3) {
			x[0].x = xin[0];
			for (int n = 1; n < N / 4; n++) {
				x[N - n].x = xin[n];
				x[N - (n + N / 2)].x = -xin[n];
			}
			for (int n = 1; n < N - n; n++) {
				x[n].x = x[N - n].x;
			}
		} else if (opts & FFT_DCT4) {
			for (int n = 0; n < N / 8; n++) {
				x[2 * n + 1].x = xin[n];
			}
			for (int n = N / 4; n < N / 2; n++) {
				x[n + N / 4].x = -x[n].x;
			}
			for (int n = 1; n < N - n; n++) {
				x[N - n].x = x[n].x;
			}
		} else if (opts & FFT_ODD) {
			for (int n = 1; n < N; n += 2) {
				x[n].x = xin[n - 1];
				x[n].y = xin[n];
			}
		} else {
			for (int n = 0; n < N; n++) {
				x[n].x = xin[2 * n];
				x[n].y = xin[2 * n + 1];
			}
		}
	}
	x = fft(x, N, opts, true);

	if (opts & FFT_INV) {
		if (opts & FFT_REAL) {
			xout.resize(N);
			for (int n = 0; n < N; n++) {
				xout[n] = x[n].x;
			}
		} else {
			xout.resize(2 * N);
			for (int n = 0; n < N; n++) {
				xout[2 * n] = x[n].x;
				xout[2 * n + 1] = x[n].y;
			}
		}
	} else {
		if (opts & FFT_SKEW) {
			xout.resize(N / 2);
			for (int n = 0; n < N / 2 - n - 1; n++) {
				xout[n] = x[2 * n + 1].x;
				xout[N / 2 - n - 1] = x[2 * n + 1].y;
			}
			if ((N / 2) % 2 == 1) {
				xout[N / 4] = x[N / 2].x;
			}
		} else if (opts & FFT_REAL) {
			xout.resize(N);
			xout[0] = x[0].x;
			if (N % 2 == 0) {
				xout[N / 2] = x[N / 2].x;
			}
			for (int n = 1; n < N - n; n++) {
				xout[n] = x[n].x;
				xout[N - n] = x[n].y;
			}
		} else if (opts & FFT_DCT1) {
			xout.resize((N + 2) / 2);
			for (int n = 0; n < xout.size(); n++) {
				xout[n] = x[n].x * math_vertex(1);
			}
		} else if (opts & FFT_DCT2) {
			xout.resize(N / 4);
			for (int n = 0; n < xout.size(); n++) {
				xout[n] = x[n].x;
			}
		} else if (opts & FFT_DCT3) {
			xout.resize(N / 4);
			for (int n = 0; n < xout.size(); n++) {
				xout[n] = x[2 * n + 1].x;
			}
		} else if (opts & FFT_DCT4) {
			xout.resize(N / 8);
			for (int n = 0; n < xout.size(); n++) {
				xout[n] = x[2 * n + 1].x;
			}
		} else if (opts & FFT_ODD) {
			xout.resize(N);
			for (int n = 0; n < N / 2; n++) {
				xout[2 * n] = x[n].x;
				xout[2 * n + 1] = x[n].y;
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

int op_count(const std::vector<cmplx>& x, int opts) {
	int N = x.size();
	std::vector<math_vertex> y;
	if (opts & FFT_INV) {
		if (opts & FFT_REAL) {
			for (int n = 0; n < N; n++) {
				y.push_back(x[n].x);
			}
		} else {
			for (int n = 0; n < N; n++) {
				y.push_back(x[n].x);
				y.push_back(x[n].y);
			}
		}
	} else {
		if (opts & FFT_REAL) {
			for (int n = 0; n < (N + 1) / 2; n++) {
				y.push_back(x[n].x);
				y.push_back(x[n].y);
			}
		} else if (opts & FFT_DCT1) {
			for (int n = 0; n < N / 2; n++) {
				y.push_back(x[n].x);
			}
		} else if ((opts & FFT_DCT2) || (opts & FFT_DCT3)) {
			for (int n = 0; n < N / 4; n++) {
				y.push_back(x[n].x);
			}
		} else if (opts & FFT_DCT4) {
			for (int n = 0; n < N / 8; n++) {
				y.push_back(x[n].x);
			}
		} else {
			for (int n = 0; n < N; n++) {
				y.push_back(x[n].x);
				y.push_back(x[n].y);
			}
		}
	}
	return math_vertex::operation_count(y).total();
}

std::set<int> radices_prime_factor(std::map<int, int> pfac, int N, int R = 1) {
	std::set<int> radices;
	if (R != 1 && R != N) {
		radices.insert(R);
	}
	for (auto i = pfac.begin(); i != pfac.end(); i++) {
		auto pfac1 = pfac;
		pfac1.erase(i->first);
		auto tmp = radices_prime_factor(pfac1, N, R * std::pow(i->first, i->second));
		radices.insert(tmp.begin(), tmp.end());
	}
	return radices;
}

std::set<int> radices_cooley_tukey(std::map<int, int> pfac, int N, int R = 1) {
	std::set<int> radices;
	if (R != 1 && R != N) {
		radices.insert(R);
	}
	for (auto i = pfac.begin(); i != pfac.end(); i++) {
		auto pfac1 = pfac;
		pfac1[i->first]--;
		if (pfac1[i->first] <= 0) {
			pfac1.erase(i->first);
		}
		auto tmp = radices_cooley_tukey(pfac1, N, R * i->first);
		radices.insert(tmp.begin(), tmp.end());
	}
	return radices;
}

struct best_y {
	int R;
	int method;
	int cnt;
};

void convolve_reset_cache();

std::map<best_x2, int> best_methods;

void fft_reset_cache() {
	best_methods = decltype(best_methods)();
	convolve_reset_cache();
}

std::string method_name(int i) {
	std::string method;
	switch (i) {
	case RADIX2:
		method += "radix-2";
		break;
	case TANGENT:
		method += "tangent";
		break;
	case RADER_BRENNER:
		method += "rader-brenner";
		break;
	case RADIX4:
		method += "split4";
		break;
	case SINGLETON:
		method += "Singleton";
		break;
	case RADERS_FFT:
		method += "Raders (fft)";
		break;
	case RADERS_FAST:
		method += "Raders (fast)";
		break;
	case RADERS_PADDED:
		method += "Raders (padded)";
		break;
	case SUZUKI:
		method += "Suzuki";
		break;
	case COOLEYTUKEY:
		method += "Cooley-Tukey";
		break;
	case PRIME_FACTOR:
		method += "Good-Thomas";
		break;
	default:
		assert(false);
	}
	return method;
}

std::string get_best_method(int N, int opts) {
	best_x2 x;
	x.opts = opts;
	x.N = N;
	int ii = best_methods[x];
	return method_name(ii);
}

static std::map<best_x, best_y> best;

void print_fft_bests() {
	for (auto i = best.begin(); i != best.end(); i++) {
		std::string opts;
		if (i->first.opts & FFT_INV) {
			opts += "inverse ";
		}
		if (i->first.opts & FFT_REAL) {
			opts += "real";
		} else if (i->first.opts & FFT_DCT1) {
			opts += "DCT-I";
		} else if (i->first.opts & FFT_DCT2) {
			opts += "DCT-II";
		} else if (i->first.opts & FFT_DCT3) {
			opts += "DCT-III";
		} else if (i->first.opts & FFT_DCT4) {
			opts += "DCT-IV";
		} else {
			opts += "complex";
		}
		opts += " N = " + std::to_string(i->first.N);
		if (i->first.opts & 0xFFFF0000) {
			int M = (i->first.opts >> 16);
			opts += " M = " + std::to_string(M);
		}
		return;
		std::string method = method_name(i->second.method);
//		if (i->second.method == SPLIT_PRIME) {
		fprintf(stderr, "%32s | %16s | %i x %i \n", opts.c_str(), method.c_str(), i->second.R, i->first.N / i->second.R);
		for (int n = 0; n < i->first.sig.size(); n++) {
			fprintf(stderr, "%i ", i->first.sig[n]);
		}
		fprintf( stderr, "\n");
//		}
	}
}

std::vector<cmplx> fft(std::vector<cmplx> xin, int N, int opts, bool root) {
	for (auto x : xin) {
		x.x.set_goal();
		x.y.set_goal();
	}
	if (N == 1) {
		return xin;
	}
	int group = next_group++;
	for (auto v : xin) {
		v.x.set_group_id(group);
		v.y.set_group_id(group);
	}
	std::vector<cmplx> xout;
	auto pfac = prime_factorization(N);
	std::vector<int> factors;
	for (auto i = pfac.begin(); i != pfac.end(); i++) {
		factors.push_back(i->first);
		if (i->second > 1) {
			factors.push_back(std::pow(i->first, i->second));
		}
		if (i->first == 2 && i->second >= 2) {
			factors.push_back(4);
		}
	}
	best_x X;
	X.N = N;
	X.sig = fft_input_signature(xin);
	bool allzero = true;
	for (int n = 0; n < X.sig.size(); n++) {
		if (X.sig[n]) {
			allzero = false;
			break;
		}
	}
	if (allzero) {
		return xin;
	}
	X.opts = opts;
	if (best.find(X) == best.end()) {
		int huge = std::numeric_limits<int>::max();
		std::vector<best_y> tries;
		best_y y;
		auto rpf = radices_prime_factor(pfac, N);
		for (auto R : rpf) {
			y.R = R;
			y.method = PRIME_FACTOR;
			y.cnt = op_count(fft_prime_factor(R, N / R, xin, opts), opts);
			tries.push_back(y);
		}
		auto rct = radices_cooley_tukey(pfac, N);
		for (auto R : rct) {
			y.R = R;
			y.method = COOLEYTUKEY;
			y.cnt = op_count(fft_cooley_tukey(R, N / R, xin, opts), opts);
			tries.push_back(y);
		}
		if (N % 2 == 0) {
			y.R = 2;
			y.method = RADIX2;
			y.cnt = op_count(fft_radix2(xin, N, opts), opts);
			tries.push_back(y);
			if (N > 2 && N % 4 == 0) {
				y.R = 4;
				y.method = RADIX4;
				y.cnt = op_count(fft_radix4(xin, N, opts), opts);
				tries.push_back(y);
				y.R = 4;
				y.method = RADER_BRENNER;
				y.cnt = op_count(fft_rader_brenner(xin, N, opts), opts);
				tries.push_back(y);
			}
			if (1 << ilogb(N) == N) {
				int cnt1, cnt2;
				y.R = 4;
				y.method = TANGENT;
				cnt1 = y.cnt = op_count(fft_modsplit(xin, N, opts), opts);
				tries.push_back(y);
			}
		} else {
			int rcnt, pcnt;
			if (is_prime(N)) {
				y.R = N;
				y.method = SINGLETON;
				y.cnt = op_count(fft_singleton(xin, N, opts), opts);
				tries.push_back(y);
			}
			if (pfac.size() == 1) {
				if (pfac.begin()->second > 1 && pfac.begin()->first == 3) {
					y.R = pfac.begin()->first;
					y.method = SUZUKI;
					y.cnt = op_count(fft_suzuki(y.R, xin, N, opts), opts);
					tries.push_back(y);
				}
				y.R = N;
				y.method = RADERS_FFT;
				rcnt = y.cnt = op_count(fft_raders_fft(xin, N, false, opts), opts);
				tries.push_back(y);
				if (N <= 48) {
					y.R = N;
					y.method = RADERS_FAST;
					rcnt = y.cnt = op_count(fft_raders_fast(xin, N, opts), opts);
					tries.push_back(y);
				}
				if (!(opts & !FFT_INV)) {
					y.R = N;
					y.method = RADERS_PADDED;
					pcnt = y.cnt = op_count(fft_raders_fft(xin, N, true, opts), opts);
					tries.push_back(y);
				}
			}
		}
		int besti = -1;
		int bestcnt = huge;
		for (int i = 0; i < tries.size(); i++) {
			if (tries[i].cnt < bestcnt) {
				bestcnt = tries[i].cnt;
				besti = i;
			}
		}
		assert(besti != -1);
		best[X] = tries[besti];
	}
	auto Y = best[X];
	if (root) {
		best_x2 x2;
		x2.N = N;
		x2.opts = opts;
		best_methods[x2] = Y.method;
	}
	auto R = Y.R;
	switch (Y.method) {
	case RADIX2:
		xout = fft_radix2(xin, N, opts);
		break;
	case RADIX4:
		xout = fft_radix4(xin, N, opts);
		break;
	case SINGLETON:
		xout = fft_singleton(xin, N, opts);
		break;
	case RADERS_FAST:
		xout = fft_raders_fast(xin, N, opts);
		break;
	case RADERS_FFT:
		xout = fft_raders_fft(xin, N, false, opts);
		break;
	case TANGENT:
		xout = fft_modsplit(xin, N, opts);
		break;
	case RADER_BRENNER:
		xout = fft_rader_brenner(xin, N, opts);
		break;
	case RADERS_PADDED:
		xout = fft_raders_fft(xin, N, true, opts);
		break;
	case SUZUKI:
		xout = fft_suzuki(R, xin, N, opts);
		break;
	case PRIME_FACTOR:
		xout = fft_prime_factor(R, N / R, xin, opts);
		break;
	case COOLEYTUKEY:
		xout = fft_cooley_tukey(R, N / R, xin, opts);
		break;
	default:
		assert(false);
	}
//	best.erase(X);
	group = next_group++;
	for (auto v : xout) {
		v.x.set_group_id(group);
		v.y.set_group_id(group);
	}
	for (auto x : xout) {
		x.x.set_goal();
		x.y.set_goal();
	}
	return std::move(xout);
}

