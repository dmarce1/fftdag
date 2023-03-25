#include "fft.hpp"

#include <cmath>
#include <complex>
#include <fftw3.h>

#define RADIX2 0
#define RADIX4 1
#define SINGLETON 2
#define RADERS 3
#define RADERS_PADDED 4
#define PRIME_POWER 5
#define PRIME_FACTOR 6
#define NFFT 7

std::vector<cmplx> fft_prime_power(int R, std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_radix4(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_singleton(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_prime_factor(int N1, int N2, std::vector<cmplx> xin, int opts);
std::vector<cmplx> fft_raders(std::vector<cmplx> xin, int N, bool padded, int opts);
std::vector<cmplx> fft(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_radix2(std::vector<cmplx> xin, int N, int opts);

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
		if (opts & FFT_REAL) {
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
		} else {
			for (int n = 0; n < N; n++) {
				x[n].x = xin[2 * n];
				x[n].y = xin[2 * n + 1];
			}
		}
	}
	x = fft(x, N, opts);

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

int op_count(const std::vector<cmplx>& xs) {
	std::vector<math_vertex> outs;
	for (auto x : xs) {
		outs.push_back(x.x);
		outs.push_back(x.y);
	}
	return math_vertex::operation_count(outs).total();
}

bool are_coprime(int a, int b) {
	auto afacs = prime_factorization(a);
	auto bfacs = prime_factorization(b);
	for (auto i : afacs) {
		for (auto j : bfacs) {
			if (i.first == j.first) {
				return false;
			}
		}
	}
	return true;
}

std::vector<int> fft_input_signature(std::vector<cmplx> xin) {
	std::vector<int> sig;
	std::vector<math_vertex> X;
	for (auto x : xin) {
		X.push_back(x.x);
		X.push_back(x.y);
	}
	for (int i = 0; i < X.size(); i++) {
		int xid;
		bool neg = X[i].is_neg();
		if (neg) {
			X[i] = X[i].get_neg();
		}
		if (X[i].is_zero()) {
			sig.push_back(0);
		} else {
			xid = X[i].get_unique_id();
			if (neg) {
				sig.push_back(-xid);
			} else {
				sig.push_back(xid);
			}
		}
	}
	int next = 1;
	std::vector<bool> done(sig.size(), false);
	auto tmp = std::move(sig);
	sig.resize(tmp.size());
	for (int i = 0; i < tmp.size(); i++) {
		if (done[i]) {
			continue;
		}
		if (tmp[i] == 0) {
			sig[i] = 0;
			done[i] = true;
			continue;
		}
		int index = next++;
		sig[i] = index;
		for (int j = i + 1; j < tmp.size(); j++) {
			if (done[j]) {
				continue;
			}
			if (tmp[j] == tmp[i]) {
				sig[j] = index;
				done[j] = true;
			} else if (tmp[j] == -tmp[i]) {
				sig[j] = -index;
				done[j] = true;
			}
		}
		done[i] = true;
	}
/*	fprintf(stderr, "%i: ", sig.size());
	for (int n = 0; n < sig.size(); n++) {
		fprintf(stderr, "%i ", sig[n]);
	}
	fprintf(stderr, "\n");*/
	return sig;
}

struct best_x {
	int N;
	int opts;
	std::vector<int> sig;
	bool operator==(const best_x& other) const {
		return N == other.N && opts == other.opts && sig == other.sig;
	}
	bool operator<(const best_x& other) const {
		if (N < other.N) {
			return true;
		} else if (N > other.N) {
			return false;
		} else if (opts < other.opts) {
			return true;
		} else if (opts > other.opts) {
			return false;
		} else {
			for (int n = 0; n < N; n++) {
				if (sig[n] < other.sig[n]) {
					return true;
				} else if (sig[n] > other.sig[n]) {
					return false;
				}
			}
		}
		return false;
	}
};

struct best_y {
	int R;
	int method;
};

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
		std::string method;
		switch (i->second.method) {
		case RADIX2:
			method += "radix-2";
			break;
		case RADIX4:
			if (1 << ilogb(i->first.N) == i->first.N) {
				method += "tangent";
			} else {
				method += "split";
			}
			break;
		case SINGLETON:
			method += "Singleton";
			break;
		case RADERS:
			method += "Raders";
			break;
		case RADERS_PADDED:
			method += "Raders (padded)";
			break;
		case PRIME_POWER:
			method += "radix-n";
			break;
		case PRIME_FACTOR:
			method += "Good-Thomas";
			break;
		default:
			assert(false);
		}
		fprintf(stderr, "%32s | %16s | %i\n", opts.c_str(), method.c_str(), i->second.R);
	}
}

std::vector<cmplx> fft(std::vector<cmplx> xin, int N, int opts) {
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
	X.opts = opts;
	if (best.find(X) == best.end()) {
		int huge = std::numeric_limits<int>::max();
		std::vector<std::array<int, NFFT>> results(factors.size());
		for (int i = 0; i < factors.size(); i++) {
			for (int j = 0; j < NFFT; j++) {
				results[i][j] = huge;
			}
		}
		for (int i = 0; i < factors.size(); i++) {
			int R = factors[i];
			int rc;
			if (R == 2) {
				results[i][RADIX2] = op_count(fft_radix2(xin, N, opts));
			}
			if (R == 4) {
				if (1 << ilogb(N) == N) {
					results[i][RADIX4] = op_count(fft_modsplit(xin, N, opts));
				} else {
					results[i][RADIX4] = op_count(fft_radix4(xin, N, opts));
				}
			}
			if (R != 2 && is_prime(R) && N == R) {
				results[i][SINGLETON] = op_count(fft_singleton(xin, N, opts));
			}
			if (R != 2 && is_prime(R) && N == R) {
				results[i][RADERS] = op_count(fft_raders(xin, N, false, opts));
			}
			if (R != 2 && is_prime(R) && N == R) {
				results[i][RADERS_PADDED] = op_count(fft_raders(xin, N, true, opts));
			}
			if (R % 2 != 0 && is_prime(R) && (pfac.size() == 1) && N != R) {
				results[i][PRIME_POWER] = op_count(fft_prime_power(R, xin, N, opts));
			}
			if (are_coprime(R, N / R) && N != R) {
				results[i][PRIME_FACTOR] = op_count(fft_prime_factor(R, N / R, xin, opts));
			}
		}
		int besti;
		int bestj = -1;
		int bestcnt = huge;
		for (int i = 0; i < factors.size(); i++) {
			for (int j = 0; j < NFFT; j++) {
				if (results[i][j] < bestcnt) {
					bestcnt = results[i][j];
					besti = i;
					bestj = j;
				}
			}
		}
		assert(bestj != -1);
		assert(factors[besti] > 1);
		best_y Y;
		Y.R = factors[besti];
		Y.method = bestj;
		best[X] = Y;
	}
	auto Y = best[X];
	auto R = Y.R;
	switch (Y.method) {
	case RADIX2:
		xout = fft_radix2(xin, N, opts);
		break;
	case RADIX4:
		if (1 << ilogb(N) == N) {
			xout = fft_modsplit(xin, N, opts);
		} else {
			xout = fft_radix4(xin, N, opts);
		}
		break;
	case SINGLETON:
		xout = fft_singleton(xin, N, opts);
		break;
	case RADERS:
		xout = fft_raders(xin, N, false, opts);
		break;
	case RADERS_PADDED:
		xout = fft_raders(xin, N, true, opts);
		break;
	case PRIME_POWER:
		xout = fft_prime_power(R, xin, N, opts);
		break;
	case PRIME_FACTOR:
		xout = fft_prime_factor(R, N / R, xin, opts);
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

