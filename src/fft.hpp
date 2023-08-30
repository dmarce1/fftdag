#pragma once

#include "math.hpp"

#define FFT_REAL (1 << 0)
#define FFT_INV  (1 << 1)
#define FFT_DCT1 (1 << 2)
#define FFT_DCT2 (1 << 3)
#define FFT_DCT3 (1 << 4)
#define FFT_DCT4 (1 << 5)
#define FFT_DST1 (1 << 6)
#define FFT_DST2 (1 << 7)
#define FFT_DST3 (1 << 8)
#define FFT_DST4 (1 << 9)
#define FFT_ODD  (1 << 10)
#define FFT_DIT  (1 << 11)
#define FFT_DIF  (1 << 12)

std::vector<math_vertex> fft(std::vector<math_vertex> xin, int N, int opts);
std::vector<cmplx> fft_modsplit(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft(std::vector<cmplx> xin, int N, int opts, bool = false);
std::vector<cmplx> fft_cooley_tukey(int N1, int N2, std::vector<cmplx> xin, int opts);
std::vector<std::complex<double>> fft_bruun(std::vector<std::complex<double>> xin, int N);
std::vector<cmplx> fft_bruun(std::vector<cmplx> xin, int N, int opts);
std::string method_name(int i);
std::string get_best_method(int N, int opts);
std::vector<cmplx> fft_raders_fast(std::vector<cmplx> xin, int N, int opts);
std::vector<cmplx> fft_raders_fft(std::vector<cmplx> xin, int N, bool padded, int opts);
void fft_reset_cache();

void fft_reset();
void print_fft_bests();
bool can_agarwal(int N);

inline std::vector<int> fft_input_signature(std::vector<cmplx> xin) {
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
			for (int n = 0; n < sig.size(); n++) {
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

struct best_x2 {
	int N;
	int opts;
	bool operator==(const best_x2& other) const {
		return N == other.N && opts == other.opts;
	}
	bool operator<(const best_x2& other) const {
		if (N < other.N) {
			return true;
		} else if (N > other.N) {
			return false;
		} else if (opts < other.opts) {
			return true;
		} else if (opts > other.opts) {
			return false;
		}
		return false;
	}
};
