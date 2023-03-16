#include "fft.hpp"
#include <fftw3.h>
#include <vector>
#include <complex>
#include <chrono>
#include <unordered_map>
#include <fenv.h>

class timer {
	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
	double time;
public:
	inline timer() {
		time = 0.0;
	}
	inline void stop() {
		std::chrono::time_point<std::chrono::high_resolution_clock> stop_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> dur = stop_time - start_time;
		time += dur.count();
	}
	inline void start() {
		start_time = std::chrono::high_resolution_clock::now();
	}
	inline void reset() {
		time = 0.0;
	}
	inline double read() {
		return time;
	}
};

void fftw(std::vector<std::complex<double>>& x) {
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, fftw_complex*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);
		plans[N] = fftw_plan_dft_1d(N, in[N], out[N], FFTW_FORWARD, FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n][0] = x[n].real();
		i[n][1] = x[n].imag();
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		x[n].real(o[n][0]);
		x[n].imag(o[n][1]);
	}

}

void fftw_dct(std::vector<double>& x) {
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, double*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (double*) malloc(sizeof(double) * N);
		out[N] = (double*) malloc(sizeof(double) * N);
		plans[N] = fftw_plan_r2r_1d(N, in[N], out[N], FFTW_REDFT10, FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = x[n];
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		x[n] = (o[n]);
	}

}

void fftw_dst(std::vector<double>& x) {
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, double*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (double*) malloc(sizeof(double) * N);
		out[N] = (double*) malloc(sizeof(double) * N);
		plans[N] = fftw_plan_r2r_1d(N, in[N], out[N], FFTW_RODFT10, FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = x[n];
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		x[n] = (o[n]);
	}

}

void fftw_real(std::vector<std::complex<double>>& xout, const std::vector<double>& xin) {
	const int N = xin.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (double*) malloc(sizeof(double) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * (N / 2 + 1));
		plans[N] = fftw_plan_dft_r2c_1d(N, in[N], out[N], FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = xin[n];
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N / 2 + 1; n++) {
		xout[n].real(o[n][0]);
		xout[n].imag(o[n][1]);
	}

}

void fftw_real_inv(const std::vector<std::complex<double>>& xout, std::vector<double>& xin) {
	const int N = xin.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (double*) malloc(sizeof(double) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * (N / 2 + 1));
		plans[N] = fftw_plan_dft_c2r_1d(N, out[N], in[N], FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N / 2 + 1; n++) {
		o[n][0] = xout[n].real();
		o[n][1] = xout[n].imag();
	}
	o[0][1] = 0;
	if (N % 2 == 0) {
		o[N / 2][1] = 0;
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		xin[n] = i[n];
	}

}

double rand1() {
	return (rand() + 0.5) / RAND_MAX;
}

void test() {
	timer tm3, tm4;

	tm3.reset();
	tm4.reset();
	printf("\ncomplex\n");
	for (int N = FFT_NMIN; N <= FFT_NMAX; N += 1) {
		timer tm1, tm2;
		double err;
		double max;
		for (int ti = 0; ti < 256; ti++) {
			err = 0.0;
			max = 0.0;
			std::vector<std::complex<double>> X(N);
			std::vector<std::complex<double>> Y(N);
			for (int n = 0; n < N; n++) {
				X[n] = std::complex<double>(rand1(), rand1());
			}
			auto X0 = X;
			Y = X;
			tm1.start();
			tm3.start();
			fft((double*) X.data(), N);
			tm1.stop();
			tm3.stop();
			tm2.start();
			tm4.start();
			fftw(Y);
			tm2.stop();
			tm4.stop();
			for (int i = 0; i < X.size(); i++) {
				Y[i] -= X[i];
			}
			for (int n = 0; n < N; n++) {
				err += std::abs(Y[n]) * std::abs(Y[n]);
				max = std::max(max, std::abs(X0[n]));
				//		printf("%i %e %e %e %e\n", n, X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag());
			}
			err = sqrt(err / N) / max;
		}
		printf("%4i %e %e %e %e %e %e %e\n", N, err, tm1.read(), tm2.read(), tm1.read() / tm2.read(), tm3.read(), tm4.read(), tm3.read() / tm4.read());
	}
}

int main() {
	feenableexcept (FE_DIVBYZERO);
	feenableexcept (FE_OVERFLOW);
	feenableexcept (FE_INVALID);
	test();
	printf("\n");
	test();
	printf("EXITING\n");
}

