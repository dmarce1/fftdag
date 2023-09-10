#include "dag.hpp"
#include "fft.hpp"
#include "util.hpp"
#include "pebble.hpp"
#include "convolve.hpp"
#include "polynomial.hpp"
#include "instructions.hpp"

#include <time.h>

constexpr int Nmin = 2;
constexpr int Nmax = 32;
void test_poly();
//#define USE_DCT
double rand1() {
	return (rand() + 0.5) / RAND_MAX;
}

#include <fenv.h>

int inverse(int a, int n) {
	int t = 0;
	int newt = 1;
	int r = n;
	int newr = a;
	while (newr != 0) {
		int quotient = r / newr;
		int t0 = newt;
		int r0 = newr;
		newt = t - quotient * newt;
		newr = r - quotient * newr;
		t = t0;
		r = r0;
		if (t < 0) {
			t += n;
		}
		printf("%i\n", t);
	}
	return t;
}

int gcd(int a, int b) {
	while (a != 0 && b != 0) {
		if (a > b) {
			a = a - b;
		} else {
			b = b - a;
		}
	}
	return std::max(a, b);
}

std::string apply_header(std::string code, std::string name) {
	std::string rc;
	rc += std::string("               .global        ") + name + "\n";
	rc += std::string("               .text        ");
	rc += "\n" + name + ":\n";
	rc += code;
	return rc;
}

int main(int argc, char **argv) {
	FILE* fp = fopen("sfft.hpp", "wt");
	fprintf(fp, "#pragma once\n\n");
	fprintf(fp, "#define FFT_NMAX %i\n", Nmax);
	fprintf(fp, "#define FFT_NMIN %i\n\n", Nmin);
	fprintf(fp, "#include <cassert>\n");
	fprintf(fp, "#include <vector>\n\n");
	fprintf(fp, "extern \"C\" {\n");
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	for (int width = 1; width <= 4; width *= 2) {
		set_simd_width(width);
		int cnt1 = 0;
		int cnt2 = 0;
		fft_reset_cache();
		fprintf(stderr, "------------------------------REAL/COMPLEX-------------------------\n");
		for (int N = Nmin; N <= Nmax; N++) {
			auto inputs = math_vertex::new_inputs(2 * N);
			auto outputs = fft(inputs, N, 0);
			auto cnt = math_vertex::operation_count(outputs);
			auto tmp = math_vertex::execute_all(std::move(inputs), outputs, true, true, DIT);
			fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, 0).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
			std::string code = apply_header(tmp.first, std::string("sfft_real_complex_w") + std::to_string(width) + "_" + std::to_string(N));
			std::string fname = "fft.real.complex.w" + std::to_string(width) + "." + std::to_string(N) + ".S";
			FILE* fp = fopen(fname.c_str(), "wt");
			fprintf(fp, "%s\n", code.c_str());
			fclose(fp);
			cnt1 += tmp.second;
			cnt2 += cnt.total();
		}

		fprintf(stderr, "------------------------------SKEW-------------------------\n");
		for (int N = Nmin; N <= Nmax; N++) {
			auto inputs = math_vertex::new_inputs(N);
			auto outputs = fft(inputs, 2 * N, FFT_SKEW);
			auto cnt = math_vertex::operation_count(outputs);
			auto tmp = math_vertex::execute_all(std::move(inputs), outputs, false, false, NONE);
			fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_SKEW).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
			std::string code = apply_header(tmp.first, std::string("sfft_skew_w") + std::to_string(width) + "_" + std::to_string(N));
			std::string fname = "fft.skew.w" + std::to_string(width) + "." + std::to_string(N) + ".S";
			FILE* fp = fopen(fname.c_str(), "wt");
			fprintf(fp, "%s\n", code.c_str());
			fclose(fp);
			cnt1 += tmp.second;
			cnt2 += cnt.total();
		}
		fft_reset_cache();
		fprintf(stderr, "------------------------------REAL-------------------------\n");
		for (int N = Nmin; N <= Nmax; N++) {
			auto inputs = math_vertex::new_inputs(N);
			auto outputs = fft(inputs, N, FFT_REAL);
			auto cnt = math_vertex::operation_count(outputs);
			auto tmp = math_vertex::execute_all(std::move(inputs), outputs, false, false, NONE);
			fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_REAL).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
			std::string code = apply_header(tmp.first, std::string("sfft_real_w") + std::to_string(width) + "_" + std::to_string(N));
			std::string fname = "fft.real.w" + std::to_string(width) + "." + std::to_string(N) + ".S";
			FILE* fp = fopen(fname.c_str(), "wt");
			fprintf(fp, "%s\n", code.c_str());
			fclose(fp);
			cnt1 += tmp.second;
			cnt2 += cnt.total();
		}

		fprintf(stderr, "------------------------------COMPLEX-------------------------\n");
		for (int N = Nmin; N <= Nmax; N++) {
			auto inputs = math_vertex::new_inputs(2 * N);
			auto outputs = fft(inputs, N, 0);
			auto cnt = math_vertex::operation_count(outputs);
			auto tmp = math_vertex::execute_all(std::move(inputs), outputs, true, false, NONE);
			fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, 0).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
			std::string code = apply_header(tmp.first, std::string("sfft_complex_w") + std::to_string(width) + "_" + std::to_string(N));
			std::string fname = "fft.complex.w" + std::to_string(width) + "." + std::to_string(N) + ".S";
			FILE* fp = fopen(fname.c_str(), "wt");
			fprintf(fp, "%s\n", code.c_str());
			fclose(fp);
			cnt1 += tmp.second;
			cnt2 += cnt.total();
		}

		fprintf(stderr, "------------------------------COMPLEX-DIT_---------------------\n");
		for (int N = Nmin; N <= Nmax; N++) {
			auto inputs = math_vertex::new_inputs(2 * N);
			auto outputs = fft(inputs, N, 0);
			auto cnt = math_vertex::operation_count(outputs);
			auto tmp = math_vertex::execute_all(std::move(inputs), outputs, true, false, DIT);
			fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, 0).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
			std::string code = apply_header(tmp.first, std::string("sfft_complex_dit_w") + std::to_string(width) + "_" + std::to_string(N));
			std::string fname = "fft.complex.dit.w" + std::to_string(width) + "." + std::to_string(N) + ".S";
			FILE* fp = fopen(fname.c_str(), "wt");
			fprintf(fp, "%s\n", code.c_str());
			fclose(fp);
			cnt1 += tmp.second;
			cnt2 += cnt.total();
		}

		fprintf(stderr, "------------------------------COMPLEX-DIF_---------------------\n");
		for (int N = Nmin; N <= Nmax; N++) {
			auto inputs = math_vertex::new_inputs(2 * N);
			auto outputs = fft(inputs, N, 0);
			auto cnt = math_vertex::operation_count(outputs);
			auto tmp = math_vertex::execute_all(std::move(inputs), outputs, true, false, DIF);
			fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, 0).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
			std::string code = apply_header(tmp.first, std::string("sfft_complex_dif_w") + std::to_string(width) + "_" + std::to_string(N));
			std::string fname = "fft.complex.dif.w" + std::to_string(width) + "." + std::to_string(N) + ".S";
			FILE* fp = fopen(fname.c_str(), "wt");
			fprintf(fp, "%s\n", code.c_str());
			fclose(fp);
			cnt1 += tmp.second;
			cnt2 += cnt.total();
		}

		/*fft_reset_cache();
		 fprintf( stderr, "------------------------------REAL--------------------------------\n");
		 for (int N = Nmin; N <= Nmax; N++) {
		 auto inputs = math_vertex::new_inputs(N);
		 auto outputs = fft(inputs, N, FFT_REAL);
		 auto cnt = math_vertex::operation_count(outputs);
		 auto tmp = math_vertex::execute_all(std::move(inputs), outputs, false, false, NONE);
		 fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_REAL).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		 std::string code = apply_header(tmp.first, std::string("sfft_real_w") + std::to_string(width) + "_" + std::to_string(N));
		 std::string fname = "fft.real.w" + std::to_string(width) + "." std::to_string(N) + ".S";
		 FILE* fp = fopen(fname.c_str(), "wt");
		 fprintf(fp, "%s\n", code.c_str());
		 fclose(fp);
		 cnt1 += tmp.second;
		 cnt2 += cnt.total();
		 }*/
		fft_reset_cache();

		fprintf(fp, "void sfft_perf_shuf_w%i(double* x, double* y);", width);
		fprintf(fp, "void sfft_inv_perf_shuf_w%i(double* x, double* y);", width);

		/*	for (int n = FFT_NMIN; n <= FFT_NMAX; n++) {
		 fprintf(fp, "void sfft_complex_w%i_%i(double*, double*, int, int);\n", width, n);
		 fprintf(fp, "void sfft_complex_dit_w%i_%i(double*, double*, double*, int, int);\n", width, n);
		 fprintf(fp, "void sfft_complex_dif_w%i_%i(double*, double*, double*, int, int);\n", width, n);
		 }*
		 for (int n = FFT_NMIN; n <= FFT_NMAX; n++) {
		 //	fprintf(fp, "void sfft_real_w%i_%i(double*, double*);\n", width, n);
		 }
		 /*                                         rdi        rsi        rdx    rcx    r8     r9 */
		fprintf(fp, "void sfft_real_w%i(double* x, int s, int N);\n", width);
		fprintf(fp, "void sfft_skew_w%i(double* x, int s, int N);\n", width);
		fprintf(fp, "void sfft_complex_w%i(double* x, double* y, int s, int N);\n", width);
		fprintf(fp, "void sfft_real_complex_w%i(double* x, double* y, int s, int N, double* w);\n", width);
		fprintf(fp, "void sfft_complex_dit_w%i(double* x, double* y, int s, int N, double* w);\n", width);
		fprintf(fp, "void sfft_complex_dif_w%i(double* x, double* y, int s, int N, double* w);\n", width);
	}
	fprintf(fp, "}\n");
	fclose(fp);

	fft_reset_cache();
	system("cp ../../gen_src/main.cpp .\n");
	system("cp ../../gen_src/test.cpp .\n");
	system("cp ../../gen_src/load.S .\n");
	system("cp ../../gen_src/store.S .\n");
	system("cp ../../gen_src/sfft_complex.S .\n");
	system("cp ../../gen_src/sfft_real.S .\n");
	system("cp ../../gen_src/sfft_real_complex.S .\n");
	system("cp ../../gen_src/sfft_skew.S .\n");
	system("cp ../../gen_src/twiddle.S .\n");
	system("cp ../../gen_src/shuffle.S .\n");
	system("cp ../../gen_src/util.cpp .\n");
	system("cp ../../gen_src/util.hpp .\n");
	system("cp ../../gen_src/types.hpp .\n");

	/* fprintf(fp, "\n");
	 fprintf(fp, "inline void sfft_complex(double* xr, double* xi, size_t sr, size_t si, size_t N) {\n");
	 fprintf(fp, "\tswitch(N) {\n");
	 for (int n = FFT_NMIN; n <= FFT_NMAX; n++) {
	 fprintf(fp, "\tcase %i:\n", n);
	 fprintf(fp, "\t\tsfft_complex_%i(xr, xi, sr, si);\n", n);
	 fprintf(fp, "\t\tbreak;\n");
	 }
	 fprintf(fp, "\tdefault:\n");
	 fprintf(fp, "\t\tassert(false);\n");
	 fprintf(fp, "\t}\n");
	 fprintf(fp, "}\n");
	 fprintf(fp, "\n");
	 fprintf(fp, "\n");
	 fprintf(fp, "inline void sfft_complex_dit(double* xr, double* xi, double* wr, double* wi, size_t sr, size_t si, size_t N) {\n");
	 fprintf(fp, "\tswitch(N) {\n");
	 for (int n = FFT_NMIN; n <= FFT_NMAX; n++) {
	 fprintf(fp, "\tcase %i:\n", n);
	 fprintf(fp, "\t\tsfft_complex_dit_%i(xr, xi, wr, wi, sr, si);\n", n);
	 fprintf(fp, "\t\tbreak;\n");
	 }
	 fprintf(fp, "\tdefault:\n");
	 fprintf(fp, "\t\tassert(false);\n");
	 fprintf(fp, "\t}\n");
	 fprintf(fp, "}\n");
	 fprintf(fp, "\n");
	 fprintf(fp, "\n");
	 fprintf(fp, "inline void sfft_complex_dif(double* xr, double* xi, double* wr, double* wi, size_t sr, size_t si, size_t N) {\n");
	 fprintf(fp, "\tswitch(N) {\n");
	 for (int n = FFT_NMIN; n <= FFT_NMAX; n++) {
	 fprintf(fp, "\tcase %i:\n", n);
	 fprintf(fp, "\t\tsfft_complex_dif_%i(xr, xi, wr, wi, sr, si);\n", n);
	 fprintf(fp, "\t\tbreak;\n");
	 }
	 fprintf(fp, "\tdefault:\n");
	 fprintf(fp, "\t\tassert(false);\n");
	 fprintf(fp, "\t}\n");
	 fprintf(fp, "}\n");
	 fprintf(fp, "\n");
	 fprintf(fp, "\n");
	 fprintf(fp, "inline void sfft_real(double* x, size_t s, size_t N) {\n");
	 fprintf(fp, "\tswitch(N) {\n");
	 for (int n = FFT_NMIN; n <= FFT_NMAX; n++) {
	 fprintf(fp, "\tcase %i:\n", n);
	 fprintf(fp, "\t\tsfft_real_%i(x, s);\n", n);
	 fprintf(fp, "\t\tbreak;\n");
	 }
	 fprintf(fp, "\tdefault:\n");
	 fprintf(fp, "\t\tassert(false);\n");
	 fprintf(fp, "\t}\n");
	 fprintf(fp, "}\n");
	 fprintf(fp, "\n");
	 fprintf(fp, "inline void sfft_dct2(double* x, size_t s, size_t N) {\n");
	 fprintf(fp, "\tswitch(N) {\n");
	 for (int n = FFT_NMIN; n <= FFT_NMAX; n++) {
	 fprintf(fp, "\tcase %i:\n", n);
	 fprintf(fp, "\t\tsfft_dct2_%i(x, s);\n", n);
	 fprintf(fp, "\t\tbreak;\n");
	 }
	 fprintf(fp, "\tdefault:\n");
	 fprintf(fp, "\t\tassert(false);\n");
	 fprintf(fp, "\t}\n");
	 fprintf(fp, "}\n");
	 fprintf(fp, "\n");
	 fclose(fp);*/

	fp = fopen("Makefile", "wt");
	fprintf(fp, "CC=g++\n");
	fprintf(fp, "CFLAGS=-I. -march=native\n");
	//fprintf(fp, "CFLAGS=-I. -g -fsanitize=address -D_GLIBCXX_DEBUG -march=native\n");
	fprintf(fp, "DEPS = sfft.hpp\n");
	fprintf(fp, "OBJ = ");
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.skew.w1.%i.o ", n);
		fprintf(fp, "fft.skew.w2.%i.o ", n);
		fprintf(fp, "fft.skew.w4.%i.o ", n);
		fprintf(fp, "fft.real.w1.%i.o ", n);
		fprintf(fp, "fft.real.w2.%i.o ", n);
		fprintf(fp, "fft.real.w4.%i.o ", n);
		fprintf(fp, "fft.complex.w1.%i.o ", n);
		fprintf(fp, "fft.complex.w2.%i.o ", n);
		fprintf(fp, "fft.complex.w4.%i.o ", n);
		fprintf(fp, "fft.real.complex.w1.%i.o ", n);
		fprintf(fp, "fft.real.complex.w2.%i.o ", n);
		fprintf(fp, "fft.real.complex.w4.%i.o ", n);
		fprintf(fp, "fft.complex.dit.w1.%i.o ", n);
		fprintf(fp, "fft.complex.dit.w2.%i.o ", n);
		fprintf(fp, "fft.complex.dit.w4.%i.o ", n);
		fprintf(fp, "fft.complex.dif.w1.%i.o ", n);
		fprintf(fp, "fft.complex.dif.w2.%i.o ", n);
		fprintf(fp, "fft.complex.dif.w4.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		//	fprintf(fp, "fft.skew.%i.o ", n);
	}
	fprintf(fp, "\n%%.o: %%.S $(DEPS)\n");
	fprintf(fp, "\t$(CC) -c -o $@ $< $(CFLAGS)\n\n");
	fprintf(fp, "\n%%.o: %%.cpp $(DEPS)\n");
	fprintf(fp, "\t$(CC) -c -o $@ $< $(CFLAGS)\n\n");
	fprintf(fp, "ffttest: $(OBJ) util.o test.o load.o sfft_skew.o sfft_complex.o sfft_real.o store.o twiddle.o\n");
	fprintf(fp, "\t$(CC) -o $@ $^ $(CFLAGS) -lfftw3\n");
	fprintf(fp, "libsfft.a: shuffle.o sfft_real.o sfft_complex.o sfft_skew.o $(OBJ)\n");
	fprintf(fp, "\tar -rcs libsfft.a shuffle.o sfft_complex.o sfft_real.o sfft_skew.o $(OBJ)\n");
	fclose(fp);
	print_fft_bests();
	return 0;
}
