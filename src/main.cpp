#include "dag.hpp"
#include "fft.hpp"
#include "util.hpp"
#include "pebble.hpp"
#include "convolve.hpp"
#include "polynomial.hpp"

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
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	int cnt1 = 0;
	int cnt2 = 0;

	fft_reset_cache();
	fprintf( stderr, "------------------------------COMPLEX-----------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(2 * N);
		auto outputs = fft(inputs, N, 0);
		auto cnt = math_vertex::operation_count(outputs);
		auto tmp = math_vertex::execute_all(std::move(inputs), outputs, true, 4);
		fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, 0).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		std::string code = apply_header(tmp.first, std::string("short_fft_complex_simd4_v4df") + std::to_string(N));
		std::string fname = "fft.complex." + std::to_string(N) + ".S";
		FILE* fp = fopen(fname.c_str(), "wt");
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
		cnt1 += tmp.second;
		cnt2 += cnt.total();
	}

	fft_reset_cache();
	fprintf( stderr, "------------------------------REAL--------------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, N, FFT_REAL);
		auto cnt = math_vertex::operation_count(outputs);
		auto tmp = math_vertex::execute_all(std::move(inputs), outputs, false, 4);
		fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_REAL).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		std::string code = apply_header(tmp.first, std::string("short_fft_real_simd4_v4df") + std::to_string(N));
		std::string fname = "fft.real." + std::to_string(N) + ".S";
		FILE* fp = fopen(fname.c_str(), "wt");
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
		cnt1 += tmp.second;
		cnt2 += cnt.total();
	}
	fft_reset_cache();
	/*	fprintf( stderr, "------------------------------REAL INVERSE-------------------------\n");
	 for (int N = Nmin; N <= Nmax; N++) {
	 auto inputs = math_vertex::new_inputs(N);
	 auto outputs = fft(inputs, N, FFT_REAL | FFT_INV);
	 auto cnt = math_vertex::operation_count(outputs);
	 auto tmp = math_vertex::execute_all(std::move(inputs), outputs);
	 fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_REAL | FFT_INV).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
	 std::string code;
	 code += std::string("               .global        ") + "fft_kernel_complex_" + std::to_string(N) + "\n";
	 code += "\nfft_kernel_complex_" + std::to_string(N) + ":\n";
	 code += tmp.first;
	 std::string fname = "fft.real_inv." + std::to_string(N) + ".S";
	 FILE* fp = fopen(fname.c_str(), "wt");
	 code = "#include \"types.hpp\"\n\n";
	 fprintf(fp, "%s\n", code.c_str());
	 fclose(fp);
	 cnt1 += tmp.second;
	 cnt2 += cnt.total();
	 }*/
	fprintf( stderr, "O: %i D: %i\n", cnt2, cnt1);
	system("cp ../../gen_src/main.cpp .\n");
	system("cp ../../gen_src/test.cpp .\n");
	system("cp ../../gen_src/util.cpp .\n");
	system("cp ../../gen_src/util.hpp .\n");
	system("cp ../../gen_src/types.hpp .\n");

	FILE* fp = fopen("fft.hpp", "wt");
	fprintf(fp, "#pragma once\n\n");
	fprintf(fp, "#define FFT_NMAX %i\n", Nmax);
	fprintf(fp, "#define FFT_NMIN %i\n\n\n", Nmin);
	fprintf(fp, "extern \"C\" {");

	for( int n = FFT_NMIN; n <= FFT_NMAX; n++) {
		fprintf(fp, "void short_fft_complex_v4df_%i(v4df*, v4df*, v4df*, size_t, size_t);\n", n);
	}
	for( int n = FFT_NMIN; n <= FFT_NMAX; n++) {
		fprintf(fp, "void short_fft_real_v4df_%i(v4df*, v4df*, size_t);\n", n);
	}
	fprintf(fp, "}\n");
	fclose(fp);

	fp = fopen("Makefile", "wt");
	fprintf(fp, "CC=g++\n");
	fprintf(fp, "CFLAGS=-I. -Ofast -march=native\n");
//	fprintf(fp, "CFLAGS=-I. -g -fsanitize=address -D_GLIBCXX_DEBUG -march=native\n");
	fprintf(fp, "DEPS = fft.hpp\n");
//	fprintf(fp, "OBJ = fft.o util.o ");
	fprintf(fp, "OBJ = ");
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.complex.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.real.%i.o ", n);
	}
	fprintf(fp, "\n%%.o: %%.cpp $(DEPS)\n");
	fprintf(fp, "\t$(CC) -c -o $@ $< $(CFLAGS)\n\n");
//	fprintf(fp, "ffttest: $(OBJ) test.o\n");
//	fprintf(fp, "\t$(CC) -o $@ $^ $(CFLAGS) -lfftw3\n");
	fprintf(fp, "sfftlib.a: $(OBJ)\n");
	fprintf(fp, "\tar -rcs sfftlib.a $(OBJ)\n");
	fclose(fp);
	print_fft_bests();
	return 0;
}
