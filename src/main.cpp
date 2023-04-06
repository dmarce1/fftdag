#include "dag.hpp"
#include "fft.hpp"
#include "util.hpp"
#include "pebble.hpp"
#include "convolve.hpp"
#include "polynomial.hpp"

#include <time.h>


constexpr int Nmin = 2;
constexpr int Nmax = 64;
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

int main(int argc, char **argv) {
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	int cnt1 = 0;
	int cnt2 = 0;
	/*fprintf( stderr, "------------------------------COMPLEX INVERSE-----------------------\n");
	 for (int N = Nmin; N <= Nmax; N++) {
	 auto inputs = math_vertex::new_inputs(2 * N);
	 auto outputs = fft(inputs, N, FFT_INV);
	 auto cnt = math_vertex::operation_count(outputs);
	 auto tmp = math_vertex::execute_all(std::move(inputs), outputs);
	 fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_INV).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
	 auto code = tmp.first;
	 std::string fname = "fft.complex_inv." + std::to_string(N) + ".cpp";
	 FILE* fp = fopen(fname.c_str(), "wt");
	 code = std::string("\n\nvoid fft_complex_inv_") + std::to_string(N) + "(fft_type* x) {\n" + code + "}\n\n";
	 code = "#include \"type.hpp\"\n" + code;
	 fprintf(fp, "%s\n", code.c_str());
	 fclose(fp);
	 cnt1 += tmp.second;
	 cnt2 += cnt.total();
	 }

	fft_reset_cache();
	fprintf( stderr, "------------------------------REDUCED-----------------------------\n");
	for (int N = Nmin; N <= 2 * Nmax; N += 2) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, N, FFT_ODD);
		auto cnt = math_vertex::operation_count(outputs);
		auto tmp = math_vertex::execute_all(std::move(inputs), outputs);
		fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_ODD).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		auto code = tmp.first;
		std::string fname = "fft.reduced." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_reduced_") + std::to_string(N) + "(fft_type* x) {\n" + code + "}\n\n";
	 	code = "#include \"type.hpp\"\n" + code;
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
		cnt1 += tmp.second;
		cnt2 += cnt.total();
	}*/

	fft_reset_cache();
	fprintf( stderr, "------------------------------REAL INVERSE-------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, N, FFT_REAL | FFT_INV);
		auto cnt = math_vertex::operation_count(outputs);
		auto tmp = math_vertex::execute_all(std::move(inputs), outputs);
		fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_REAL | FFT_INV).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		auto code = tmp.first;
		std::string fname = "fft.real_inv." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_real_inv_") + std::to_string(N) + "(fft_type* x) {\n" + code + "}\n\n";
		code = "#include \"type.hpp\"\n" + code;
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
		cnt1 += tmp.second;
		cnt2 += cnt.total();
	}

	fft_reset_cache();
	fprintf( stderr, "------------------------------COMPLEX-----------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(2 * N);
		auto outputs = fft(inputs, N, 0);
		auto cnt = math_vertex::operation_count(outputs);
		auto tmp = math_vertex::execute_all(std::move(inputs), outputs);
		fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, 0).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		auto code = tmp.first;
		std::string fname = "fft.complex." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_complex_") + std::to_string(N) + "(fft_type* x) {\n" + code + "}\n\n";
		code = "#include \"type.hpp\"\n" + code;
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
		auto tmp = math_vertex::execute_all(std::move(inputs), outputs);
		fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_REAL).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		auto code = tmp.first;
		std::string fname = "fft.real." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_real_") + std::to_string(N) + "(fft_type* x) {\n" + code + "}\n\n";
		code = "#include \"type.hpp\"\n" + code;
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
		cnt1 += tmp.second;
		cnt2 += cnt.total();
	}
#ifdef USE_DCT
	fprintf( stderr, "------------------------------DCT-I--------------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, 2 * (N - 1), FFT_DCT1);
		auto cnt = math_vertex::operation_count(outputs);
		auto tmp = math_vertex::execute_all(std::move(inputs), outputs);
		fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_DCT1).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		auto code = tmp.first;
		std::string fname = "fft.dct1." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_dct1_") + std::to_string(N) + "(fft_type* x) {\n" + code + "}\n\n";
		code = "#include \"type.hpp\"\n" + code;
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
		cnt1 += tmp.second;
		cnt2 += cnt.total();
	}
	fprintf( stderr, "------------------------------DCT-II--------------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, 4 * N, FFT_DCT2);
		auto cnt = math_vertex::operation_count(outputs);
		auto tmp = math_vertex::execute_all(std::move(inputs), outputs);
		fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_DCT2).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		auto code = tmp.first;
		std::string fname = "fft.dct2." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_dct2_") + std::to_string(N) + "(fft_type* x) {\n" + code + "}\n\n";
		code = "#include \"type.hpp\"\n" + code;
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
		cnt1 += tmp.second;
		cnt2 += cnt.total();
	}

	fprintf( stderr, "------------------------------DCT-III--------------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, 4 * N, FFT_DCT3);
		auto cnt = math_vertex::operation_count(outputs);
		auto tmp = math_vertex::execute_all(std::move(inputs), outputs);
		fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_DCT3).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		auto code = tmp.first;
		std::string fname = "fft.dct3." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_dct3_") + std::to_string(N) + "(fft_type* x) {\n" + code + "}\n\n";
		code = "#include \"type.hpp\"\n" + code;
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
		cnt1 += tmp.second;
		cnt2 += cnt.total();
	}
	fprintf( stderr, "------------------------------DCT-IV--------------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, 8 * N, FFT_DCT4);
		auto cnt = math_vertex::operation_count(outputs);
		auto tmp = math_vertex::execute_all(std::move(inputs), outputs);
		fprintf(stderr, "N = %4i | %16s | tot = %4i | add = %4i | mul = %4i | neg = %4i | decls = %i\n", N, get_best_method(N, FFT_DCT4).c_str(), cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg, tmp.second);
		auto code = tmp.first;
		std::string fname = "fft.dct4." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_dct4_") + std::to_string(N) + "(fft_type* x) {\n" + code + "}\n\n";
		code = "#include \"type.hpp\"\n" + code;
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
		cnt1 += tmp.second;
		cnt2 += cnt.total();
	}
#endif
	/*	fprintf( stderr, "------------------------------CONVOLVE-----------------------------\n");
	 for (int N = Nmin; N <= std::min(Nmax, 32); N++) {
	 std::vector<std::complex<fft_type>> h(N);
	 srand(42);
	 for (int n = 0; n < N; n++) {
	 h[n] = rand() / (RAND_MAX + 0.5);
	 }
	 auto x = math_vertex::new_inputs(2 * N);
	 std::vector<cmplx> X(N);
	 for (int n = 0; n < N; n++) {
	 X[n].x = x[2 * n];
	 X[n].y = x[2 * n + 1];
	 }
	 auto Z = convolve_fft(X, h);
	 auto Y = convolve_fast(X, h);
	 std::vector<math_vertex> y;
	 for (int n = 0; n < N; n++) {
	 y.push_back(Y[n].x);
	 y.push_back(Y[n].y);
	 }
	 auto cnt = math_vertex::operation_count(y);
	 auto cnt3 = math_vertex::operation_count(Z);
	 if (cnt.total() == 0 || cnt.total() > cnt3.total()) {
	 Y = std::move(Z);
	 }
	 if (cnt.total() == 0) {
	 cnt.add = cnt.mul = cnt.neg = 0;
	 }
	 auto tmp = math_vertex::execute_all(std::move(x), y);
	 fprintf(stderr, "N = %4i%2s | tot = %4i(%4i) | add = %4i(%4i) | mul = %4i(%4i) | neg = %4i(%4i) | decls = %i\n", N, cnt.total() < cnt3.total() ? "^*" : "", cnt.total(), cnt3.total(), cnt.add, cnt3.add, cnt.mul, cnt3.mul, cnt.neg, cnt3.neg,
	 tmp.second);
	 auto code = tmp.first;
	 std::string fname = "convolve." + std::to_string(N) + ".cpp";
	 FILE* fp = fopen(fname.c_str(), "wt");
	 code = std::string("\n\nvoid convolve_") + std::to_string(N) + "(fft_type* x) {\n" + code + "}\n\n";
	 fprintf(fp, "%s\n", code.c_str());
	 fclose(fp);
	 cnt1 += tmp.second;
	 cnt2 += cnt.total();
	 }*/
	fprintf( stderr, "O: %i D: %i\n", cnt2, cnt1);
	system("cp ../../gen_src/main.cpp .\n");
	system("cp ../../gen_src/util.cpp .\n");
	system("cp ../../gen_src/util.hpp .\n");
	system("cp ../../gen_src/type.hpp .\n");

	FILE* fp = fopen("fft.hpp", "wt");
	fprintf(fp, "#pragma once\n\n");
	fprintf(fp, "#include \"type.hpp\"\n\n\n");
	fprintf(fp, "#define FFT_NMAX %i\n", Nmax);
	fprintf(fp, "#define FFT_NMIN %i\n\n\n", Nmin);
	/*	for (int N = Nmin; N <= Nmax; N++) {
	 fprintf(fp, "void convolve_%i(fft_type*);\n", N);
	 }*/
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_complex_%i(fft_type*);\n", N);
	}
//	for (int N = Nmin; N <= 2 * Nmax; N += 2) {
//		fprintf(fp, "void fft_reduced_%i(fft_type*);\n", N);
//	}
//	for (int N = Nmin; N <= Nmax; N++) {
//		fprintf(fp, "void fft_complex_inv_%i(fft_type*);\n", N);
//	}
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_real_%i(fft_type*);\n", N);
	}
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_real_inv_%i(fft_type*);\n", N);
	}
	//fprintf(fp, "\n\nvoid convolve(fft_type*, int);\n");
	fprintf(fp, "\n\nvoid fft_complex(fft_type*, int);\n");
//	fprintf(fp, "\n\nvoid fft_reduced(fft_type*, int);\n");
//	fprintf(fp, "\n\nvoid fft_complex_inv(fft_type*, int);\n");
	fprintf(fp, "\n\nvoid fft_real(fft_type*, int);\n");
	fprintf(fp, "\n\nvoid fft_real_inv(fft_type*, int);\n");
#ifdef USE_DCT
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_dct1_%i(fft_type*);\n", N);
	}
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_dct2_%i(fft_type*);\n", N);
	}
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_dct3_%i(fft_type*);\n", N);
	}
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_dct4_%i(fft_type*);\n", N);
	}
	fprintf(fp, "\n\nvoid fft_dct1(fft_type*, int);\n");
	fprintf(fp, "\n\nvoid fft_dct2(fft_type*, int);\n");
	fprintf(fp, "\n\nvoid fft_dct3(fft_type*, int);\n");
	fprintf(fp, "\n\nvoid fft_dct4(fft_type*, int);\n");
	fprintf(fp, "\n\n");
#endif
	fclose(fp);

	fp = fopen("fft.cpp", "wt");
	fprintf(fp, "#include \"fft.hpp\"\n\n\n");
	fprintf(fp, "using fft_func_type = void(*)(fft_type*);\n\n");
	/*	fprintf(fp, "fft_func_type convolve_pointer[FFT_NMAX + 1] = {");
	 for (int N = 0; N <= Nmax; N++) {
	 if (N % 8 == 0) {
	 fprintf(fp, "\n\t");
	 }
	 if (N < Nmin) {
	 fprintf(fp, "nullptr");
	 } else {
	 fprintf(fp, "&convolve_%i", N);
	 }
	 if (N != Nmax) {
	 fprintf(fp, ", ");
	 }
	 }
	 fprintf(fp, "\n};\n\n");*/
	fprintf(fp, "fft_func_type fft_complex_pointer[FFT_NMAX + 1] = {");
	for (int N = 0; N <= Nmax; N++) {
		if (N % 8 == 0) {
			fprintf(fp, "\n\t");
		}
		if (N < Nmin) {
			fprintf(fp, "nullptr");
		} else {
			fprintf(fp, "&fft_complex_%i", N);
		}
		if (N != Nmax) {
			fprintf(fp, ", ");
		}
	}
	fprintf(fp, "\n};\n\n");
/*	fprintf(fp, "fft_func_type fft_reduced_pointer[2 * FFT_NMAX + 1] = {");
	for (int N = 0; N <= 2 * Nmax; N++) {
		if (N % 8 == 0) {
			fprintf(fp, "\n\t");
		}
		if (N < Nmin || N % 2 != 0) {
			fprintf(fp, "nullptr");
		} else {
			fprintf(fp, "&fft_reduced_%i", N);
		}
		if (N != 2 * Nmax) {
			fprintf(fp, ", ");
		}
	}
	fprintf(fp, "\n};\n\n");*/
	/*	fprintf(fp, "fft_func_type fft_complex_inv_pointer[FFT_NMAX + 1] = {");
	 for (int N = 0; N <= Nmax; N++) {
	 if (N % 8 == 0) {
	 fprintf(fp, "\n\t");
	 }
	 if (N < Nmin) {
	 fprintf(fp, "nullptr");
	 } else {
	 fprintf(fp, "&fft_complex_inv_%i", N);
	 }
	 if (N != Nmax) {
	 fprintf(fp, ", ");
	 }
	 }
	 fprintf(fp, "\n};\n\n");*/
	fprintf(fp, "fft_func_type fft_real_pointer[FFT_NMAX + 1] = {");
	for (int N = 0; N <= Nmax; N++) {
		if (N % 8 == 0) {
			fprintf(fp, "\n\t");
		}
		if (N < Nmin) {
			fprintf(fp, "nullptr");
		} else {
			fprintf(fp, "&fft_real_%i", N);
		}
		if (N != Nmax) {
			fprintf(fp, ", ");
		}
	}
	fprintf(fp, "\n};\n\n");
	fprintf(fp, "fft_func_type fft_real_inv_pointer[FFT_NMAX + 1] = {");
	for (int N = 0; N <= Nmax; N++) {
		if (N % 8 == 0) {
			fprintf(fp, "\n\t");
		}
		if (N < Nmin) {
			fprintf(fp, "nullptr");
		} else {
			fprintf(fp, "&fft_real_inv_%i", N);
		}
		if (N != Nmax) {
			fprintf(fp, ", ");
		}
	}
	fprintf(fp, "\n};\n\n");
	/*	fprintf(fp, "void convolve(fft_type* x, int N) {\n");
	 fprintf(fp, "\t(*convolve_pointer[N])(x);\n");
	 fprintf(fp, "}\n\n\n");*/
	fprintf(fp, "void fft_complex(fft_type* x, int N) {\n");
	fprintf(fp, "\t(*fft_complex_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
//	fprintf(fp, "void fft_reduced(fft_type* x, int N) {\n");
//	fprintf(fp, "\t(*fft_reduced_pointer[N])(x);\n");
//	fprintf(fp, "}\n\n\n");
	/*	fprintf(fp, "void fft_complex_inv(fft_type* x, int N) {\n");
	 fprintf(fp, "\t(*fft_complex_inv_pointer[N])(x);\n");
	 fprintf(fp, "}\n\n\n");*/
	fprintf(fp, "void fft_real(fft_type* x, int N) {\n");
	fprintf(fp, "\t(*fft_real_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
	fprintf(fp, "void fft_real_inv(fft_type* x, int N) {\n");
	fprintf(fp, "\t(*fft_real_inv_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
#ifdef USE_DCT
	fprintf(fp, "fft_type fft_dct1_pointer[FFT_NMAX + 1] = {");
	for (int N = 0; N <= Nmax; N++) {
		if (N % 8 == 0) {
			fprintf(fp, "\n\t");
		}
		if (N < Nmin) {
			fprintf(fp, "nullptr");
		} else {
			fprintf(fp, "&fft_dct1_%i", N);
		}
		if (N != Nmax) {
			fprintf(fp, ", ");
		}
	}
	fprintf(fp, "\n};\n\n");
	fprintf(fp, "fft_type fft_dct2_pointer[FFT_NMAX + 1] = {");
	for (int N = 0; N <= Nmax; N++) {
		if (N % 8 == 0) {
			fprintf(fp, "\n\t");
		}
		if (N < Nmin) {
			fprintf(fp, "nullptr");
		} else {
			fprintf(fp, "&fft_dct2_%i", N);
		}
		if (N != Nmax) {
			fprintf(fp, ", ");
		}
	}
	fprintf(fp, "\n};\n\n");
	fprintf(fp, "fft_type fft_dct3_pointer[FFT_NMAX + 1] = {");
	for (int N = 0; N <= Nmax; N++) {
		if (N % 8 == 0) {
			fprintf(fp, "\n\t");
		}
		if (N < Nmin) {
			fprintf(fp, "nullptr");
		} else {
			fprintf(fp, "&fft_dct3_%i", N);
		}
		if (N != Nmax) {
			fprintf(fp, ", ");
		}
	}
	fprintf(fp, "\n};\n\n");
	fprintf(fp, "fft_type fft_dct4_pointer[FFT_NMAX + 1] = {");
	for (int N = 0; N <= Nmax; N++) {
		if (N % 8 == 0) {
			fprintf(fp, "\n\t");
		}
		if (N < Nmin) {
			fprintf(fp, "nullptr");
		} else {
			fprintf(fp, "&fft_dct4_%i", N);
		}
		if (N != Nmax) {
			fprintf(fp, ", ");
		}
	}
	fprintf(fp, "\n};\n\n");
	fprintf(fp, "void fft_dct1(fft_type* x, int N) {\n");
	fprintf(fp, "\t(*fft_dct1_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
	fprintf(fp, "void fft_dct2(fft_type* x, int N) {\n");
	fprintf(fp, "\t(*fft_dct2_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
	fprintf(fp, "void fft_dct3(fft_type* x, int N) {\n");
	fprintf(fp, "\t(*fft_dct3_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
	fprintf(fp, "void fft_dct4(fft_type* x, int N) {\n");
	fprintf(fp, "\t(*fft_dct4_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
#endif

	fclose(fp);

	fp = fopen("Makefile", "wt");
	fprintf(fp, "CC=g++\n");
	fprintf(fp, "CFLAGS=-I. -Ofast -march=native\n");
//	fprintf(fp, "CFLAGS=-I. -g -D_GLIBCXX_DEBUG\n");
	fprintf(fp, "DEPS = fft.hpp\n");
	fprintf(fp, "OBJ = fft.o util.o ");
	/*	for (int n = Nmin; n <= Nmax; n++) {
	 fprintf(fp, "convolve.%i.o ", n);
	 }*/
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.complex.%i.o ", n);
	}
	//for (int n = Nmin; n <= 2 * Nmax; n += 2) {
	//	fprintf(fp, "fft.reduced.%i.o ", n);
//	}
	/*	for (int n = Nmin; n <= Nmax; n++) {
	 fprintf(fp, "fft.complex_inv.%i.o ", n);
	 }*/
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.real.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.real_inv.%i.o ", n);
	}
#ifdef USE_DCT
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.dct1.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.dct2.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.dct3.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.dct4.%i.o ", n);
	}
#endif
	fprintf(fp, "\n%%.o: %%.cpp $(DEPS)\n");
	fprintf(fp, "\t$(CC) -c -o $@ $< $(CFLAGS)\n\n");
	fprintf(fp, "ffttest: $(OBJ) main.o\n");
	fprintf(fp, "\t$(CC) -o $@ $^ $(CFLAGS) -lfftw3\n");
	fprintf(fp, "sfftlib.a: $(OBJ)\n");
	fprintf(fp, "\tar -rcs sfftlib.a $(OBJ)\n");
	fclose(fp);
	print_fft_bests();
	return 0;
}
