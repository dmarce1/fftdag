#include "dag.hpp"
#include "fft.hpp"
#include "util.hpp"
#include "pebble.hpp"

#include <time.h>

constexpr int Nmin = 2;
constexpr int Nmax = 64;

int main(int argc, char **argv) {
	int N = 8;
	auto inputs = math_vertex::new_inputs(2 * N);
	auto outputs = fft(inputs, N, 0);
	std::vector<dag_vertex<math_vertex::properties>> vin, vout;
	for (auto i : inputs) {
		vin.push_back(i);
	}
	for (auto o : outputs) {
		vout.push_back(o);
	}
	pebble_game<math_vertex::properties> game(vin, vout);
	int score = game.search();
	fprintf(stderr, "score = %i\n", score);

	abort();

	fprintf( stderr, "------------------------------COMPLEX-----------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(2 * N);
		auto outputs = fft(inputs, N, 0);
		auto cnt = math_vertex::operation_count(outputs);
		inputs.clear();
		fprintf(stderr, "N = %4i | tot = %4i | add = %4i | mul = %4i | neg = %4i\n", N, cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg);
		auto code = math_vertex::execute_all(outputs);
		std::string fname = "fft.complex." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_complex_") + std::to_string(N) + "(double* x) {\n" + code + "}\n\n";
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
	}
	fprintf( stderr, "------------------------------COMPLEX INVERSE-----------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(2 * N);
		auto outputs = fft(inputs, N, FFT_INV);
		auto cnt = math_vertex::operation_count(outputs);
		inputs.clear();
		fprintf(stderr, "N = %4i | tot = %4i | add = %4i | mul = %4i | neg = %4i\n", N, cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg);
		auto code = math_vertex::execute_all(outputs);
		std::string fname = "fft.complex_inv." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_complex_inv_") + std::to_string(N) + "(double* x) {\n" + code + "}\n\n";
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
	}
	fprintf( stderr, "------------------------------REAL--------------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, N, FFT_REAL);
		auto cnt = math_vertex::operation_count(outputs);
		inputs.clear();
		fprintf(stderr, "N = %4i | tot = %4i | add = %4i | mul = %4i | neg = %4i\n", N, cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg);
		auto code = math_vertex::execute_all(outputs);
		std::string fname = "fft.real." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_real_") + std::to_string(N) + "(double* x) {\n" + code + "}\n\n";
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
	}
	fprintf( stderr, "------------------------------REAL INVERSE-------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, N, FFT_REAL | FFT_INV);
		auto cnt = math_vertex::operation_count(outputs);
		inputs.clear();
		fprintf(stderr, "N = %4i | tot = %4i | add = %4i | mul = %4i | neg = %4i\n", N, cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg);
		auto code = math_vertex::execute_all(outputs);
		std::string fname = "fft.real_inv." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_real_inv_") + std::to_string(N) + "(double* x) {\n" + code + "}\n\n";
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
	}
	fprintf( stderr, "------------------------------DCT-II--------------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, 4 * N, FFT_DCT2);
		auto cnt = math_vertex::operation_count(outputs);
		inputs.clear();
		fprintf(stderr, "N = %4i | tot = %4i | add = %4i | mul = %4i | neg = %4i\n", N, cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg);
		auto code = math_vertex::execute_all(outputs);
		std::string fname = "fft.dct2." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_dct2_") + std::to_string(N) + "(double* x) {\n" + code + "}\n\n";
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
	}

	fprintf( stderr, "------------------------------DCT-II INVERSE--------------------------\n");
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(N);
		auto outputs = fft(inputs, 4 * N, FFT_DCT2 | FFT_INV);
		auto cnt = math_vertex::operation_count(outputs);
		inputs.clear();
		fprintf(stderr, "N = %4i | tot = %4i | add = %4i | mul = %4i | neg = %4i\n", N, cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg);
		auto code = math_vertex::execute_all(outputs);
		std::string fname = "fft.dct2_inv." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_dct2_inv_") + std::to_string(N) + "(double* x) {\n" + code + "}\n\n";
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
	}

	system("cp ../../gen_src/main.cpp .\n");

	FILE* fp = fopen("fft.hpp", "wt");
	fprintf(fp, "#pragma once\n\n");
	fprintf(fp, "#define FFT_NMAX %i\n", Nmax);
	fprintf(fp, "#define FFT_NMIN %i\n\n\n", Nmin);
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_complex_%i(double*);\n", N);
	}
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_complex_inv_%i(double*);\n", N);
	}
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_real_%i(double*);\n", N);
	}
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_real_inv_%i(double*);\n", N);
	}
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_dct2_%i(double*);\n", N);
	}
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_dct2_inv_%i(double*);\n", N);
	}
	fprintf(fp, "\n\nvoid fft_complex(double*, int);\n");
	fprintf(fp, "\n\nvoid fft_complex_inv(double*, int);\n");
	fprintf(fp, "\n\nvoid fft_real(double*, int);\n");
	fprintf(fp, "\n\nvoid fft_real_inv(double*, int);\n");
	fprintf(fp, "\n\nvoid fft_dct2(double*, int);\n");
	fprintf(fp, "\n\nvoid fft_dct2_inv(double*, int);\n");
	fprintf(fp, "\n\n");
	fclose(fp);

	fp = fopen("fft.cpp", "wt");
	fprintf(fp, "#include \"fft.hpp\"\n\n\n");
	fprintf(fp, "using fft_type = void(*)(double*);\n\n");
	fprintf(fp, "fft_type fft_complex_pointer[FFT_NMAX + 1] = {");
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
	fprintf(fp, "fft_type fft_complex_inv_pointer[FFT_NMAX + 1] = {");
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
	fprintf(fp, "\n};\n\n");
	fprintf(fp, "fft_type fft_real_pointer[FFT_NMAX + 1] = {");
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
	fprintf(fp, "fft_type fft_real_inv_pointer[FFT_NMAX + 1] = {");
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
	fprintf(fp, "fft_type fft_dct2_inv_pointer[FFT_NMAX + 1] = {");
	for (int N = 0; N <= Nmax; N++) {
		if (N % 8 == 0) {
			fprintf(fp, "\n\t");
		}
		if (N < Nmin) {
			fprintf(fp, "nullptr");
		} else {
			fprintf(fp, "&fft_dct2_inv_%i", N);
		}
		if (N != Nmax) {
			fprintf(fp, ", ");
		}
	}
	fprintf(fp, "\n};\n\n");
	fprintf(fp, "void fft_complex(double* x, int N) {\n");
	fprintf(fp, "\t(*fft_complex_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
	fprintf(fp, "void fft_complex_inv(double* x, int N) {\n");
	fprintf(fp, "\t(*fft_complex_inv_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
	fprintf(fp, "void fft_real(double* x, int N) {\n");
	fprintf(fp, "\t(*fft_real_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
	fprintf(fp, "void fft_real_inv(double* x, int N) {\n");
	fprintf(fp, "\t(*fft_real_inv_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
	fprintf(fp, "void fft_dct2(double* x, int N) {\n");
	fprintf(fp, "\t(*fft_dct2_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");
	fprintf(fp, "void fft_dct2_inv(double* x, int N) {\n");
	fprintf(fp, "\t(*fft_dct2_inv_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");

	fclose(fp);

	fp = fopen("Makefile", "wt");
	fprintf(fp, "CC=g++\n");
	fprintf(fp, "CFLAGS=-I. -Ofast -march=native\n");
//	fprintf(fp, "CFLAGS=-I. -g -D_GLIBCXX_DEBUG\n");
	fprintf(fp, "DEPS = fft.hpp\n");
	fprintf(fp, "OBJ = main.o fft.o ");
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.complex.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.complex_inv.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.real.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.real_inv.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.dct2.%i.o ", n);
	}
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.dct2_inv.%i.o ", n);
	}
	fprintf(fp, "\n%%.o: %%.cpp $(DEPS)\n");
	fprintf(fp, "\t$(CC) -c -o $@ $< $(CFLAGS)\n\n");
	fprintf(fp, "ffttest: $(OBJ)\n");
	fprintf(fp, "\t$(CC) -o $@ $^ $(CFLAGS) -lfftw3\n");
	fclose(fp);
	return 0;
}
