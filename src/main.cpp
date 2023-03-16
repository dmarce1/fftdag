#include "dag.hpp"
#include "fft.hpp"

#include <time.h>
#include "util.hpp"

constexpr int Nmin = 2;
constexpr int Nmax = 64;

int main(int argc, char **argv) {
	for (int N = Nmin; N <= Nmax; N++) {
		auto inputs = math_vertex::new_inputs(2 * N);
		auto outputs = fft(inputs, N);
		auto cnt = math_vertex::operation_count(outputs);
		inputs.clear();
		fprintf(stderr, "N = %4i | tot = %4i | add = %4i | mul = %4i | neg = %4i\n", N, cnt.add + cnt.mul + cnt.neg, cnt.add, cnt.mul, cnt.neg);
		auto code = math_vertex::execute_all(outputs);
		std::string fname = "fft." + std::to_string(N) + ".cpp";
		FILE* fp = fopen(fname.c_str(), "wt");
		code = std::string("\n\nvoid fft_") + std::to_string(N) + "(double* x) {\n" + code + "}\n\n";
		fprintf(fp, "%s\n", code.c_str());
		fclose(fp);
	}

	system("cp ../../gen_src/main.cpp .\n");

	FILE* fp = fopen("fft.hpp", "wt");
	fprintf(fp, "#pragma once\n\n");
	fprintf(fp, "#define FFT_NMAX %i\n", Nmax);
	fprintf(fp, "#define FFT_NMIN %i\n\n\n", Nmin);
	for (int N = Nmin; N <= Nmax; N++) {
		fprintf(fp, "void fft_%i(double*);\n", N);
	}
	fprintf(fp, "\n\nvoid fft(double*, int);\n");
	fprintf(fp, "\n\n");
	fclose(fp);

	fp = fopen("fft.cpp", "wt");
	fprintf(fp, "#include \"fft.hpp\"\n\n\n");
	fprintf(fp, "using fft_type = void(*)(double*);\n\n");
	fprintf(fp, "fft_type fft_pointer[FFT_NMAX + 1] = {");
	for (int N = 0; N <= Nmax; N++) {
		if (N % 8 == 0) {
			fprintf(fp, "\n\t");
		}
		if (N < Nmin) {
			fprintf(fp, "nullptr");
		} else {
			fprintf(fp, "&fft_%i", N);
		}
		if (N != Nmax) {
			fprintf(fp, ", ");
		}
	}
	fprintf(fp, "\n};\n\n");
	fprintf(fp, "void fft(double* x, int N) {\n");
	fprintf(fp, "\t(*fft_pointer[N])(x);\n");
	fprintf(fp, "}\n\n\n");

	fclose(fp);

	fp = fopen("Makefile", "wt");
	fprintf(fp, "CC=g++\n");
	fprintf(fp, "CFLAGS=-I. -Ofast -march=native\n");
	fprintf(fp, "DEPS = fft.hpp\n");
	fprintf(fp, "OBJ = main.o fft.o ");
	for (int n = Nmin; n <= Nmax; n++) {
		fprintf(fp, "fft.%i.o ", n);
	}
	fprintf(fp, "\n%%.o: %%.cpp $(DEPS)\n");
	fprintf(fp, "\t$(CC) -c -o $@ $< $(CFLAGS)\n\n");
	fprintf(fp, "ffttest: $(OBJ)\n");
	fprintf(fp, "\t$(CC) -o $@ $^ $(CFLAGS) -lfftw3\n");
	fclose(fp);
	return 0;
}
