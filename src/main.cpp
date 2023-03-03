#include "dag.hpp"
#include "fft.hpp"

void print_code1(int N) {
	printf("\n"
			"\n"
			"\n"
			"#include <fftw3.h>\n"
			"#include <cmath>\n"
			"#include <vector>\n"
			"#include <complex>\n"
			"#include <unordered_map>\n"
			"\n");
	printf("void test(double* x) {\n");
}

void print_code2(int N) {
	printf("}\n");
	printf("\n"
			"\n"
			"void fftw(std::vector<std::complex<double>>& x) {\n"
			"\tconst int N = x.size();\n"
			"\tstatic std::unordered_map<int, fftw_plan> plans;\n"
			"\tstatic std::unordered_map<int, fftw_complex*> in;\n"
			"\tstatic std::unordered_map<int, fftw_complex*> out;\n"
			"\tif (plans.find(N) == plans.end()) {\n"
			"\t\tin[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);\n"
			"\t\tout[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);\n"
			"\t\tplans[N] = fftw_plan_dft_1d(N, in[N], out[N], FFTW_FORWARD, FFTW_ESTIMATE);\n"
			"\t}\n"
			"\tauto* i = in[N];\n"
			"\tauto* o = out[N];\n"
			"\tfor (int n = 0; n < N; n++) {\n"
			"\t\ti[n][0] = x[n].real();\n"
			"\t\ti[n][1] = x[n].imag();\n"
			"\t}\n"
			"\tfftw_execute(plans[N]);\n"
			"\tfor (int n = 0; n < N; n++) {\n"
			"\t\tx[n].real(o[n][0]);\n"
			"\tx[n].imag(o[n][1]);\n"
			"\t}\n"
			"}\n"
			"\n"
			"\n"
			"double rand1() {\n"
			"\treturn (rand() + 0.5) / RAND_MAX;\n"
			"}\n"
			"\n");
	printf("#include <chrono>\n"
			"class timer {\n"
			"\tstd::chrono::time_point<std::chrono::high_resolution_clock> start_time;\n"
			"\tdouble time;\n"
			"public:\n"
			"\tinline timer() {\n"
			"\t\ttime = 0.0;\n"
			"\t}\n"
			"\tinline void stop() {\n"
			"\t\tstd::chrono::time_point<std::chrono::high_resolution_clock> stop_time = std::chrono::high_resolution_clock::now();\n"
			"\t\tstd::chrono::duration<double> dur = stop_time - start_time;\n"
			"\t\ttime += dur.count();\n"
			"\t}\n"
			"\tinline void start() {\n"
			"\t\tstart_time = std::chrono::high_resolution_clock::now();\n"
			"\t}\n"
			"\tinline void reset() {\n"
			"\t\ttime = 0.0;\n"
			"\t}\n"
			"\tinline double read() {\n"
			"\t\treturn time;\n"
			"\t}\n"
			"};\n"
			"\n"
			"\n"
			"\n"
			"int main() {\n"
			"\tconstexpr int N = %i;\n"
			"\tstd::vector<double> xin(2 * N);\n"
			"\tstd::vector<std::complex<double>> y(N);\n"
			"\ttimer tm1, tm2;\n"
			"\tfor( int i = 0; i < 256; i++) {\n"
			"\t\tfor( int n = 0; n < 2 * N; n++) {\n"
			"\t\t\txin[n] = rand1();\n"
			"\t\t\tif( n %% 2 == 0 ) {\n"
			"\t\t\t\ty[n / 2].real(xin[n]);\n"
			"\t\t\t} else {\n"
			"\t\t\t\ty[n / 2].imag(xin[n]);\n"
			"\t\t\t}\n"
			"\t\t}\n"
			"\t\ttm1.start();\n"
			"\t\ttest(xin.data());\n"
			"\t\tauto xout = xin;\n"
			"\t\ttm1.stop();\n"
			"\t\ttm2.start();\n"
			"\t\tfftw(y);\n"
			"\t\ttm2.stop();\n"
			"\t\tdouble error = 0.0;\n"
			"\t\tfor( int n = 0; n < N; n++) {\n"
			"\t\t\terror += std::pow(xout[2 * n] - y[n].real(), 2);\n"
			"\t\t\terror += std::pow(xout[2 * n + 1] - y[n].imag(), 2);\n"
			"\t\t\t//printf( \"%%i %%e %%e %%e %%e;\\n\", n, xout[2*n], xout[2*n+1], y[n].real(), y[n].imag())\n;\n"
			"\t\t}\n"
			"\t\terror = error / (2.0 * N);\n"
			"\t\tif( i == 255 ) {\n"
			"\t\t\tprintf( \"Error = %%e\\n\", error );\n"
			"\t\t}\n"
			"\t}\n"
			"\tprintf( \"%%e %%e %%e\\n\", tm1.read(), tm2.read(), tm2.read() / tm1.read() );\n"
			"\t\n"
			"}\n"
			"", N);

}

int test() {
	constexpr int N = 3;
	auto inputs = dag_node::create_inputs(2 * N);
	auto outputs = fft_prime_power(3, inputs, N);
	dag_node::set_outputs(outputs);
	dag_node::optimize();
	//dag_node::print_list();
	print_code1(N);
	dag_node::print_code();
	print_code2(N);
	auto opcnt = dag_node::get_operation_count();
	fprintf(stderr, "tot = %i | add = %i | mul = %i | neg = %i\n", opcnt.tot, opcnt.add, opcnt.mul, opcnt.neg);
	return 0;
}
int main(int argc, char **argv) {
	test();
}
