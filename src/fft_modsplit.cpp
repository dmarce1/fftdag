#include "fft.hpp"

std::vector<cmplx> fft_modsplitS(std::vector<cmplx> xin, int N);
std::vector<cmplx> fft_modsplitS2(std::vector<cmplx> xin, int N);
std::vector<cmplx> fft_modsplitS4(std::vector<cmplx> xin, int N);

double scale_factor(int N, int k) {
	if (N <= 4) {
		return 1.0;
	} else {
		int k4 = k % (N / 4);
		if (k4 <= N / 8) {
			return cos(2.0 * M_PI * k4 / N) * scale_factor(N / 4, k4);
		} else {
			return sin(2.0 * M_PI * k4 / N) * scale_factor(N / 4, k4);
		}
	}
}


static std::vector<cmplx> split_even(const std::vector<cmplx>& xin, int N) {
	std::vector<cmplx> xout;
	for (int n = 0; n < N; n += 2) {
		xout.push_back(xin[n]);
	}
	return std::move(xout);
}

static std::vector<cmplx> split_podd(const std::vector<cmplx>& xin, int N) {
	std::vector<cmplx> xout;
	for (int n = 1; n < N; n += 4) {
		xout.push_back(xin[n]);
	}
	return std::move(xout);
}

static std::vector<cmplx> split_nodd(const std::vector<cmplx>& xin, int N) {
	std::vector<cmplx> xout;
	for (int n = -1; n < N - 4; n += 4) {
		xout.push_back(xin[(n + N) % N]);
	}
	return std::move(xout);
}

std::vector<cmplx> fft_modsplit(std::vector<cmplx> xin, int N, int opts) {
	cmplx I( { 0.0, 1.0 });
	if (N == 0) {
	} else if (N == 1) {
		return std::move(xin);
	} else if (N == 2) {
		auto tmp = xin[0];
		xin[0] = xin[0] + xin[1];
		xin[1] = tmp - xin[1];
		return std::move(xin);
	}
	assert(N%2==0);
	std::vector<cmplx> xout(N);
	auto u0 = fft_modsplit(split_even(xin, N), N / 2, opts);
	auto zp = fft_modsplitS(split_podd(xin, N), N / 4);
	auto zn = fft_modsplitS(split_nodd(xin, N), N / 4);
	for (int k = 0; k < N / 4; k++) {
		double s = scale_factor(N / 4, k);
		cmplx wp;
		cmplx wn;
		double theta = -2.0 * M_PI * k / N;
		wp.x = s * cos(theta);
		wp.y = s * sin(theta);
		wn = wp.conj();
		xout[k + 0 * N / 4] = u0[k + 0 / 4] + 1 * (wp * zp[k] + wn * zn[k]);
		xout[k + 1 * N / 4] = u0[k + N / 4] - I * (wp * zp[k] - wn * zn[k]);
		xout[k + 2 * N / 4] = u0[k + 0 / 4] - 1 * (wp * zp[k] + wn * zn[k]);
		xout[k + 3 * N / 4] = u0[k + N / 4] + I * (wp * zp[k] - wn * zn[k]);
	}
	for (auto x : xout) {
//		x.x.set_goal();
//		x.y.set_goal();
	}
	return std::move(xout);
}

std::vector<cmplx> fft_modsplitS(std::vector<cmplx> xin, int N) {
	cmplx I( { 0.0, 1.0 });
	if (N == 0) {
	} else if (N == 1) {
		return std::move(xin);
	} else if (N == 2) {
		auto tmp = xin[0];
		xin[0] = xin[0] + xin[1];
		xin[1] = tmp - xin[1];
		return std::move(xin);
	}
	assert(N%2==0);
	std::vector<cmplx> xout(N);
	auto u0 = fft_modsplitS2(split_even(xin, N), N / 2);
	auto zp = fft_modsplitS(split_podd(xin, N), N / 4);
	auto zn = fft_modsplitS(split_nodd(xin, N), N / 4);
	for (int k = 0; k < N / 4; k++) {
		double s = scale_factor(N / 4, k) / scale_factor(N, k);
		cmplx tp;
		cmplx tn;
		double theta = -2.0 * M_PI * k / N;
		tp.x = s * cos(theta);
		tp.y = s * sin(theta);
		tn = tp.conj();
		xout[k + 0 * N / 4] = u0[k + 0 / 4] + 1 * (tp * zp[k] + tn * zn[k]);
		xout[k + 1 * N / 4] = u0[k + N / 4] - I * (tp * zp[k] - tn * zn[k]);
		xout[k + 2 * N / 4] = u0[k + 0 / 4] - 1 * (tp * zp[k] + tn * zn[k]);
		xout[k + 3 * N / 4] = u0[k + N / 4] + I * (tp * zp[k] - tn * zn[k]);
	}
	for (auto x : xout) {
//		x.x.set_goal();
//		x.y.set_goal();
	}
return std::move(xout);
}

std::vector<cmplx> fft_modsplitS2(std::vector<cmplx> xin, int N) {
	cmplx I( { 0.0, 1.0 });
	if (N == 0) {
	} else if (N == 1) {
		return std::move(xin);
	} else if (N == 2) {
		auto tmp = xin[0];
		xin[0] = xin[0] + xin[1];
		xin[1] = tmp - xin[1];
		return std::move(xin);
	}
	assert(N%2==0);
	std::vector<cmplx> xout(N);
	auto u0 = fft_modsplitS4(split_even(xin, N), N / 2);
	auto zp = fft_modsplitS(split_podd(xin, N), N / 4);
	auto zn = fft_modsplitS(split_nodd(xin, N), N / 4);
	for (int k = 0; k < N / 4; k++) {
		double s = scale_factor(N / 4, k) / scale_factor(N, k);
		cmplx tp;
		cmplx tn;
		double theta = -2.0 * M_PI * k / N;
		tp.x = s * cos(theta);
		tp.y = s * sin(theta);
		tn = tp.conj();
		xout[k + 0 * N / 4] = u0[k + 0 / 4] + 1 * (tp * zp[k] + tn * zn[k]) * (scale_factor(N, k) / scale_factor(2 * N, k + 0 / 4));
		xout[k + 1 * N / 4] = u0[k + N / 4] - I * (tp * zp[k] - tn * zn[k]) * (scale_factor(N, k) / scale_factor(2 * N, k + N / 4));
		xout[k + 2 * N / 4] = u0[k + 0 / 4] - 1 * (tp * zp[k] + tn * zn[k]) * (scale_factor(N, k) / scale_factor(2 * N, k + 0 / 4));
		xout[k + 3 * N / 4] = u0[k + N / 4] + I * (tp * zp[k] - tn * zn[k]) * (scale_factor(N, k) / scale_factor(2 * N, k + N / 4));
	}
	for (auto x : xout) {
//		x.x.set_goal();
//		x.y.set_goal();
	}
	return std::move(xout);
}

std::vector<cmplx> fft_modsplitS4(std::vector<cmplx> xin, int N) {
	cmplx I( { 0.0, 1.0 });
	if (N == 0) {
	} else if (N == 1) {
		return std::move(xin);
	} else if (N == 2) {
		auto tmp = xin[0];
		auto s0 = (scale_factor(2, 0) / scale_factor(8, 0));
		auto s1 = (scale_factor(2, 1) / scale_factor(8, 1));
		auto x0 = s0 * (xin[0] + xin[1]);
		auto x1 = s1 * (xin[0] - xin[1]);
		xin[0] = x0;
		xin[1] = x1;
		return std::move(xin);
	}
	assert(N%2==0);
	std::vector<cmplx> xout(N);
	auto u0 = fft_modsplitS2(split_even(xin, N), N / 2);
	auto zp = fft_modsplitS(split_podd(xin, N), N / 4);
	auto zn = fft_modsplitS(split_nodd(xin, N), N / 4);
	for (int k = 0; k < N / 4; k++) {
		double s = scale_factor(N / 4, k) / scale_factor(N, k);
		cmplx tp;
		cmplx tn;
		double theta = -2.0 * M_PI * k / N;
		tp.x = s * cos(theta);
		tp.y = s * sin(theta);
		tn = tp.conj();
		xout[k + 0 * N / 4] = (u0[k + 0 / 4] + 1 * (tp * zp[k] + tn * zn[k])) * (scale_factor(N, k) / scale_factor(4 * N, k + 0 * N / 4));
		xout[k + 1 * N / 4] = (u0[k + N / 4] - I * (tp * zp[k] - tn * zn[k])) * (scale_factor(N, k) / scale_factor(4 * N, k + 1 * N / 4));
		xout[k + 2 * N / 4] = (u0[k + 0 / 4] - 1 * (tp * zp[k] + tn * zn[k])) * (scale_factor(N, k) / scale_factor(4 * N, k + 2 * N / 4));
		xout[k + 3 * N / 4] = (u0[k + N / 4] + I * (tp * zp[k] - tn * zn[k])) * (scale_factor(N, k) / scale_factor(4 * N, k + 3 * N / 4));
	}
	for (auto x : xout) {
//		x.x.set_goal();
//		x.y.set_goal();
	}
	return std::move(xout);
}
