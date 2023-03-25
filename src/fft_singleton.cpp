#include "fft.hpp"


std::vector<cmplx> fft_singleton(std::vector<cmplx> xin, int N, int opts) {
	assert(is_prime(N));
	const int M = (N - 1) / 2;
	std::vector<cmplx> xout(N);
	std::vector<cmplx> tp(M + 1);
	std::vector<cmplx> tm(M + 1);
	std::vector<math_vertex> up(M + 1);
	std::vector<math_vertex> um(M + 1);
	std::vector<math_vertex> ap(M + 1);
	std::vector<math_vertex> am(M + 1);
	std::vector<math_vertex> bp(M + 1);
	std::vector<math_vertex> bm(M + 1);
	std::vector<std::vector<bool>> pc(M + 1, std::vector<bool>(M + 1, false));
	std::vector<bool> xc(M + 1, false);
	for (int j = 1; j <= M; j++) {
		tp[j] = xin[j] + xin[N - j];
		tm[j] = xin[j] - xin[N - j];
	}
	for (int j = 1; j <= M; j++) {
		tp[j].x.set_goal();
		tm[j].x.set_goal();
		tp[j].y.set_goal();
		tm[j].y.set_goal();
	}
	for (int i = 1; i <= M; i++) {
		am[i] = 0.0;
		bm[i] = 0.0;
		ap[i] = 0.0;
		bp[i] = 0.0;
	}
	xout[0] = cmplx( { 0.0, 0.0 });
	std::vector<bool> first(M + 1, true);
	for (int i0 = 1; i0 <= M; i0++) {
		for (int j0 = 1; j0 <= M; j0++) {
			for (int i1 = i0 + 1; i1 <= M; i1++) {
				for (int j1 = j0 + 1; j1 <= M; j1++) {
					{
						double c00 = cos(2.0 * M_PI * i0 * j0 / N);
						double c10 = cos(2.0 * M_PI * i1 * j0 / N);
						double c01 = cos(2.0 * M_PI * i0 * j1 / N);
						double c11 = cos(2.0 * M_PI * i1 * j1 / N);
						if (close2(c00, c11) && close2(c10, c01)) {
							double c0 = c00;
							double c1 = c01;
							math_vertex a = math_vertex(0.5 * (c0 + c1));
							math_vertex b = math_vertex(0.5 * (c0 - c1));
							cmplx s0 = (tp[j0] + tp[j1]);
							cmplx s1 = (tp[j0] - tp[j1]);
							bool f = first[i0] && first[i1];
							auto as0 = (f ? xin[0] : cmplx( { 0.0, 0.0 })) + a * s0;
							first[i0] = false;
							first[i1] = false;
							ap[i0] += as0.x + b * s1.x;
							ap[i1] += as0.x - b * s1.x;
							bp[i0] += as0.y + b * s1.y;
							bp[i1] += as0.y - b * s1.y;
							if (!xc[j0] && !xc[j1]) {
								xout[0] += s0;
								xc[j0] = true;
								xc[j1] = true;
							}
							pc[i0][j0] = true;
							pc[i0][j1] = true;
							pc[i1][j0] = true;
							pc[i1][j1] = true;
						}
					}
				}
			}
		}
	}
	for (int i = 1; i <= M; i++) {
		for (int j = 1; j <= M; j++) {
			auto co = cos(2.0 * M_PI * j * i / N);
			if (!pc[i][j]) {
				ap[i] += tp[j].x * co;
				bp[i] += tp[j].y * co;
			}
		}
	}
	for (int i = 1; i <= M; i++) {
		if (!xc[i]) {
			xout[0] += tp[i];
		}
	}
	xout[0] += xin[0];
	for (int i = 1; i <= M; i++) {
		if (first[i]) {
			ap[i] += xin[0].x;
			bp[i] += xin[0].y;
		}
	}
	for (int j = 1; j <= M; j++) {
		ap[j].set_goal();
		bp[j].set_goal();
		xout[0].x.set_goal();
		xout[0].y.set_goal();
	}
	for (int i = 1; i <= M; i++) {
		for (int j = 1; j <= M; j++) {
			auto si = sin(2.0 * M_PI * j * i / N);
			am[i] -= tm[j].y * si;
			bm[i] -= tm[j].x * si;
		}
	}
	for (int j = 1; j <= M; j++) {
		am[j].set_goal();
		bm[j].set_goal();
	}
	for (int i = 1; i <= M; i++) {
		xout[i].x = ap[i] - am[i];
		xout[i].y = bp[i] + bm[i];
		xout[N - i].x = ap[i] + am[i];
		xout[N - i].y = bp[i] - bm[i];
	}
	return std::move(xout);
}
