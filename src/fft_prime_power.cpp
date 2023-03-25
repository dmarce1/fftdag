#include "fft.hpp"



std::vector<cmplx> fft_prime_power(int R, std::vector<cmplx> xin, int N, int opts) {
	assert(is_prime(R) && prime_factorization(N).size() == 1);
	cmplx I( { 0.0, 1.0 });
	if (N == R) {
		assert(false);
	}
	const int N1 = R;
	const int N2 = N / R;
	std::vector<cmplx> xout(N);
	std::vector<std::vector<cmplx>> sub(N1, std::vector<cmplx>(N2));
	int begin = -(N1 / 2);
	int end = begin + N1;
	const auto I0 = [N]( int n ) {
		while( n < 0 ) {
			n += N;
		}
		return n % N;
	};
	const auto I1 = [N1]( int n ) {
		while( n < 0 ) {
			n += N1;
		}
		return n % N1;
	};
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = begin; n1 < end; n1++) {
			sub[I1(n1)][n2] = xin[I0(n2 * N1 + n1)];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		sub[n1] = fft(sub[n1], N2, opts);
	}
	{
		for (int k1 = begin; k1 < end; k1++) {
			xout[I0(k1 * N2)] = sub[I1(0)][0];
		}
		for (int n1 = 1; n1 < end; n1++) {
			for (int k1 = begin; k1 < end; k1++) {
				const auto w = twiddle(n1 * k1, N1);
				const auto t_0 = sub[I1(n1)][0];
				const auto t_1 = sub[I1(-n1)][0];
				xout[I0(k1 * N2)].x += w.x * (t_0.x + t_1.x);
				xout[I0(k1 * N2)].y += w.x * (t_0.y + t_1.y);
				xout[I0(k1 * N2)].x += w.y * (t_1.y - t_0.y);
				xout[I0(k1 * N2)].y += w.y * (t_0.x - t_1.x);
			}
		}
	}
	for (int k2 = 1; k2 < N2; k2++) {
		for (int k1 = begin; k1 < end; k1++) {
			xout[I0(k2 + k1 * N2)] = sub[I1(0)][k2];
		}
		for (int k1 = begin + 1; k1 < end; k1++) {
			cmplx x = cmplx( { 0.0, 0.0 });
			for (int n1 = 1; n1 < end; n1++) {
				const auto w = twiddle(n1 * (N2 * k1 + k2), N);
				const auto t_0 = sub[I1(n1)][k2];
				const auto t_1 = sub[I1(-n1)][k2];
				x.x = x.x + w.x * (t_0.x + t_1.x) + w.y * (t_1.y - t_0.y);
				;
				x.y = x.y + w.x * (t_0.y + t_1.y) + w.y * (t_0.x - t_1.x);
			}
			xout[I0(k2 + k1 * N2)].x += x.x;
			xout[I0(k2 + k1 * N2)].y += x.y;
			xout[I0(k2 + begin * N2)].x -= x.x;
			xout[I0(k2 + begin * N2)].y -= x.y;
		}
	}
	return xout;
}
