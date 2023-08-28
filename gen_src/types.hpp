/*
 * types.hpp
 *
 *  Created on: Apr 6, 2023
 *      Author: dmarce1
 */

#ifndef TYPES333_HPP_
#define TYPES333_HPP_

#include <immintrin.h>

typedef double fft_v2df __attribute__ ((vector_size (16)));
typedef double fft_v4df __attribute__ ((vector_size (32)));

class fft_simd2 {
	fft_v2df a;
public:
	fft_simd2& operator=(const fft_simd2& other) = default;
	inline fft_simd2() = default;
	inline fft_simd2& operator=(double b) {
		a = b - (fft_v2df ) { 0.0, 0.0 };
		return *this;
	}
	inline fft_simd2(double b) {
		a = b - (fft_v2df ) { 0.0, 0.0 };
	}
	inline fft_simd2 operator*(const fft_simd2& b) const {
		fft_simd2 c;
		c.a = a * b.a;
		return c;
	}
	inline fft_simd2 operator+(const fft_simd2& b) const {
		fft_simd2 c;
		c.a = a + b.a;
		return c;
	}
	inline fft_simd2 operator-(const fft_simd2& b) const {
		fft_simd2 c;
		c.a = a - b.a;
		return c;
	}
	inline fft_simd2 operator-() const {
		fft_simd2 c;
		c.a = (fft_v2df ) { 0.0, 0.0 } - a;
		return c;
	}
	inline fft_simd2& operator*=(const fft_simd2& b) {
		*this = *this * b;
		return *this;
	}
	inline fft_simd2& operator+=(const fft_simd2& b) {
		*this = *this + b;
		return *this;
	}
	inline fft_simd2& operator-=(const fft_simd2& b) {
		*this = *this - b;
		return *this;
	}
	inline double operator[](int i) const {
		return a[i];
	}
	inline double& operator[](int i) {
		return a[i];
	}
};

#include <limits>

class fft_simd4 {
	fft_v4df a;
public:
	fft_simd4& operator=(const fft_simd4& other) = default;
	inline fft_simd4() {
	}
	inline fft_simd4& operator=(double b) {
		a = b - (fft_v4df ) { 0.0, 0.0, 0.0, 0.0 };
		return *this;
	}
	inline fft_simd4(double b) {
		*this = b;
	}
	inline fft_simd4 operator*(const fft_simd4& b) const {
		fft_simd4 c;
		c.a = a * b.a;
		return c;
	}
	inline fft_simd4 operator+(const fft_simd4& b) const {
		fft_simd4 c;
		c.a = a + b.a;
		return c;
	}
	inline fft_simd4 operator-(const fft_simd4& b) const {
		fft_simd4 c;
		c.a = a - b.a;
		return c;
	}
	inline fft_simd4 operator-() const {
		fft_simd4 c;
		c = fft_simd4(0.0) - *this;
		return c;
	}
	inline fft_simd4& operator*=(const fft_simd4& b) {
		*this = *this * b;
		return *this;
	}
	inline fft_simd4& operator+=(const fft_simd4& b) {
		*this = *this + b;
		return *this;
	}
	inline void gather(double* base, long long* indices) {
		const __m256i i = *((__m256i *) indices);
		(__m256d &) a = _mm256_i64gather_pd(base, i, 1);
	}
	inline fft_simd4& operator-=(const fft_simd4& b) {
		*this = *this - b;
		return *this;
	}
	inline double operator[](int i) const {
		return a[i];
	}
	inline double& operator[](int i) {
		return a[i];
	}
	fft_simd4 reverse() const {
		constexpr int flags = (3 << 0) | (2 << 2) | (1 << 4) | (0 << 6);
		auto b = _mm256_permute4x64_pd(a, flags);
		fft_simd4 B;
		B.a = b;
		return B;
	}
	fft_simd4 gather(const double* base, const __m256i& indices) {
		a = _mm256_i64gather_pd(base, indices, sizeof(double));
		return *this;
	}
	fft_simd4 load(const double* ptr) {
		a = _mm256_load_pd(ptr);
		return *this;
	}
	void store(double* ptr) {
		_mm256_store_pd(ptr, a);
	}
};

#include <array>

template<int W>
class simd {
	static constexpr int nwidth = 4;
	static constexpr int N = W / nwidth;
	typedef double v4d __attribute__ ((vector_size (32)));
	std::array<v4d, N> q;
public:
	simd& operator=(const simd& other) = default;
	inline simd() {
	}
	inline simd& operator=(double b) {
		for (int n = 0; n < N; n++) {
			q[n] = b - (v4d ) { 0.0, 0.0, 0.0, 0.0 };
		}
		return *this;
	}
	inline simd(double b) {
		*this = b;
	}
	inline simd operator*(const simd& b) const {
		simd c;
		for (int n = 0; n < N; n++) {
			c.q[n] = q[n] * b.q[n];
		}
		return c;
	}
	inline simd operator+(const simd& b) const {
		simd c;
		for (int n = 0; n < N; n++) {
			c.q[n] = q[n] + b.q[n];
		}
		return c;
	}
	inline simd operator-(const simd& b) const {
		simd c;
		for (int n = 0; n < N; n++) {
			c.q[n] = q[n] - b.q[n];
		}
		return c;
	}
	inline simd operator-() const {
		simd c;
		c = simd(0.0) - *this;
		return c;
	}
	inline simd& operator*=(const simd& b) {
		*this = *this * b;
		return *this;
	}
	inline simd& operator+=(const simd& b) {
		*this = *this + b;
		return *this;
	}
	inline simd& operator-=(const simd& b) {
		*this = *this - b;
		return *this;
	}
	inline double operator[](int i) const {
		return q[i / N][i % nwidth];
	}
	inline double& operator[](int i) {
		return q[i / N][i % nwidth];
	}
	simd load(const double* ptr) {
		for (int n = 0; n < N; n++) {
			q[n] = _mm256_loadu_pd(ptr + nwidth * n);
		}
		return *this;
	}
	simd reverse() {
		constexpr int flags = (3 << 0) | (2 << 2) | (1 << 4) | (0 << 6);
		for (int n = 0; n < N - n - 1; n++) {
			const auto tmp = q[n];
			q[n] = q[N - n - 1];
			q[N - n - 1] = tmp;
		}
		for (int n = 0; n < N; n++) {
			q[n] = _mm256_permute4x64_pd(q[n], flags);
		}
		return *this;
	}
};

template<class T>
class complex {
	T x;
	T y;
public:
	inline T real() const {
		return x;
	}
	inline T& real() {
		return x;
	}
	inline T imag() const {
		return y;
	}
	inline T& imag() {
		return y;
	}
	complex() {
	}
	inline complex(T a) {
		x = a;
		y = 0.0;
	}
	inline complex conj() const {
		complex C;
		C.x = x;
		C.y = -y;
		return C;
	}
	inline friend complex Ix(complex A) {
		complex C;
		C.x = -A.y;
		C.y = A.x;
		return C;
	}
	complex(const complex&) = default;
	template<class U>
	friend inline complex<T> operator*(const complex<T>& A, const complex<U>& B) {
		complex<T> C;
		C.real() = A.real() * B.real() - A.imag() * B.imag();
		C.imag() = A.real() * B.imag() + A.imag() * B.real();
		return C;
	}
	friend inline complex<T> operator+(const complex<T>& A, const complex<T>& B) {
		complex<T> C;
		C.x = A.x + B.x;
		C.y = A.y + B.y;
		return C;
	}
	friend inline complex<T> operator-(const complex<T>& A, const complex<T>& B) {
		complex<T> C;
		C.x = A.x - B.x;
		C.y = A.y - B.y;
		return C;
	}
	friend inline complex<T> operator-(const complex<T>& A) {
		complex<T> C;
		C.x = -A.x;
		C.y = -A.y;
		return C;
	}
	friend inline complex<T> operator*(const complex<T>& A, T& B) {
		complex<T> C;
		C.x = A.x * B;
		C.y = A.y * B;
		return C;
	}
	friend inline complex<T> operator*(T A, const complex<T>& B) {
		return B * A;
	}
	friend inline T abs(const complex<T>& A) {
		return sqrt(A.x * A.x + A.y * A.y);
	}
	template<class U>
	inline complex<T>& operator*=(const complex<U>& other) {
		*this = *this * other;
		return *this;
	}
	inline complex<T>& operator+=(const complex<T>& other) {
		*this = *this + other;
		return *this;
	}
	inline complex<T>& operator-=(const complex<T>& other) {
		*this = *this - other;
		return *this;
	}
	inline complex<T>& operator*=(T other) {
		*this = *this * other;
		return *this;
	}
};

inline complex<double> operator*(const complex<double>& A, const complex<double>& B) {
	complex<double> C;
	const __m128d a = (const __m128d &) A;
	const __m128d b = (const __m128d &) B;
	__m128d& c = (__m128d&) C;
	const __m128d r = _mm_mul_pd(a, b);
	const __m128d i = _mm_mul_pd(a, _mm_permute_pd(b, 0x1));
	const __m128d u = _mm_shuffle_pd(r, i, 0x0);
	const __m128d v = _mm_shuffle_pd(r, i, 0x3);
	c = _mm_addsub_pd(u, v);
	return C;
}

inline __m256d mul(const __m256d& a, const __m256d& b) {
	__m256d c;
	const __m256d r = _mm256_mul_pd(a, b);
	const __m256d i = _mm256_mul_pd(a, _mm256_permute_pd(b, 0x5));
	const __m256d u = _mm256_shuffle_pd(r, i, 0x0);
	const __m256d v = _mm256_shuffle_pd(r, i, 0xf);
	c = _mm256_addsub_pd(u, v);
	return c;
}

#endif /* TYPES_HPP_ */
