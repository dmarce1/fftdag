#include "fft.hpp"

template<class T>
class polynomial {
	std::vector<T> a;
public:
	T operator[](int i) const {
		return a[i];
	}
	T& operator[](int i) {
		return a[i];
	}
	int degree() const {
		return a.size() - 1;
	}
	std::string to_str() const {
		std::string str = "(";
		for (int n = degree(); n >= 0; n--) {
			str += std::to_string(a[n]) + ",";
		}
		str += ")";
		return str;
	}
	polynomial operator+(polynomial B) const {
		auto C = *this;
		int deg = std::max(degree(), B.degree());
		B.a.resize(deg + 1, T(0));
		C.a.resize(deg + 1, T(0));
		for (int n = 0; n <= deg; n++) {
			C[n] += B[n];
		}
		return C;
	}
	polynomial operator-(polynomial B) const {
		auto C = *this;
		int deg = std::max(degree(), B.degree());
		B.resize(deg + 1, T(0));
		C.resize(deg + 1, T(0));
		for (int n = 0; n <= deg; n++) {
			C[n] -= B[n];
		}
		return C;
	}
	polynomial& operator+=(const polynomial& B) {
		*this = *this + B;
		return *this;
	}
	polynomial& operator-=(const polynomial& B) {
		*this = *this - B;
		return *this;
	}
	polynomial operator*(const T& b) const {
		auto A = *this;
		for (int n = 0; n <= A.degree(); n++) {
			A[n] *= b;
		}
	}
	polynomial& operator*=(const T& b) {
		*this = *this * b;
		return *this;
	}
	void resize(int deg) {
		a.resize(deg + 1, T(0));
	}
	polynomial operator*(const polynomial& B) const {
		auto A = *this;
		polynomial C;
		int deg = A.degree() + B.degree();
		C.resize(deg);
		for (int n = 0; n <= A.degree(); n++) {
			for (int m = 0; m <= B.degree(); m++) {
				C[n + m] += A[n] * B[m];
			}
		}
		return C;
	}
	polynomial operator%(polynomial D) {
		auto Q = *this;
		int deg = Q.degree();
		for (int d = deg; d >= D.degree(); d--) {
			const T a = Q[d] / D[D.degree()];
			for (int n = 0; D.degree() - n >= 0 ; n++) {
				Q[d - n] -= a * D[D.degree() - n];
			}
		}
		return Q;
	}
};

template<class T>
polynomial<T> operator*(const T& b, const polynomial<T>& A) {
	return A * b;
}

