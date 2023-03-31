#include "fft.hpp"

#include <map>

bool close2(double, double);

template<class T>
class polynomial {
	std::map<int, T> a;
public:
	bool zero(int i) const {
		return a.find(i) == a.end();
	}
	T operator[](int i) const {
		auto j = a.find(i);
		assert(j != a.end());
		return j->second;
	}
	T& operator[](int i) {
		return a[i];
	}
	int degree() const {
		return a.rbegin()->first;
	}
	std::string to_str() const {
		std::string str = "";
		for (auto i = a.rbegin(); i != a.rend(); i++) {
			int n = i->first;
			auto b = i->second;
			if (!close2(std::abs(b), 0.0)) {
				if (n != degree()) {
					str += " + ";
				}
				str += "(" + std::to_string(b.real());
				if (!close2(b.imag(), 0.0)) {
					str += " + i" + std::to_string(b.imag());
				}
				str += ")";
				if (n > 0) {
					str += " z";
					if (n > 1) {
						str += "^" + std::to_string(n);
					}
				}
			}
		}
		return str;
	}
	polynomial operator+(polynomial B) const {
		auto A = *this;
		for (auto i = A.begin(); i != A.end(); i++) {
			auto j = B.find(i->first);
			if (j != B.end()) {
				i->second += j->second;
			}
		}
		for (auto i = B.begin(); i != B.end(); i++) {
			auto j = A.find(i->first);
			if (j == A.end()) {
				A[i->first] = i->second;
			}
		}
		return A;
	}
	polynomial operator-(polynomial B) const {
		auto A = *this;
		for (auto i = A.begin(); i != A.end(); i++) {
			auto j = B.find(i->first);
			if (j != B.end()) {
				i->second -= j->second;
			}
		}
		for (auto i = B.begin(); i != B.end(); i++) {
			auto j = A.find(i->first);
			if (j == A.end()) {
				A[i->first] = -i->second;
			}
		}
		return A;
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
	polynomial operator*(const polynomial& B) const {
		auto A = *this;
		polynomial C;
		int deg = A.degree() + B.degree();
		for (auto i = A.a.begin(); i != A.a.end(); i++) {
			for (auto j = B.a.begin(); j != B.a.end(); j++) {
				int nm = i->first + j->first;
				if (C.a.find(nm) == C.a.end()) {
					C[nm] = 0.0;
				}
				C[nm] += i->second * j->second;
			}
		}
		return C;
	}
	template<class U>
	polynomial operator%(polynomial<U> D) const {
		auto Q = *this;
		int deg = Q.degree();
		for (int d = deg; d >= D.degree(); d--) {
			if (!Q.zero(d)) {
				const T a = Q[d] * (1.0 / D[D.degree()]);
				Q.a.erase(d);
				for (int n = 1; D.degree() - n >= 0; n++) {
					if (!D.zero(D.degree() - n)) {
						if (Q.zero(d - n)) {
							Q[d - n] = 0.0;
						}
						Q[d - n] -= a * D[D.degree() - n];
					}
				}
			}
		}
		return Q;
	}
	template<class U>
	polynomial operator/(polynomial<U> D) const {
		auto Q = *this;
		int qd = Q.degree();
		int dd = D.degree();
		polynomial P;
		for (int d = qd; d >= D.degree(); d--) {
			if (!Q.zero(d)) {
				P[d - dd] = Q[d] * (1.0 / D[dd]);
				Q.a.erase(d);
				for (int n = 1; dd - n >= 0; n++) {
					if (!D.zero(D.degree() - n)) {
						if (Q.zero(d - n)) {
							Q[d - n] = 0.0;
						}
						Q[d - n] -= P[d - dd] * D[dd - n];
					}
				}
			}
		}
		return P;
	}
}
;

template<class T>
polynomial<T> operator*(const T& b, const polynomial<T>& A) {
	return A * b;
}

