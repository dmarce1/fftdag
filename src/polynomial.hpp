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
	bool zero() const {
		for (auto i = a.begin(); i != a.end(); i++) {
			if (!close2(std::abs(i->second), 0.0)) {
				return false;
			}
		}
		return true;
	}
	T operator[](int i) const {
		auto j = a.find(i);
		assert(j != a.end());
		return j->second;
	}
	T& operator[](int i) {
		return a[i];
	}
	int degree() {
		while (a.size() && close2(std::abs(a.rbegin()->second), 0.0)) {
			a.erase(a.rbegin()->first);
		}
		if (a.size()) {
			return a.rbegin()->first;
		} else {
			return -1;
		}
	}
	std::string to_str() {
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
		for (auto i = A.a.begin(); i != A.a.end(); i++) {
			auto j = B.a.find(i->first);
			if (j != B.a.end()) {
				i->second += j->second;
			}
		}
		for (auto i = B.a.begin(); i != B.a.end(); i++) {
			auto j = A.a.find(i->first);
			if (j == A.a.end()) {
				A[i->first] = i->second;
			}
		}
		return A;
	}
	polynomial operator-(polynomial B) const {
		auto A = *this;
		for (auto i = A.a.begin(); i != A.a.end(); i++) {
			auto j = B.a.find(i->first);
			if (j != B.a.end()) {
				i->second -= j->second;
			}
		}
		for (auto i = B.a.begin(); i != B.a.end(); i++) {
			auto j = A.a.find(i->first);
			if (j == A.a.end()) {
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
	polynomial& operator*=(const polynomial & B) {
		*this = *this * B;
		return *this;
	}
	polynomial& operator/=(const polynomial & B) {
		*this = *this / B;
		return *this;
	}

	template<class U>
	friend class polynomial;

	template<class U>
	polynomial operator*(const polynomial<U>& B) const {
		auto A = *this;
		polynomial C;
		for (auto i = A.a.begin(); i != A.a.end(); i++) {
			for (auto j = B.a.begin(); j != B.a.end(); j++) {
				int nm = i->first + j->first;
				if (C.a.find(nm) == C.a.end()) {
					C[nm] = T(0);
				}
				C[nm] += i->second * j->second;
			}
		}
		return C;
	}
	T leading() {
		while (close2(std::abs(a.rbegin()->second), 0.0)) {
			a.erase(a.rbegin()->first);
		}
		return a.rbegin()->second;
	}
	template<class U>
	polynomial operator%(polynomial<U> b) const {
		auto a = *this;
		polynomial q;
		polynomial r = a;
		int d = b.degree();
		U c = b.leading();
		while (r.degree() >= d) {
			polynomial s;
			s[r.degree() - d] = r.leading() * (U(1) / c);
			q += s;
			r = r - s * b;
			r.a.erase(r.a.rbegin()->first);
		}
		return r;
	}
	template<class U>
	polynomial operator/(polynomial<U> b) const {
		auto a = *this;
		polynomial q;
		polynomial r = a;
		int d = b.degree();
		U c = b.leading();
		if (r.degree() <= d) {
			printf("? %e\n", b.leading());
		} else {
			printf("!\n");
		}
		while (r.degree() >= d) {
			polynomial s;
			s[r.degree() - d] = r.leading() * (U(1) / c);
			q += s;
			r = r - s * b;
			r.a.erase(r.a.rbegin()->first);
		}
		return q;
	}
}
;

template<class T>
polynomial<T> operator*(const T& b, const polynomial<T>& A) {
	return A * b;
}

