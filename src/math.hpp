/*
 * math.hpp
 *
 *  Created on: Mar 8, 2023
 *      Author: dmarce1
 */

#ifndef MATH_HPP_
#define MATH_HPP_

#include "dag.hpp"
#include "names.hpp"
#include "util.hpp"
#include "assoc_set.hpp"

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <stack>

typedef enum {
	INVALID, SUB, NEG, MUL, IN, CON, ADD, FMA, NFMA, FMS, NFMS
} operation_t;

bool is_arithmetic(operation_t op);

class cse_database;
class cmplx;
class math_vertex;

class math_vertex {
public:
	class value_number {
		operation_t op;
		int sgn;
		assoc_set x;
		friend class math_vertex;
	public:
		bool operator==(const value_number& other) const;
		value_number operator-() const;
		value_number();
		std::string to_string() const;
	};
	struct properties {
		std::shared_ptr<name_server> names;
		operation_t op;
		name_server::name_ptr name;
		double value;
		int out_num;
		bool cse;
		int depth;
		int group_id;
		bool goal;
		std::shared_ptr<value_number> vnum;
		std::string print_code(const std::vector<math_vertex>& edges);
		properties();
	};
	class weak_ref {
		dag_vertex<properties>::weak_ref ptr;
	public:
		weak_ref() = default;
		weak_ref(const math_vertex& v);
		bool operator<(const weak_ref&) const;
		bool operator==(const weak_ref&) const;
		int use_count() const {
			return ptr.use_count();
		}
		friend class math_vertex;
	};
	struct value_key {
		size_t operator()(const value_number& value) const;
	};
	struct op_cnt_t {
		int add;
		int mul;
		int neg;
		int total() const {
			return add + mul + neg;
		}
	};
	struct distrib_t {
		double c;
		dag_vertex<properties> v;
	};
	struct key {
		size_t operator()(const math_vertex&) const;
	};
	operator dag_vertex<properties>() const;
	void set_group_id(int id);
	int get_group_id() const;
	math_vertex(const weak_ref& ref);
	bool operator<(const math_vertex& other) const;
	bool operator==(const math_vertex& other) const;
	math_vertex optimize();
	int get_unique_id() const {
		return v.get_unique_id();
	}
	~math_vertex();
	void set_goal() {
		v.properties().goal = true;
	}
	static int cse_size() {
		return cse.size();
	}
	bool is_neg() const;
	math_vertex get_neg() const;
	bool is_zero() const;
	bool is_one() const;
	bool is_neg_one() const;
	bool is_constant() const;
	math_vertex();
	bool valid() const;
	math_vertex(const math_vertex&v) = default;
	math_vertex(math_vertex&& v) = default;
	math_vertex(const dag_vertex<properties>& v);
	math_vertex(dag_vertex<properties> && v);
	math_vertex(double constant);
	math_vertex optimize_fma();
	math_vertex& operator=(const math_vertex& other) {
		if (get_unique_id() != other.get_unique_id()) {
			v = other.v;
		}
		return *this;
	}
	math_vertex& operator=(math_vertex&& other) {
		if (get_unique_id() != other.get_unique_id()) {
			v = std::move(other.v);
		}
		return *this;
	}
	math_vertex& operator=(double constant);
//	math_vertex get_edge_out(int i) const;
//	int get_edge_out_count() const;
	math_vertex get_edge_in(int i) const;
	int get_edge_in_count() const;
	void set_value_number(value_number&&);
	assoc_set associative_adds() const;
	std::pair<assoc_set, int> associative_muls() const;
	void replace_edge(const math_vertex&, math_vertex&&);
	math_vertex& operator+=(const math_vertex& other);
	math_vertex& operator-=(const math_vertex& other);
	math_vertex& operator*=(const math_vertex& other);
	operation_t get_op() const;
	double get_value() const;
	op_cnt_t operation_count(dag_vertex<properties>::executor&);
	static std::vector<math_vertex> available2execute(dag_vertex<properties>::executor& exe, dag_vertex<properties>::executor& exe1, std::vector<math_vertex>&);
	std::vector<math_vertex> available2execute(dag_vertex<properties>::executor& exe, dag_vertex<properties>::executor& exe1);
	std::string execute(dag_vertex<properties>::executor& exe);
	static math_vertex new_input(std::shared_ptr<name_server> db);
	static op_cnt_t operation_count(std::vector<math_vertex>);
	static op_cnt_t operation_count(std::vector<cmplx>);
	static std::vector<math_vertex> new_inputs(int cnt);
	static std::pair<std::string, int> execute_all(std::vector<math_vertex>&&, std::vector<math_vertex>& vertices, bool cmpx, int simdsz);
	static void optimize(std::vector<math_vertex>& vertices);
	friend math_vertex operator+(const math_vertex& A, const math_vertex& B);
	friend math_vertex operator-(const math_vertex& A, const math_vertex& B);
	friend math_vertex operator*(const math_vertex& A, const math_vertex& B);
	friend math_vertex operator-(const math_vertex& A);
	struct cse_entry {
		math_vertex::weak_ref ptr;
		bool vacate;
		cse_entry();
		cse_entry& operator=(const math_vertex&);
	};
	static void print_cse();
	static void reset();
private:
	dag_vertex<properties> v;
	bool check_cse();
	static bool first_init;
	static std::map<double, math_vertex> consts;
	static std::unordered_map<value_number, cse_entry, value_key> cse;
	static math_vertex tri_op(operation_t op, math_vertex A, math_vertex B, math_vertex C);
	static math_vertex binary_op(operation_t op, math_vertex A, math_vertex B);
	static math_vertex unary_op(operation_t op, math_vertex A);
	static std::vector<math_vertex> essential_constants;
	void set_database(const std::shared_ptr<name_server>& db);
};

struct cmplx {
	math_vertex x;
	math_vertex y;
	cmplx() = default;
	cmplx(const cmplx &) = default;
	void set_goal() {
		x.set_goal();
		y.set_goal();
	}
	cmplx operator+(cmplx other) const {
		cmplx C;
		C.x = x + other.x;
		C.y = y + other.y;
		return C;
	}
	cmplx operator-(cmplx other) const {
		cmplx C;
		C.x = x - other.x;
		C.y = y - other.y;
		return C;
	}
	cmplx operator-() const {
		cmplx C;
		C.x = -x;
		C.y = -y;
		return C;
	}
	cmplx operator*(math_vertex other) const {
		cmplx C;
		C.x = x * other;
		C.y = y * other;
		return C;
	}
	cmplx conj() const {
		cmplx C;
		C.x = x;
		C.y = -y;
		return C;
	}
	cmplx operator*(cmplx other) const {
		cmplx C;
		C.x = x * other.x - y * other.y;
		C.y = x * other.y + y * other.x;
		return C;
	}
	cmplx operator+=(cmplx other) {
		*this = *this + other;
		return *this;
	}
	cmplx operator-=(cmplx other) {
		*this = *this - other;
		return *this;
	}
	cmplx operator*=(cmplx other) {
		*this = *this * other;
		return *this;
	}
	cmplx operator*=(math_vertex other) {
		*this = *this * other;
		return *this;
	}
	cmplx& operator=(const cmplx&) = default;
	cmplx& operator=(const std::complex<double>& a) {
		x = a.real();
		y = a.imag();
		return *this;
	}
	cmplx(std::complex<double> a) {
		x = a.real();
		y = a.imag();
	}
};

inline cmplx operator*(math_vertex a, cmplx b) {
	return b * a;
}

inline cmplx twiddle(int nk, int N) {
	cmplx C;
	const double theta = -2.0 * M_PI * nk / N;
	C.x = cos(theta);
	C.y = sin(theta);
	return C;
}

inline bool is_arithmetic(operation_t op) {
	switch (op) {
	case ADD:
	case SUB:
	case NEG:
	case MUL:
	case FMA:
	case FMS:
	case NFMA:
	case NFMS:
		return true;
	default:
		return false;
	}
}

inline bool is_fma(operation_t op) {
	switch (op) {
	case FMA:
	case FMS:
	case NFMA:
	case NFMS:
		return true;
	default:
		return false;
	}
}

namespace std {
inline double abs(cmplx z) {
	if (z.x.get_op() == CON && z.y.get_op() == CON) {
		double x = z.x.get_value();
		double y = z.y.get_value();
		return sqrt(x * x + y * y);
	}
	return 1.0;
}
}

#endif /* MATH_HPP_ */
