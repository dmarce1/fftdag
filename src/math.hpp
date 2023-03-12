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
	ADD, SUB, NEG, MUL, IN, CON, INVALID
} operation_t;

bool is_arithmetic(operation_t op);

class cse_database;

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
		std::shared_ptr<value_number> vnum;
		std::string print_code(const std::vector<properties>& edges);
		properties();
	};
	class weak_ref {
		dag_vertex<properties>::weak_ref ptr;
	public:
		weak_ref() = default;
		weak_ref(const math_vertex& v);
		bool operator<(const weak_ref&) const;
		bool operator==(const weak_ref&) const;
		friend class math_vertex;
	};
	struct value_key {
		size_t operator()(const value_number& value) const;
	};
	struct op_cnt_t {
		int add;
		int mul;
		int neg;
	};
	struct distrib_t {
		double c;
		dag_vertex<properties> v;
	};
	struct key {
		size_t operator()(const math_vertex&) const;
	};
	math_vertex(const weak_ref& ref);
	bool operator<(const math_vertex& other) const;
	bool operator==(const math_vertex& other) const;
	math_vertex optimize();
	math_vertex distribute_muls();
	math_vertex associate_adds();
	~math_vertex();
	bool is_neg() const;
	math_vertex get_neg() const;
	bool is_zero() const;
	bool is_one() const;
	bool is_neg_one() const;
	math_vertex() = default;
	math_vertex(const math_vertex&v) = default;
	math_vertex(math_vertex&& v) = default;
	math_vertex& operator=(const math_vertex&v) = default;
	math_vertex& operator=(math_vertex&& v) = default;
	math_vertex(const dag_vertex<properties>& v);
	math_vertex(dag_vertex<properties> && v);
	math_vertex(double constant);
	math_vertex& operator=(double constant);
	math_vertex get_edge_out(int i) const;
	int get_edge_out_count() const;
	math_vertex get_edge_in(int i) const;
	int get_edge_in_count() const;
	void set_value_number(value_number&&);
	assoc_set associative_adds() const;
	std::vector<distrib_t> distributive_muls();
	std::pair<assoc_set, int> associative_muls() const;
	void replace_edge(const math_vertex&, math_vertex&&);
	math_vertex& operator+=(const math_vertex& other);
	math_vertex& operator-=(const math_vertex& other);
	math_vertex& operator*=(const math_vertex& other);
	std::unordered_set<math_vertex, key> collect_additive_terms(std::unordered_set<math_vertex, key>&);
	std::unordered_set<math_vertex, key> collect_additive_terms_up(std::unordered_set<math_vertex, key>& path);
	std::unordered_set<math_vertex, key> collect_additive_terms_down(std::unordered_set<math_vertex, key>& path);
	operation_t get_op() const;
	math_vertex post_optimize();
	double get_value() const;
	op_cnt_t operation_count(dag_vertex<properties>::executor&);
	std::string execute(dag_vertex<properties>::executor& exe);
	static math_vertex new_input(std::shared_ptr<name_server> db, std::string&& name);
	static op_cnt_t operation_count(std::vector<math_vertex>&);
	static std::vector<math_vertex> new_inputs(int cnt);
	static std::string execute_all(std::vector<math_vertex>& vertices);
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
private:
	dag_vertex<properties> v;
	bool check_cse();
	static std::map<double, math_vertex> consts;
	static std::unordered_map<value_number, cse_entry, value_key> cse;
	static math_vertex binary_op(operation_t op, math_vertex A, math_vertex B);
	static math_vertex unary_op(operation_t op, math_vertex A);
	void set_database(const std::shared_ptr<name_server>& db);
};

#endif /* MATH_HPP_ */
