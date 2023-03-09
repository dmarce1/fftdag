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
#include <stack>

typedef enum {
	ADD, SUB, NEG, MUL, IN, CON
} operation_t;

bool is_arithmetic(operation_t op);

class cse_database;

class math_vertex {
public:
	struct properties {
		std::shared_ptr<name_server> names;
		operation_t op;
		name_server::name_ptr name;
		double value;
		int out_num;
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
	class value_number {
		operation_t op;
		assoc_set x;
		friend class math_vertex;
	public:
		bool operator==(const value_number& other) const;
		value_number();
	};
	struct value_key {
		size_t operator()( const value_number& value ) const;
	};
	math_vertex(const weak_ref& ref);
	bool operator<(const math_vertex& other) const;
	bool operator==(const math_vertex& other) const;
	math_vertex optimize();
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
	math_vertex get_edge_in(int i) const;
	int get_edge_count() const;
	void replace_edge(const math_vertex&, math_vertex&&);
	math_vertex& operator+=(const math_vertex& other);
	math_vertex& operator-=(const math_vertex& other);
	math_vertex& operator*=(const math_vertex& other);
	std::string execute(dag_vertex<properties>::executor& exe);
	static math_vertex new_input(std::shared_ptr<name_server> db, std::string&& name);
	static std::vector<math_vertex> new_inputs(int cnt);
	static std::string execute_all(std::vector<math_vertex>& vertices);
	friend math_vertex operator+(const math_vertex& A, const math_vertex& B);
	friend math_vertex operator-(const math_vertex& A, const math_vertex& B);
	friend math_vertex operator*(const math_vertex& A, const math_vertex& B);
	friend math_vertex operator-(const math_vertex& A);
private:
	dag_vertex<properties> v;
	value_number vnum;
	value_number compute_value_number();
	static std::map<double, math_vertex> consts;
	static std::unordered_map<value_number, math_vertex::weak_ref, value_key> cse;
	static math_vertex binary_op(operation_t op, math_vertex A, math_vertex B);
	static math_vertex unary_op(operation_t op, math_vertex A);
	void set_database(const std::shared_ptr<name_server>& db);
};

#endif /* MATH_HPP_ */
