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

#include <map>
#include <unordered_map>
#include <stack>

typedef enum {
	ADD, SUB, NEG, MUL, IN, CON
} operation_t;

bool is_arithmetic(operation_t op);

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
	math_vertex optimize();
	~math_vertex();
	bool is_neg() const;
	math_vertex get_neg() const;
	void swap(math_vertex& v);
	bool is_zero() const;
	bool is_one() const;
	bool is_neg_one() const;
	int get_value_number() const;
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
	static std::map<double, math_vertex> consts;
	static math_vertex binary_op(operation_t op, math_vertex A, math_vertex B);
	static math_vertex unary_op(operation_t op, math_vertex A);
	void set_database(const std::shared_ptr<name_server>& db);
};

#endif /* MATH_HPP_ */
