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

#include <unordered_map>
#include <stack>



typedef enum {
	ADD, SUB, NEG, MUL, IN, CON
} operation_t;

bool is_arithmetic(operation_t op);


class math_vertex {
	struct properties {
		operation_t op;
		name_server::name_ptr name;
		std::shared_ptr<name_server> names;
		double value;
		std::string print_code(const std::vector<properties>& edges);
	};
	dag_vertex<properties> v;
	static math_vertex binary_op(operation_t op, const math_vertex& A, const math_vertex& B);
	static math_vertex unary_op(operation_t op, const math_vertex& A);
public:
	math_vertex() = default;
	math_vertex(const math_vertex&v) = default;
	math_vertex(math_vertex&& v) = default;
	math_vertex& operator=(const math_vertex&v) = default;
	math_vertex& operator=(math_vertex&& v) = default;
	math_vertex(double constant);
	math_vertex& operator=(double constant);
	static math_vertex new_input(std::shared_ptr<name_server> names, std::string&& name);
	static std::vector<math_vertex> new_inputs(std::string&& base_name, int cnt);
	friend math_vertex operator+(const math_vertex& A, const math_vertex& B);
	friend math_vertex operator-(const math_vertex& A, const math_vertex& B);
	friend math_vertex operator*(const math_vertex& A, const math_vertex& B);
	friend math_vertex operator-(const math_vertex& A);
	math_vertex& operator+=(const math_vertex& other);
	math_vertex& operator-=(const math_vertex& other);
	math_vertex& operator*=(const math_vertex& other);
	std::string execute(dag_vertex<properties>::executor& exe);
	static std::string execute_all(std::vector<math_vertex>& vertices);
};




#endif /* MATH_HPP_ */
