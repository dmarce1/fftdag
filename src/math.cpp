#include "math.hpp"

bool is_arithmetic(operation_t op) {
	switch (op) {
	case ADD:
	case SUB:
	case NEG:
	case MUL:
		return true;
	default:
		return false;
	}
}

std::string math_vertex::properties::print_code(const std::vector<properties>& edges) {
	std::string code;
	if (name == nullptr) {
		assert(names != nullptr);
		name = names->generate_name();
	}
	if (is_arithmetic(op)) {
		char* ptr;
		switch (op) {
		case ADD:
			asprintf(&ptr, "\t%s = %s + %s;\n", name->c_str(), edges[0].name->c_str(), edges[1].name->c_str());
			break;
		case SUB:
			asprintf(&ptr, "\t%s = %s - %s;\n", name->c_str(), edges[0].name->c_str(), edges[1].name->c_str());
			break;
		case MUL:
			asprintf(&ptr, "\t%s = %s * %s;\n", name->c_str(), edges[0].name->c_str(), edges[1].name->c_str());
			break;
		case NEG:
			asprintf(&ptr, "\t%s = -%s;\n", name->c_str(), edges[0].name->c_str());
			break;
		default:
			break;
		}
		fprintf( stderr, "%s", ptr);
		code = ptr;
		free(ptr);
	}
	return code;
}

math_vertex math_vertex::binary_op(operation_t op, const math_vertex& A, const math_vertex& B) {
	properties props;
	props.op = op;
	math_vertex C;
	C.v = dag_vertex<properties>::new_(std::move(props));
	C.v.add_edge_in(A.v);
	C.v.add_edge_in(B.v);
	if( A.v.properties().names != nullptr) {
		C.v.properties().names = A.v.properties().names;
	} else if( B.v.properties().names != nullptr) {
		C.v.properties().names = B.v.properties().names;
	} else {
		C.v.properties().names = nullptr;
	}
	return C;
}

math_vertex math_vertex::unary_op(operation_t op, const math_vertex& A) {
	properties props;
	props.op = op;
	math_vertex C;
	C.v = dag_vertex<properties>::new_(std::move(props));
	C.v.add_edge_in(A.v);
	if( A.v.properties().names != nullptr) {
		C.v.properties().names = A.v.properties().names;
	} else {
		C.v.properties().names = nullptr;
	}
	return C;
}

math_vertex::math_vertex(double constant) {
	properties props;
	props.op = CON;
	props.value = constant;
	props.name = std::make_shared<std::string>(std::to_string(constant));
	props.names = nullptr;
	v = dag_vertex<properties>::new_(std::move(props));
}

math_vertex& math_vertex::operator=(double constant) {
	*this = math_vertex(constant);
	return *this;
}

math_vertex math_vertex::new_input(std::shared_ptr<name_server> names, std::string&& name) {
	properties props;
	math_vertex C;
	props.op = IN;
	auto tmp = names->reserve_name(std::move(name));
	props.name = tmp.first;
	assert(tmp.second == "");
	props.names = names;
	C.v = dag_vertex<properties>::new_(std::move(props));
	return std::move(C);
}

std::vector<math_vertex> math_vertex::new_inputs(std::string&& base_name, int cnt) {
	std::vector<math_vertex> inputs;
	auto names = std::make_shared<name_server>();
	for (int i = 0; i < cnt; i++) {
		inputs.push_back(new_input(names, base_name + "[" + std::to_string(i) + "]"));
	}
	return std::move(inputs);
}

math_vertex operator+(const math_vertex& A, const math_vertex& B) {
	return math_vertex::binary_op(ADD, A, B);
}

math_vertex operator-(const math_vertex& A, const math_vertex& B) {
	return math_vertex::binary_op(SUB, A, B);
}

math_vertex operator*(const math_vertex& A, const math_vertex& B) {
	return math_vertex::binary_op(MUL, A, B);
}

math_vertex operator-(const math_vertex& A) {
	return math_vertex::unary_op(NEG, A);
}

math_vertex& math_vertex::operator+=(const math_vertex& other) {
	*this = *this + other;
	return *this;
}

math_vertex& math_vertex::operator-=(const math_vertex& other) {
	*this = *this - other;
	return *this;
}

math_vertex& math_vertex::operator*=(const math_vertex& other) {
	*this = *this * other;
	return *this;
}

std::string math_vertex::execute(dag_vertex<properties>::executor& exe) {
	std::string code;
	v.execute(exe, [&code](properties& in, const std::vector<properties>& edges) {
		code += in.print_code(edges);
	});
	return code;
}

std::string math_vertex::execute_all(std::vector<math_vertex>& vertices) {
	std::string code;
	dag_vertex<properties>::executor exe;
	for (auto& v : vertices) {
		code += v.execute(exe);
	}
	return std::move(code);
}
