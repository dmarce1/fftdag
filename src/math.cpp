#include "math.hpp"
#include "util.hpp"

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

bool is_additive(operation_t op) {
	switch (op) {
	case ADD:
	case SUB:
		return true;
	default:
		return false;
	}
}

int edge_count(operation_t op) {
	switch (op) {
	case ADD:
	case SUB:
	case MUL:
		return 2;
	case NEG:
		return 1;
	default:
		return 0;
	}
}

math_vertex::weak_ref::weak_ref(const math_vertex& v) {
	ptr = dag_vertex<properties>::weak_ref(v.v);
}

math_vertex::math_vertex(const weak_ref& ref) {
	v = dag_vertex<properties>(ref.ptr);
}

bool math_vertex::weak_ref::operator<(const math_vertex::weak_ref& other) const {
	return math_vertex(*this) < math_vertex(other);
}

bool math_vertex::weak_ref::operator==(const math_vertex::weak_ref& other) const {
	return math_vertex(*this) == math_vertex(other);
}

bool math_vertex::value_number::operator==(const value_number& other) const {
	return (op == other.op) && (x == other.x);
}

math_vertex::value_number::value_number() {
	op = ADD;
}

size_t math_vertex::value_key::operator()(const value_number& value) const {
	int rc = (int) value.op;
	for (auto i = value.x.begin(); i != value.x.end(); i++) {
		rc = rc * 1664525 + 1013904223;
	}
	return rc;
}

std::string math_vertex::value_number::to_string() const {
	std::string rc = std::to_string((int) op) + std::string("|");
	for (auto i = x.begin(); i != x.end(); i++) {
		rc += ((i->second > 0) ? std::string("+") : std::string("-")) + std::to_string(std::abs(i->second)) + "*(" + std::to_string(i->first) + ") ";
	}
	return rc;
}

math_vertex::properties::properties() {
	out_num = -1;
}

std::string math_vertex::properties::print_code(const std::vector<properties>& edges) {
	std::string code;
	if (out_num != -1) {
		std::string nm = std::string("x[") + std::to_string(out_num) + std::string("]");
		auto tmp = names->reserve_name(std::move(nm));
		name = tmp.first;
		code += tmp.second;
	} else if (name == nullptr) {
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
		code += ptr;
		free(ptr);
	}
	return code;
}

void math_vertex::replace_edge(const math_vertex& a, math_vertex&& b) {
	v.replace_edge_in(a.v, std::move(b.v));
}

void math_vertex::set_database(const std::shared_ptr<name_server>& db) {
	if (v.properties().names == nullptr && db != nullptr) {
		if (is_arithmetic(v.properties().op)) {
			v.properties().names = db;
			for (int i = 0; i < v.get_edge_count(); i++) {
				get_edge_in(i).set_database(db);
			}
		}
	}
}

math_vertex math_vertex::binary_op(operation_t op, math_vertex A, math_vertex B) {
	properties props;
	props.op = op;
	math_vertex C;
	C.v = dag_vertex<properties>::new_(std::move(props));
	auto db = A.v.properties().names;
	db = (db == nullptr) ? B.v.properties().names : db;
	A.set_database(db);
	B.set_database(db);
	C.set_database(db);
	C.v.add_edge_in(A.v);
	C.v.add_edge_in(B.v);
	C = C.optimize();
	C.check_cse();
	return std::move(C);
}

math_vertex math_vertex::unary_op(operation_t op, math_vertex A) {
	properties props;
	props.op = op;
	math_vertex C;
	C.v = dag_vertex<properties>::new_(std::move(props));
	auto db = A.v.properties().names;
	A.set_database(db);
	C.set_database(db);
	C.v.add_edge_in(A.v);
	C = C.optimize();
	C.check_cse();
	return std::move(C);
}

math_vertex::~math_vertex() {

}

math_vertex::math_vertex(double constant) {
	if (constant >= 0.0) {
		auto iter = consts.upper_bound(constant);
		int flag = false;
		if (iter != consts.end()) {
			if (close2(iter->first, constant)) {
				flag = true;
			} else if (iter != consts.begin()) {
				iter--;
				if (close2(iter->first, constant)) {
					flag = true;
				}
			}
		}
		if (!flag) {
			properties props;
			props.op = CON;
			props.value = constant;
			props.name = std::make_shared<std::string>(std::to_string(constant));
			props.names = nullptr;
			v = dag_vertex<properties>::new_(std::move(props));
			consts[constant] = *this;
		} else {
			*this = iter->second;
		}
	} else {
		*this = -math_vertex(-constant);
	}
	check_cse();
}

math_vertex& math_vertex::operator=(double constant) {
	*this = math_vertex(constant);
	return *this;
}

math_vertex math_vertex::new_input(std::shared_ptr<name_server> db, std::string&& name) {
	properties props;
	math_vertex C;
	props.op = IN;
	auto tmp = db->reserve_name(std::move(name));
	props.name = tmp.first;
	assert(tmp.second == "");
	props.names = db;
	C.v = dag_vertex<properties>::new_(std::move(props));
	C.check_cse();
	return std::move(C);
}

std::vector<math_vertex> math_vertex::new_inputs(int cnt) {
	std::vector<math_vertex> inputs;
	auto db = std::make_shared<name_server>();
	for (int i = 0; i < cnt; i++) {
		inputs.push_back(new_input(db, "x[" + std::to_string(i) + "]"));
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

std::string math_vertex::execute_all(std::vector<math_vertex>& outputs) {
	std::string code;
	dag_vertex<properties>::executor exe;
	for (int n = 0; n < outputs.size(); n++) {
		outputs[n].v.properties().out_num = n;
	}
	auto db = outputs[0].v.properties().names;
	for (auto& o : outputs) {
		code += o.execute(exe);
	}
	auto decls = db->get_declarations();
	code = decls + code;
	return std::move(code);
}

math_vertex::math_vertex(const dag_vertex<properties>& v0) {
	v = v0;
}

math_vertex::math_vertex(dag_vertex<properties> && v0) {
	v = std::move(v0);
}

math_vertex math_vertex::get_edge_in(int i) const {
	return math_vertex(v.get_edge_in(i));
}

int math_vertex::get_edge_count() const {
	return v.get_edge_count();
}

bool math_vertex::check_cse() {
	value_number vn;
	const auto op = v.properties().op;
	vn.op = op;
	if (is_additive(op)) {
		vn.x = associative_adds();
	} else if (op == MUL) {
		vn.x = associative_muls();
	} else {
		vn.x.insert(v.get_unique_id());
	}
	if (cse.find(vn) == cse.end()) {
		value_number* ptr = new value_number(vn);
		v.properties().vnum = std::shared_ptr<value_number>(ptr, [](value_number* vptr) {
			cse.erase(*vptr);
			delete vptr;
		});
		cse[vn] = *this;
		return false;
	} else {
		fprintf(stderr, "%s\n", vn.to_string().c_str());
		*this = cse[vn];
		return true;
	}
}

bool math_vertex::is_zero() const {
	if (v.properties().op == CON) {
		return close2(v.properties().value, 0.0);
	}
	return false;
}

bool math_vertex::is_one() const {
	if (v.properties().op == CON) {
		return close2(v.properties().value, 1.0);
	}
	return false;
}

bool math_vertex::is_neg_one() const {
	if (v.properties().op == CON) {
		return close2(v.properties().value, -1.0);
	}
	return false;
}

bool math_vertex::is_neg() const {
	return v.properties().op == NEG;
}

math_vertex math_vertex::get_neg() const {
	if (is_neg()) {
		return get_edge_in(0);
	} else {
		return *this;
	}
}

math_vertex math_vertex::optimize() {
	const auto op = v.properties().op;
	for (int ei = 0; ei < edge_count(op); ei++) {
		replace_edge(get_edge_in(ei), get_edge_in(ei).optimize());
	}
	math_vertex c = *this;
	math_vertex a;
	math_vertex b;
	if (edge_count(op) >= 1) {
		a = get_edge_in(0);
	}
	if (edge_count(op) >= 2) {
		b = get_edge_in(1);
	}
	switch (op) {
	case ADD:
		if (a.is_neg() && b.is_neg()) {
			c = -(a.get_neg() + b.get_neg());
		} else if (a.is_neg()) {
			c = b - a.get_neg();
		} else if (b.is_neg()) {
			c = a - b.get_neg();
		} else if (a.is_zero()) {
			c = b;
		} else if (b.is_zero()) {
			c = a;
		}
		break;
	case SUB:
		if (a.is_neg() && b.is_neg()) {
			c = b.get_neg() - a.get_neg();
		} else if (a.is_neg()) {
			c = -(a.get_neg() + b);
		} else if (b.is_neg()) {
			c = a + b.get_neg();
		} else if (a.is_zero()) {
			c = -b;
		} else if (b.is_zero()) {
			c = a;
		}
		break;
	case MUL:
		if (a.is_neg() && b.is_neg()) {
			c = a.get_neg() * b.get_neg();
		} else if (a.is_neg()) {
			c = -(a.get_neg() * b);
		} else if (b.is_neg()) {
			c = -(b.get_neg() * a);
		} else if (a.is_zero() || b.is_zero()) {
			c = 0.0;
		} else if (a.is_one()) {
			c = b;
		} else if (b.is_one()) {
			c = a;
		} else if (a.is_neg_one()) {
			c = -b;
		} else if (b.is_neg_one()) {
			c = -a;
		}
		break;
	case NEG:
		if (a.is_neg()) {
			c = a.get_neg();
		}
		break;
	default:
		break;
	}

	return c;
}

bool math_vertex::operator<(const math_vertex& other) const {
	return v < other.v;
}

bool math_vertex::operator==(const math_vertex& other) const {
	return v == other.v;
}

assoc_set math_vertex::associative_adds() const {
	assoc_set adds;
	const auto op = v.properties().op;
	if (is_additive(op)) {
		adds = get_edge_in(0).associative_adds();
		if (op == ADD) {
			adds = adds + get_edge_in(0).associative_adds();
		} else {
			adds = adds - get_edge_in(1).associative_adds();
		}
	} else if (op == NEG) {
		adds = -get_edge_in(0).associative_adds();
	} else {
		adds.insert(v.get_unique_id());
	}
	return std::move(adds);
}

assoc_set math_vertex::associative_muls() const {
	assoc_set muls;
	const auto op = v.properties().op;
	if (op == MUL) {
		muls = get_edge_in(0).associative_muls() + get_edge_in(1).associative_muls();
	} else {
		muls.insert(v.get_unique_id());
	}
	return std::move(muls);
}

std::unordered_map<math_vertex::value_number, math_vertex::weak_ref, math_vertex::value_key> math_vertex::cse;
std::map<double, math_vertex> math_vertex::consts;
