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
	return (op == other.op) && (x == other.x) && (sgn == other.sgn);
}

math_vertex::value_number math_vertex::value_number::operator-() const {
	value_number n;
	n.op = op;
	n.sgn = -n.sgn;
	n.x = -x;
	return n;
}

math_vertex::value_number::value_number() {
	op = ADD;
	sgn = 0;
}

size_t math_vertex::value_key::operator()(const value_number& value) const {
	int rc = (int) value.op;
	rc = rc ^ (value.sgn * 1664525 + 1013904223);
	for (auto i = value.x.begin(); i != value.x.end(); i++) {
		rc = rc ^ (i->first * 1664525 + 1013904223);
	}
	for (auto i = value.x.begin(); i != value.x.end(); i++) {
		rc = rc ^ (i->second * 1664525 + 1013904223);
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
	cse = false;
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
			for (int i = 0; i < v.get_edge_in_count(); i++) {
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
	C = C.distribute_muls();
	C = C.optimize();
	C.check_cse();
	return std::move(C);
}

math_vertex math_vertex::post_optimize() {
	for (int i = 0; i < get_edge_in_count(); i++) {
		auto new_edge = get_edge_in(i).post_optimize();
		v.replace_edge_in(get_edge_in(i).v, std::move(new_edge.v));
	}
	math_vertex rc;
	if (is_additive(get_op())) {
		rc = associate_adds();
	} else {
		rc = *this;
	}
	return rc;
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
	C = C.distribute_muls();
	C = C.optimize();
	C.check_cse();
	return std::move(C);
}


void math_vertex::optimize(std::vector<math_vertex>& vertices) {
	for (auto& v : vertices) {
		v = v.post_optimize();
	}
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

math_vertex::op_cnt_t math_vertex::operation_count(dag_vertex<properties>::executor& exe) {
	op_cnt_t cnt;
	cnt.add = cnt.mul = cnt.neg = 0;
	v.execute(exe, [&cnt](properties& in, const std::vector<properties>& edges) {
		if(is_additive(in.op)) {
			cnt.add++;
		} else if( in.op == MUL) {
			cnt.mul++;
		} else if( in.op == NEG) {
			cnt.neg++;
		}
	});
	return cnt;
}

math_vertex::op_cnt_t math_vertex::operation_count(std::vector<math_vertex>& outputs) {
	op_cnt_t cnt;
	cnt.add = cnt.mul = cnt.neg = 0;
	dag_vertex<properties>::executor exe(false);
	for (int n = 0; n < outputs.size(); n++) {
		auto this_cnt = outputs[n].operation_count(exe);
		cnt.add += this_cnt.add;
		cnt.mul += this_cnt.mul;
		cnt.neg += this_cnt.neg;
	}
	return cnt;
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

int math_vertex::get_edge_in_count() const {
	return v.get_edge_in_count();
}

math_vertex math_vertex::get_edge_out(int i) const {
	return math_vertex(v.get_edge_out(i));
}

int math_vertex::get_edge_out_count() const {
	return v.get_edge_out_count();
}

bool math_vertex::check_cse() {
	value_number vn;
	const auto op = v.properties().op;
	vn.op = op;
	auto add_set = associative_adds();
	if (is_additive(op)) {
		vn.x = add_set;
	} else if (op == MUL) {
		vn.x = associative_muls().first;
	} else {
		vn.x.insert(v.get_unique_id());
	}
	bool pos = false;
	bool neg = false;
	auto piter = cse.find(vn);
	decltype(piter) niter;
	math_vertex vacate;
	if (piter != cse.end()) {
		if (piter->second.vacate) {
			vacate = piter->second.ptr;
		} else {
			pos = true;
		}
	} else {
		niter = cse.find(-vn);
		if (niter != cse.end()) {
			if (niter->second.vacate) {
				vacate = niter->second.ptr;
			} else {
				neg = true;
			}
		}
	}
	if (vacate.v != math_vertex().v) {
		vacate.v.properties().cse = false;
	}
	if (!pos && !neg) {
		value_number* ptr = new value_number(vn);
		v.properties().cse = true;
		bool* cse_ptr = &(v.properties().cse);
		v.properties().vnum = std::shared_ptr<value_number>(ptr, [cse_ptr](value_number* vptr) {
			if( *cse_ptr ) {
				cse.erase(*vptr);
			}
			delete vptr;
		});
		cse[vn] = *this;
		return false;
	} else {
		if (pos) {
			*this = cse[vn].ptr;
		} else {
			*this = -cse[-vn].ptr;
		}
		return true;
	}
}

math_vertex::cse_entry& math_vertex::cse_entry::operator=(const math_vertex& v) {
	ptr = v;
	return *this;
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

math_vertex math_vertex::associate_adds() {


	return optimize();


	math_vertex rc = *this;
	if (is_additive(get_op())) {
		std::unordered_set<math_vertex, key> path;
		const auto basis_vertices = collect_additive_terms(path);
		std::vector<assoc_set> basis;
		std::vector<math_vertex> basis_vs;
		std::vector<int> basis_sgns;
		auto target = v.properties().vnum->x;
		for (auto v : basis_vertices) {
			auto vnum = (v.v.properties().vnum);
			if (vnum) {
				if (!is_additive(v.get_op())) {
					vnum->x.clear();
					vnum->x.insert(v.v.get_unique_id());
				}
				if (vnum->x.is_subset_of(target)) {
					basis.push_back(vnum->x);
					basis.push_back(-vnum->x);
					basis_vs.push_back(v);
					basis_vs.push_back(v);
					basis_sgns.push_back(1);
					basis_sgns.push_back(-1);
				}
			}
		}
		std::vector<int> best_set;
		int sgn = 0;
		int n = basis.size();
		for (auto v : path) {
			if (v.v.properties().vnum) {
				cse[*v.v.properties().vnum].vacate = true;
			}
		}
		for (int k = 1; k < n; k++) {
			const auto tries = nchoosek(n, k);
			for (auto try_ : tries) {
				assoc_set set;
				for (auto t : try_) {
					set = set + basis[t];
				}
				if (set == target) {
					best_set = std::move(try_);
					sgn = 1;
					break;
				} else if (set == -target) {
					sgn = -1;
					best_set = std::move(try_);
					break;
				}
			}
			if (sgn) {
				break;
			}
		}
		if (sgn && best_set.size() < target.size()) {
			fprintf( stderr, "\n\n---------------------------------\n\n");
			fprintf( stderr, "%s\n", v.properties().vnum->to_string().c_str());
			fprintf( stderr, "---------------------------------\n");
			for (auto b : best_set) {
				fprintf( stderr, "%s\n", basis[b].to_string().c_str());
			}
			fprintf( stderr, "---------------------------------\n\n");
			rc = 0.0;
			for (auto i : best_set) {
				auto this_sgn = sgn * basis_sgns[i];
				if (this_sgn > 0) {
					rc = rc + basis_vs[i];
				} else {
					rc = rc - basis_vs[i];
				}
			}
		}
	}
	return rc;
}

std::unordered_set<math_vertex, math_vertex::key> math_vertex::collect_additive_terms(std::unordered_set<math_vertex, key>& path) {
	return collect_additive_terms_down(path);
}

size_t math_vertex::key::operator()(const math_vertex& v) const {
	return v.v.get_unique_id();
}

std::unordered_set<math_vertex, math_vertex::key> math_vertex::collect_additive_terms_down(std::unordered_set<math_vertex, key>& path) {
	std::unordered_set<math_vertex, key> terms;
	if (!is_additive(get_op())) {
		terms = collect_additive_terms_up(path);
		terms.insert(*this);
		for (int i = 0; i < get_edge_out_count(); i++) {
			if (path.find(get_edge_out(i)) == path.end()) {
				auto tmp = get_edge_out(i).collect_additive_terms_up(path);
				terms.insert(tmp.begin(), tmp.end());
			}
		}
	} else {
		for (int i = 0; i < get_edge_in_count(); i++) {
			path.insert(v);
			auto tmp = get_edge_in(i).collect_additive_terms_down(path);
			terms.insert(tmp.begin(), tmp.end());
		}
	}
	return std::move(terms);
}

std::unordered_set<math_vertex, math_vertex::key> math_vertex::collect_additive_terms_up(std::unordered_set<math_vertex, key>& path) {
	std::unordered_set<math_vertex, key> terms;
	if (is_additive(get_op())) {
		terms.insert(*this);
		for (int i = 0; i < get_edge_out_count(); i++) {
			if (path.find(get_edge_out(i)) == path.end()) {
				auto tmp = get_edge_out(i).collect_additive_terms_up(path);
				terms.insert(tmp.begin(), tmp.end());
			}
		}
	}
	return std::move(terms);
}

math_vertex math_vertex::distribute_muls() {
	auto dmuls = distributive_muls();
	std::unordered_map<double, std::vector<dag_vertex<properties>>>map_add;
	std::unordered_map<double, std::vector<dag_vertex<properties>>>map_sub;
	std::set<double, std::greater<double>> coeffs;
	std::map<dag_vertex<properties>, double> cosums;
	for (const auto& d : dmuls) {
		if (cosums.find(d.v) != cosums.end()) {
			cosums[d.v] = 0.0;
		}
		cosums[d.v] += d.c;
	}
	for (auto i = cosums.begin(); i != cosums.end(); i++) {
		if (close2(i->second, 0.0)) {
			continue;
		} else if (i->second > 0.0) {
			map_add[i->second].push_back(i->first);
		} else {
			map_sub[-i->second].push_back(i->first);
		}
		coeffs.insert(consts[std::abs(i->second)].get_value());
	}
	int max_cnt = 0;
	bool flag = false;
	for (auto c : coeffs) {
		if (!close2(c, 1.0)) {
			const int cnt = map_add[c].size() + map_sub[c].size();
			max_cnt = std::max(cnt, max_cnt);
		}
		if (max_cnt >= 2) {
			flag = true;
			break;
		}
	}
	if (flag) {
		math_vertex result = 0.0;
		for (auto c : coeffs) {
			auto c0 = math_vertex(c);
			math_vertex sum = 0.0;
			for (int i = 0; i < map_add[c].size(); i++) {
				sum = sum + map_add[c][i];
			}
			for (int i = 0; i < map_sub[c].size(); i++) {
				sum = sum - map_sub[c][i];
			}
			if (close2(c, 1.0)) {
				result = result + sum;
			} else {
				result = result + c0 * sum;
			}
		}
		return result;
	}
	return *this;
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
		if (a == b) {
			c = math_vertex(2.0) * a;
		} else if (a.is_neg() && b.is_neg()) {
			c = -(a.get_neg() + b.get_neg());
		} else if (a.is_neg()) {
			c = b - a.get_neg();
		} else if (b.is_neg()) {
			c = a - b.get_neg();
		} else if (a.is_zero()) {
			c = b;
		} else if (b.is_zero()) {
			c = a;
		} else if (a.get_op() == CON && b.get_op() == CON) {
			c = math_vertex(a.get_value() + b.get_value());
		}
		break;
	case SUB:
		if (a == b) {
			c = math_vertex(0.0);
		} else if (a.is_neg() && b.is_neg()) {
			c = b.get_neg() - a.get_neg();
		} else if (a.is_neg()) {
			c = -(a.get_neg() + b);
		} else if (b.is_neg()) {
			c = b.get_neg() + a;
		} else if (a.is_zero()) {
			c = -b;
		} else if (b.is_zero()) {
			c = a;
		} else if (a.get_op() == CON && b.get_op() == CON) {
			c = math_vertex(a.get_value() - b.get_value());
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
		} else if (a.get_op() == CON && b.get_op() == CON) {
			c = math_vertex(a.get_value() * b.get_value());
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
			adds = adds + get_edge_in(1).associative_adds();
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

operation_t math_vertex::get_op() const {
	return v.properties().op;
}

double math_vertex::get_value() const {
	return v.properties().value;
}

std::vector<math_vertex::distrib_t> math_vertex::distributive_muls() {
	std::vector<math_vertex::distrib_t> muls;
	const auto op = v.properties().op;
	if (is_additive(op)) {
		auto tmp = get_edge_in(0).distributive_muls();
		muls.insert(muls.end(), tmp.begin(), tmp.end());
		if (op == ADD) {
			auto tmp = get_edge_in(1).distributive_muls();
			muls.insert(muls.end(), tmp.begin(), tmp.end());
		} else {
			auto tmp = get_edge_in(1).distributive_muls();
			for (auto& t : tmp) {
				t.c = -t.c;
			}
			muls.insert(muls.end(), tmp.begin(), tmp.end());
		}
	} else if (op == NEG) {
		muls = get_edge_in(0).distributive_muls();
		for (auto& t : muls) {
			t.c = -t.c;
		}
	} else if (op == MUL) {
		auto A = get_edge_in(0);
		auto B = get_edge_in(1);
		if (B.get_op() == CON) {
			std::swap(A, B);
		}
		distrib_t d;
		if (A.get_op() == CON) {
			d.c = A.get_value();
			d.v = B.v;
		} else {
			d.c = 1.0;
			d.v = v;
		}
		muls.push_back(d);
	} else {
		distrib_t d;
		d.c = 1.0;
		d.v = v;
		muls.push_back(d);
	}
	return std::move(muls);
}

math_vertex::math_vertex() {
	if( !first_init) {
		first_init = true;
		essential_constants.resize(3);
		essential_constants[0] = (math_vertex(0.0));
		essential_constants[1] = (math_vertex(1.0));
		essential_constants[2] = (math_vertex(-1.0));
	}
}

math_vertex::cse_entry::cse_entry() {
	vacate = false;
}

std::pair<assoc_set, int> math_vertex::associative_muls() const {
	std::pair<assoc_set, int> muls;
	const auto op = v.properties().op;
	if (op == MUL) {
		const auto mulsa = get_edge_in(0).associative_muls();
		const auto mulsb = get_edge_in(1).associative_muls();
		muls.first = mulsa.first + mulsb.first;
		muls.second = mulsa.second * mulsb.second;
	} else if (op == NEG) {
		const auto mulsa = get_edge_in(0).associative_muls();
		muls.first = mulsa.first;
		muls.second = -mulsa.second;
	} else {
		muls.first.insert(v.get_unique_id());
		muls.second = 1;
	}
	return std::move(muls);
}

std::unordered_map<math_vertex::value_number, math_vertex::cse_entry, math_vertex::value_key> math_vertex::cse;
std::vector<math_vertex> math_vertex::essential_constants;
std::map<double, math_vertex> math_vertex::consts;
bool math_vertex::first_init = false;
