#include "math.hpp"
#include "util.hpp"

#include <algorithm>
#include <numeric>
#include <queue>

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
	return ptr < other.ptr;
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
	depth = 0;
	group_id = -1;
	goal = false;
}

bool math_vertex::valid() const {
	return v.valid();
}

std::unordered_map<std::string, std::string> mem2reg;
std::unordered_map<std::string, std::string> reg2mem;
std::queue<std::string> regq;
int regcnt = 0;

void clear_registers() {
	mem2reg.clear();
	reg2mem.clear();
	regq = decltype(regq)();
	regcnt = 0;
}

std::pair<std::string, std::string> get_register(std::string mem, bool noload, std::set<std::string> onlystore = std::set<std::string>()) {
	std::pair<std::string, std::string> rc;
	std::string& reg = rc.first;
	if (mem2reg.find(mem) == mem2reg.end()) {
		if (regcnt < 16) {
			regcnt = 16;
			for (int n = 0; n < 16; n++) {
				regq.push(std::string("%ymm") + std::to_string(n));
			}
		}
		reg = regq.front();
		regq.pop();
		if (reg2mem.find(reg) != reg2mem.end()) {
			if (reg2mem[reg][0] != 'C') {
				if (onlystore.find(reg2mem[reg]) == onlystore.end() || onlystore.size() == 0) {
					rc.second = std::string("               vmovupd        ") + reg + ", " + reg2mem[reg] + "\n";
				}
			}
			mem2reg.erase(reg2mem[reg]);
			reg2mem.erase(reg);
		}
		reg2mem[reg] = mem;
		mem2reg[mem] = reg;
		regq.push(reg);
		if (!noload) {
			rc.second += std::string("               vmovupd        ") + reg2mem[reg] + ", " + reg + "\n";
		}
	} else {
		reg = mem2reg[mem];
	}
	return rc;
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
	C.v.properties().depth = std::max(A.v.properties().depth, B.v.properties().depth) + 1;
	A.set_database(db);
	B.set_database(db);
	C.set_database(db);
	C.set_group_id(A.get_group_id());
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
	C.v.properties().depth = A.v.properties().depth + 1;
	auto db = A.v.properties().names;
	A.set_database(db);
	C.set_group_id(A.get_group_id());
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
			}
		}
		if (!flag && (iter != consts.begin())) {
			iter--;
			if (close2(iter->first, constant)) {
				flag = true;
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
}

math_vertex& math_vertex::operator=(double constant) {
	*this = math_vertex(constant);
	return *this;
}

math_vertex math_vertex::new_input(std::shared_ptr<name_server> db) {
	properties props;
	math_vertex C;
	props.op = IN;
	auto tmp = db->generate_name();
	props.name = tmp;
	props.names = db;
	C.v = dag_vertex<properties>::new_(std::move(props));
	return C;
}

std::vector<math_vertex> math_vertex::new_inputs(int cnt) {
	std::vector<math_vertex> inputs;
	auto db = std::make_shared<name_server>();
	for (int i = 0; i < cnt; i++) {
		inputs.push_back(new_input(db));
	}
	return inputs;
}

math_vertex operator+(const math_vertex& A, const math_vertex& B) {
	assert(A.valid());
	assert(B.valid());
	return math_vertex::binary_op(ADD, A, B);
}

math_vertex operator-(const math_vertex& A, const math_vertex& B) {
	assert(A.valid());
	assert(B.valid());
	return math_vertex::binary_op(SUB, A, B);
}

math_vertex operator*(const math_vertex& A, const math_vertex& B) {
	assert(A.valid());
	assert(B.valid());
	return math_vertex::binary_op(MUL, A, B);
}

math_vertex operator-(const math_vertex& A) {
	assert(A.valid());
	return math_vertex::unary_op(NEG, A);
}

math_vertex& math_vertex::operator+=(const math_vertex& other) {
	assert(valid());
	*this = *this + other;
	return *this;
}

math_vertex& math_vertex::operator-=(const math_vertex& other) {
	assert(valid());
	*this = *this - other;
	return *this;
}

math_vertex& math_vertex::operator*=(const math_vertex& other) {
	assert(valid());
	*this = *this * other;
	return *this;
}

std::vector<math_vertex> math_vertex::available2execute(dag_vertex<properties>::executor& exe1, dag_vertex<properties>::executor& exe2, std::vector<math_vertex>& outputs) {
	std::vector<math_vertex> rc;
	for (auto o : outputs) {
		std::vector<math_vertex> tmp = o.available2execute(exe1, exe2);
		rc.insert(rc.end(), tmp.begin(), tmp.end());
	}
	return rc;
}

std::vector<math_vertex> math_vertex::available2execute(dag_vertex<properties>::executor& exe1, dag_vertex<properties>::executor& exe2) {
	std::vector<math_vertex> rc;
	if (exe1.touched.find(v.get_unique_id()) == exe1.touched.end()) {
		auto n = get_edge_in_count();
		for (int i = 0; i < n; i++) {
			auto tmp = get_edge_in(i).available2execute(exe1, exe2);
			rc.insert(rc.end(), tmp.begin(), tmp.end());
		}
		if (!rc.size()) {
			if (exe2.touched.find(v.get_unique_id()) == exe2.touched.end()) {
				rc.push_back(*this);
			}
		}
		exe1.touched.insert(v.get_unique_id());
	}
	return rc;
}

std::string math_vertex::execute(dag_vertex<properties>::executor& exe) {

	std::string code;
	return code;
}

math_vertex::op_cnt_t math_vertex::operation_count(dag_vertex<properties>::executor& exe) {
	assert(valid());
	op_cnt_t cnt;
	cnt.add = cnt.mul = cnt.neg = 0;
	v.execute(exe, [&cnt](dag_vertex<properties>& in, const std::vector<dag_vertex<properties>>& edges) {
		if(is_additive(in.properties().op)) {
			cnt.add++;
		} else if( in.properties().op == MUL) {
			cnt.mul++;
		} else if( in.properties().op == NEG) {
			cnt.neg++;
		}
	});
	return cnt;
}

math_vertex::op_cnt_t math_vertex::operation_count(std::vector<math_vertex> outputs) {
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

math_vertex::op_cnt_t math_vertex::operation_count(std::vector<cmplx> x) {
	std::vector<math_vertex> X;
	for (int i = 0; i < x.size(); i++) {
		X.push_back(x[i].x);
		X.push_back(x[i].y);
	}
	return operation_count(X);
}

void math_vertex::set_group_id(int id) {
	v.properties().group_id = id;
}

math_vertex::operator dag_vertex<math_vertex::properties>() const {
	return v;
}

int math_vertex::get_group_id() const {
	return v.properties().group_id;
}

std::pair<std::string, int> math_vertex::execute_all(std::vector<math_vertex>&& inputs, std::vector<math_vertex>& outputs) {
	clear_registers();
	std::string code;
	int ncnt = 5;
	std::pair<std::string, int> rc;
	std::vector<std::string> onames;
	for (int n = 0; n < outputs.size(); n++) {
		outputs[n].v.properties().out_num = n;
		outputs[n].set_goal();
	}
	auto db = outputs[0].v.properties().names;

	std::string loads;
	for (int i = 0; i < inputs.size(); i++) {
		auto reg = get_register(*inputs[i].v.properties().name, true);
		onames.push_back(std::to_string(-32 * (1 + i)) + std::string("(%ebp)"));
		loads += reg.second;
		char* ptr;
		const char* basename = i % 2 == 0 ? "%rdi" : "%rsi";
		if (i == 0) {
			asprintf(&ptr, "%15s%-15s%s, %s\n", "", "vmovupd", (std::string("(") + basename + ")").c_str(), reg.first.c_str());
			loads += ptr;
			free(ptr);
		} else {
			static bool sw = true;
			const char* snm = sw ? "%rax" : "%r8";
			sw = !sw;
			asprintf(&ptr, "%15s%-15s%s, %s, %s\n", "", "imul", (std::string("$") + std::to_string(4 * i)).c_str(), "%rdx", snm);
			loads += ptr;
			free(ptr);
			asprintf(&ptr, "%15s%-15s%s, %s\n", "", "vmovupd", (std::string("(") + basename + ", " + snm + ", 8)").c_str(), reg.first.c_str());

			loads += ptr;
			free(ptr);
		}

	}
	dag_vertex<properties>::executor exe;
	std::vector<math_vertex::weak_ref> nodes;
	std::vector<dag_vertex<properties>> dags;
	for (auto o : outputs) {
		dags.push_back(o.v);
	}
	dags = dag_vertex<properties>::sort(exe, dags);
	std::set<math_vertex::weak_ref> done;
	for (auto d : dags) {
		nodes.push_back(math_vertex(d));
	}
	dags.clear();

	std::stack<math_vertex> S;
	std::unordered_set<int> visited;
	std::map<int, math_vertex> goals;
	for (auto n : nodes) {
		auto node = math_vertex(n);
		auto id = node.get_unique_id();
		if (node.v.properties().goal) {
			goals[id] = node;
			visited.insert(id);
		}
		if (node.get_op() == IN) {
			visited.insert(id);
		}
	}
	if (goals.size()) {
		auto i = goals.end();
		do {
			i--;
			S.push(i->second);
		} while (i != goals.begin());
	}
	nodes.resize(0);
	goals.clear();
	while (S.size()) {
		auto V = S.top();
		bool flag = true;
		for (int i = 0; i < V.get_edge_in_count(); i++) {
			auto U = V.get_edge_in(i);
			auto id = U.v.get_unique_id();
			if (visited.find(id) == visited.end()) {
				flag = false;
				S.push(U);
				visited.insert(id);
			}
		}
		if (flag) {
			nodes.push_back(V);
			S.pop();
		}
	}

	inputs.clear();

	int index = 0;
	std::string constants;
	std::string nm = std::string("NEG_ONE");
	char* ptr;
	asprintf(&ptr, "%-15s%-15s%-25.17e\n", (nm + ":").c_str(), ".double", -1.0);
	constants += ptr;
	free(ptr);
	asprintf(&ptr, "%15s%-15s%-25.17e\n", "", ".double", -1.0);
	constants += ptr;
	free(ptr);
	asprintf(&ptr, "%15s%-15s%-25.17e\n", "", ".double", -1.0);
	constants += ptr;
	free(ptr);
	asprintf(&ptr, "%15s%-15s%-25.17e\n", "", ".double", -1.0);
	constants += ptr;
	free(ptr);
	for (auto n : nodes) {
		auto node = math_vertex(n);
		if (node.get_op() == CON && done.find(node) == done.end()) {
			if (!(node.is_zero() || node.is_one() || node.is_neg_one())) {
				std::string nm = std::string("C") + std::to_string(index++);
				char* ptr;
				asprintf(&ptr, "%-15s%-15s%-25.17e\n", (nm + ":").c_str(), ".double", node.v.properties().value);
				constants += ptr;
				free(ptr);
				asprintf(&ptr, "%15s%-15s%-25.17e\n", "", ".double", node.v.properties().value);
				constants += ptr;
				free(ptr);
				asprintf(&ptr, "%15s%-15s%-25.17e\n", "", ".double", node.v.properties().value);
				constants += ptr;
				free(ptr);
				asprintf(&ptr, "%15s%-15s%-25.17e\n", "", ".double", node.v.properties().value);
				constants += ptr;
				free(ptr);
				node.v.properties().name = std::make_shared<std::string>(nm);
			}
			done.insert(node);
		}
	}
	for (auto n : nodes) {
		auto v = math_vertex(n);
		if (done.find(v) != done.end()) {
			continue;
		}
		int ncnt = v.get_edge_in_count();
		for (int j = 0; j < ncnt; j++) {
			auto edge = v.get_edge_in(j);
			if (edge.v.use_count() <= 2 && !edge.is_constant()) {
				v.v.properties().name = edge.v.properties().name;
				break;
			}
		}
		std::vector<math_vertex> in;
		for (int k = 0; k < v.get_edge_in_count(); k++) {
			in.push_back(v.get_edge_in(k));
		}
		code += v.v.properties().print_code(in);
		v.v.free_edges();
		done.insert(v);
	}
	std::vector<int> sorted(onames.size());
	std::iota(sorted.begin(), sorted.end(), 0);
	std::sort(sorted.begin(), sorted.end(), [&](int i, int j) {
		auto a = onames[i];
		auto b = onames[j];
		return mem2reg.find(a) != mem2reg.end() && mem2reg.find(b) == mem2reg.end();
	});
	std::set<std::string> onlystore(onames.begin(), onames.end());
	for (int k = 0; k < onames.size(); k++) {
		const int i = sorted[k];
		char* ptr;
		auto reg = get_register(onames[i], false, onlystore);
		code += reg.second;
		const char* basename = i % 2 == 0 ? "%rdi" : "%rsi";
		static bool sw = false;
		sw = !sw;
		const char* snm = sw ? "%rax" : "%r8";
		std::string memloc;
		if (i == 0) {
			memloc = std::string("(") + basename + ")";
		} else {
			memloc = std::string("(") + basename + ", " + snm + ", 8)";
		}
		asprintf(&ptr, "%15s%-15s%s, %s, %s\n", "", "imul", (std::string("$") + std::to_string(4 * i)).c_str(), "%rdx", snm);
		code += ptr;
		free(ptr);
		asprintf(&ptr, "%15s%-15s%s, %s\n", "", "vmovupd", reg.first.c_str(), memloc.c_str());
		code += ptr;
		free(ptr);
	}

	auto decls = db->get_declarations();
	code = decls + loads + code;
	rc.first = code;
	rc.second = 0;
	for (int i = 0; i < decls.size(); i++) {
		if (decls[i] == '\n') {
			rc.second++;
		}
	}
	//rc.first += clear_registers();
	rc.first += "               ret\n";
	rc.first += "               .align         32\n";
	rc.first += constants;
	return std::move(rc);
}

std::string math_vertex::properties::print_code(const std::vector<math_vertex>& edges) {
	std::string code;
//	if (out_num != -1) {
//		std::string nm = std::string("x[") + std::to_string(out_num) + std::string("]");
//		auto tmp = names->reserve_name(std::move(nm));
//		name = tmp.first;
//		code += tmp.second;
	if (name == nullptr) {
		assert(names != nullptr);
		name = names->generate_name();
	}
	if (is_arithmetic(op)) {
		const auto a = get_register(*name, true);
		const auto b = get_register(*edges[0].v.properties().name, false);
//		const auto c = get_register(*edges[op == NEG ? 0 : 1].v.properties().name);
		const auto A = a.first;
		const auto B = b.first;
		auto C = *edges[op == NEG ? 0 : 1].v.properties().name;
		if (mem2reg.find(C) != mem2reg.end()) {
			C = mem2reg[C];
		}
		char* ptr;
		if (a.second != "") {
			code += a.second;
		}
		if (b.second != "") {
			code += b.second;
		}
		switch (op) {
		case ADD:
			asprintf(&ptr, "%15s%-15s%s, %s, %s\n", "", "vaddpd", C.c_str(), B.c_str(), A.c_str());
			break;
		case SUB:
			asprintf(&ptr, "%15s%-15s%s, %s, %s\n", "", "vsubpd", C.c_str(), B.c_str(), A.c_str());
			break;
		case MUL:
			asprintf(&ptr, "%15s%-15s%s, %s, %s\n", "", "vmulpd", C.c_str(), B.c_str(), A.c_str());
			break;
		case NEG:
			asprintf(&ptr, "%15s%-15s%s, %s, %s\n", "", "vmulpd", "NEG_ONE", B.c_str(), A.c_str());
			break;
		default:
			break;
		}
		code += ptr;
		free(ptr);
	}
	return code;
}

math_vertex::math_vertex(const dag_vertex<properties>& v0) {
	v = v0;
}

math_vertex::math_vertex(dag_vertex<properties> && v0) {
	v = std::move(v0);
}

math_vertex math_vertex::get_edge_in(int i) const {
	assert(valid());
	return math_vertex(v.get_edge_in(i));
}

int math_vertex::get_edge_in_count() const {
	assert(valid());
	return v.get_edge_in_count();
}
/*
 math_vertex math_vertex::get_edge_out(int i) const {
 assert(valid());
 return math_vertex(v.get_edge_out(i));
 }

 int math_vertex::get_edge_out_count() const {
 assert(valid());
 return v.get_edge_out_count();
 }
 */
bool math_vertex::check_cse() {
	assert(valid());
	value_number vn;
	const auto op = v.properties().op;
	if (op == CON || op == IN) {
		return false;
	}
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
	assert(!(neg && pos));
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
		auto op = get_op();
		assert((int ) op < 6);
		assert((int ) op >= 0);
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
	assert(valid());
	if (v.properties().op == CON) {
		return close2(v.properties().value, 0.0);
	}
	return false;
}

bool math_vertex::is_one() const {
	assert(valid());
	if (v.properties().op == CON) {
		return close2(v.properties().value, 1.0);
	}
	return false;
}

bool math_vertex::is_neg_one() const {
	assert(valid());
	if (v.properties().op == CON) {
		return close2(v.properties().value, -1.0);
	}
	return false;
}

bool math_vertex::is_neg() const {
	assert(valid());
	return v.properties().op == NEG;
}

math_vertex math_vertex::get_neg() const {
	assert(valid());
	if (is_neg()) {
		auto edge = get_edge_in(0);
		if (v.properties().goal) {
			edge.set_goal();
		}
		return edge;
	} else {
		return *this;
	}
}

size_t math_vertex::key::operator()(const math_vertex& v) const {
	return v.v.get_unique_id();
}

bool math_vertex::is_constant() const {
	return get_op() == CON;
}

math_vertex math_vertex::optimize() {
	assert(valid());
	auto op = v.properties().op;
	math_vertex c = *this;
	math_vertex a;
	math_vertex b;
	if (edge_count(op) >= 1) {
		a = get_edge_in(0);
	}
	if (edge_count(op) >= 2) {
		b = get_edge_in(1);
	}
	op = c.v.properties().op;
	assert(c.valid());
	if (is_arithmetic(op)) {
		assert(a.valid());
		assert(c.v.get_edge_in_count());
		if (op != NEG) {
			assert(b.valid());
			assert(c.v.get_edge_in_count() == 2);
		}
	}
	bool flag = false;
	int opt = 0;
	switch (op) {
	case ADD:
		if (a.is_constant() && b.is_constant()) {
			c = math_vertex(a.get_value() + b.get_value());
			flag = true;
			opt = 1;
		} else if (a == b) {
			c = a * math_vertex(2.0);
			flag = true;
			opt = 2;
		} else if (a.is_neg() && b.is_neg()) {
			c = -(a.get_neg() + b.get_neg());
			flag = true;
			opt = 3;
		} else if (a.is_neg()) {
			c = b - a.get_neg();
			flag = true;
			opt = 4;
		} else if (b.is_neg()) {
			c = a - b.get_neg();
			flag = true;
			opt = 5;
		} else if (a.is_zero()) {
			c = b;
			flag = true;
			opt = 6;
		} else if (b.is_zero()) {
			c = a;
			flag = true;
			opt = 7;
		}
		break;
	case SUB:
		if (a.is_constant() && b.is_constant()) {
			c = math_vertex(a.get_value() - b.get_value());
			flag = true;
			opt = 8;
		} else if (a == b) {
			c = math_vertex(0.0);
			flag = true;
			opt = 9;
		} else if (a.is_neg() && b.is_neg()) {
			c = b.get_neg() - a.get_neg();
			flag = true;
			opt = 10;
		} else if (a.is_neg()) {
			c = -(a.get_neg() + b);
			flag = true;
			opt = 11;
		} else if (b.is_neg()) {
			c = b.get_neg() + a;
			flag = true;
			opt = 12;
		} else if (a.is_zero()) {
			c = -b;
			flag = true;
			opt = 13;
		} else if (b.is_zero()) {
			c = a;
			flag = true;
			opt = 14;
		}
		break;
	case MUL:
		if (a.is_constant() && !b.is_constant()) {
			return b * a;
		} else if (a.is_constant() && b.is_constant()) {
			c = math_vertex(a.get_value() * b.get_value());
			flag = true;
			opt = 15;
		} else if (a.is_zero() || b.is_zero()) {
			c = 0.0;
			flag = true;
			opt = 16;
		} else if (a.is_one()) {
			c = b;
			flag = true;
			opt = 17;
		} else if (b.is_one()) {
			c = a;
			flag = true;
			opt = 18;
		} else if (a.is_neg_one()) {
			c = -b;
			flag = true;
			opt = 19;
		} else if (b.is_neg_one()) {
			c = -a;
			flag = true;
			opt = 20;
		} else if (a.is_neg() && b.is_neg()) {
			c = a.get_neg() * b.get_neg();
			flag = true;
			opt = 21;
		} else if (b.is_neg()) {
			c = -(b.get_neg() * a);
			flag = true;
			opt = 22;
		} else if (a.is_neg()) {
			c = -(a.get_neg() * b);
			flag = true;
			opt = 23;
		}
		break;
	case NEG:
		if (a.is_neg()) {
			//		c = a.get_neg();
			//		flag = true;
			opt = 24;
		}
		break;
	default:
		break;
	}
	if (!flag) {
		opt = 26;
		if (is_additive(op)) {
			if (a.get_op() == MUL && b.get_op() != MUL) {
				std::swap(a, b);
			}
			if (b.get_op() == MUL) {
				auto b0 = b.get_edge_in(0);
				auto b1 = b.get_edge_in(1);
				if (b1.get_op() == CON) {
					std::swap(b0, b1);
				}
				if (a.get_op() != MUL) {
					if (a.get_op() != CON && b0.get_op() == CON && b1.get_op() != CON) {
						if (a == b1) {
							if (op == ADD) {
								c = (1.0 + b0.get_value()) * a;
							} else {
								c = (1.0 - b0.get_value()) * a;
							}
						}
					}
				} else {
					auto a0 = a.get_edge_in(0);
					auto a1 = a.get_edge_in(1);
					if (a1.get_op() == CON) {
						std::swap(a0, a1);
					}
					if (a0.get_op() == CON && b0.get_op() == CON && a1.get_op() != CON && b1.get_op() != CON) {
						if (a1 == b1) {
							if (op == ADD) {
								c = (a0.get_value() + b0.get_value()) * a1;
							} else {
								c = (a0.get_value() - b0.get_value()) * a1;
							}
						} else if (a0.get_op() == CON && b0.get_op() == CON) {
							if ((close2(a0.get_value(), b0.get_value()) && (op == ADD)) && (close2(a0.get_value(), -b0.get_value()) && (op == SUB))) {
								c = a0.get_value() * (a1 + b1);
							} else if ((close2(a0.get_value(), -b0.get_value()) && (op == SUB)) && (close2(a0.get_value(), b0.get_value()) && (op == ADD))) {
								c = a0.get_value() * (a1 - b1);
							}
						}
					}
				}
			}
		}
	}
	op = c.v.properties().op;
	assert(c.valid());
	if (is_arithmetic(op)) {
		assert(a.valid());
		assert(c.v.get_edge_in_count());
		if (op != NEG) {
			assert(b.valid());
			assert(c.v.get_edge_in_count() == 2);
		}
	}
	if (flag) {
		c = c.optimize();
	}
	return c;
}

bool math_vertex::operator<(const math_vertex& other) const {
	assert(valid());
	return v < other.v;
}

bool math_vertex::operator==(const math_vertex& other) const {
	assert(valid());
	return v == other.v;
}

assoc_set math_vertex::associative_adds() const {
	assert(valid());
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
	assert(valid());
	return v.properties().op;
}

double math_vertex::get_value() const {
	assert(valid());
	return v.properties().value;
}

math_vertex::math_vertex() {
	if (!first_init) {
		first_init = true;
		essential_constants.resize(3);
		essential_constants[0] = (math_vertex(0.0));
		essential_constants[1] = (math_vertex(1.0));
		essential_constants[2] = (math_vertex(-1.0));
	}
}

void math_vertex::print_cse() {
	for (auto i = cse.begin(); i != cse.end(); i++) {
		math_vertex v = i->second.ptr;
		assert(v.valid());
		printf("-----------------------\n");
		printf("%i %i %i\n", i->second.ptr.use_count(), i->first.op, (int) v.get_op());
	}
}

math_vertex::cse_entry::cse_entry() {
	vacate = false;
}

std::pair<assoc_set, int> math_vertex::associative_muls() const {
	assert(valid());
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

/*void math_vertex::reset() {
 cse.clear();
 essential_constants.clear();
 consts.clear();
 first_init = vacate_all = false;
 dag_vertex<properties>::reset();
 }*/
