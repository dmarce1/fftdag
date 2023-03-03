/*
 * fftdag.hpp
 *
 *  Created on: Mar 2, 2023
 *      Author: dmarce1
 */

#ifndef FFTDAG_HPP_
#define FFTDAG_HPP_

#include <cassert>
#include <memory>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <vector>
#include <deque>
#include <map>

enum vertex_type {
	INVALID, ADD, NEG, MUL, IN, CON, SUB, NSUB
};

inline bool close2(double a, double b) {
	return std::abs(a - b) < 1.0e-10;
}
;

inline bool is_arithmetic(vertex_type t) {
	switch (t) {
	case ADD:
	case SUB:
	case NSUB:
	case MUL:
	case NEG:
		return true;
	default:
		return false;
	}
}

class dag;

struct opcnt_t {
	int add;
	int mul;
	int neg;
	int tot;
};

class dag {
	struct dag_entry {
		vertex_type type;
		std::set<int> out;
		std::set<int> in;
		std::string name;
		int cnt;
		dag_entry() {
			cnt = 0;
			type = INVALID;
		}
	};
	int next_id;
	std::unordered_map<int, dag_entry> map;
	std::unordered_map<int, dag_entry> backup_map;
	std::set<int> inputs;
	std::set<int> outputs;
	void inc(int i) {
		map[i].cnt++;
	}
	void dec(int j) {
		map[j].cnt--;
		assert(map[j].cnt >= 0);
		if (map[j].cnt == 0) {
			auto in = map[j].in;
			auto out = map[j].out;
			for (auto i : in) {
				map[i].out.erase(j);
				dec(i);
			}
			inputs.erase(j);
			outputs.erase(j);
			map.erase(j);
		}
	}
public:
	~dag() {
		std::vector<int> vs;
		for (auto v : map) {
			vs.push_back(v.first);
		}
		for (auto v : vs) {
			auto eout = map[v].out;
			auto ein = map[v].in;
			for (auto e : eout) {
				map[v].out.erase(e);
				map[e].in.erase(v);
				dec(v);
			}
			for (auto e : ein) {
				map[e].out.erase(v);
				map[v].in.erase(e);
				dec(e);
			}
		}
	}
	class vertex {
		dag* graph_ptr;
		int id;
		friend class dag;
		operator int() const {
			return id;
		}
		void destroy() {
			if (graph_ptr) {
				assert(graph_ptr->map.find(id) != graph_ptr->map.end());
				graph_ptr->dec(id);
				graph_ptr = nullptr;
			}
		}
	public:
		vertex() {
			id = -1;
			graph_ptr = nullptr;
		}
		vertex(const vertex& other) {
			graph_ptr = nullptr;
			id = -1;
			*this = other;
		}
		vertex(vertex&& other) {
			graph_ptr = nullptr;
			id = -1;
			*this = std::move(other);
		}
		vertex& operator=(const vertex& other) {
			destroy();
			id = other.id;
			graph_ptr = other.graph_ptr;
			if (graph_ptr) {
				graph_ptr->inc(id);
			}
			return *this;
		}
		vertex& operator=(vertex&& other) {
			std::swap(id, other.id);
			std::swap(graph_ptr, other.graph_ptr);
			return *this;
		}
		~vertex() {
			destroy();
		}
		bool operator<(const vertex& other) const {
			return id < other.id;
		}
		bool operator==(const vertex& other) const {
			return id == other.id;
		}
	};
	friend class vertex;
	dag() {
		next_id = 0;
	}
	void backup() {
		backup_map = map;
	}
	void restore() {
		map = std::move(backup_map);
	}
	vertex make_vertex(vertex_type type) {
		vertex rc;
		rc.id = next_id++;
		rc.graph_ptr = this;
		map[rc.id].type = type;
		inc(rc.id);
		if (type == IN) {
			inputs.insert(rc.id);
		}
		return rc;
	}
	void set_output(vertex v) {
		assert(map.find(v) != map.end());
		outputs.insert(v);
	}
	void add_edge(vertex from, vertex to) {
		assert(map[from].out.find(to) == map[from].out.end());
		assert(map[to].in.find(from) == map[to].in.end());
		map[from].out.insert(to);
		map[to].in.insert(from);
		inc(from);
	}
	void remove_edge(vertex from, vertex to) {
		assert(map[from].out.find(to) != map[from].out.end());
		assert(map[to].in.find(from) != map[to].in.end());
		map[from].out.erase(to);
		map[to].in.erase(from);
		dec(from);
	}
	std::vector<vertex> get_edges_in(vertex v) {
		assert(map.find(v) != map.end());
		std::vector<vertex> edges;
		for (auto e : map[v].in) {
			vertex v;
			v.id = e;
			v.graph_ptr = this;
			inc(e);
			edges.push_back(v);
		}
		return edges;
	}
	std::vector<vertex> get_edges_out(vertex v) {
		assert(map.find(v) != map.end());
		std::vector<vertex> edges;
		for (auto e : map[v].out) {
			vertex v;
			v.id = e;
			v.graph_ptr = this;
			inc(e);
			edges.push_back(v);
		}
		return edges;
	}
	std::vector<vertex> get_inputs() {
		std::vector<vertex> ins;
		for (auto i : inputs) {
			vertex v;
			v.id = i;
			v.graph_ptr = this;
			inc(i);
			ins.push_back(v);
		}
		return ins;
	}
	std::set<vertex> get_outputs() {
		std::set<vertex> outs;
		for (auto o : outputs) {
			vertex v;
			v.id = o;
			v.graph_ptr = this;
			inc(o);
			outs.insert(v);
		}
		return outs;
	}
	void set_name(vertex v, std::string nm) {
		assert(map.find(v) != map.end());
		map[v].name = nm;
	}
	std::string get_name(vertex v) {
		return map[v].name;
	}
	std::vector<std::string> get_names(std::vector<vertex> vs) {
		std::vector<std::string> rc;
		for (auto v : vs) {
			rc.push_back(map[v].name);
		}
		return rc;
	}
	vertex_type get_type(vertex v) {
		auto entry = map[v];
		return entry.type;
	}
	bool sanity_check() {
		for (auto entry : map) {
			if (entry.second.type == ADD) {
				if (entry.second.in.size() != 2) {
					return false;
				}
			}
		}
		return true;
	}
};

class dag_node {
	using vertex = dag::vertex;
	vertex id;
//s	std::vector<vertex> children;
	static std::shared_ptr<dag> graph;
public:
	dag_node operator=(dag_node other) {
		id = other.id;
//		children = other.children;
		return *this;
	}

	friend dag_node binary_op(vertex_type type, dag_node a, dag_node b) {
		assert(sanity_check());
		dag_node c;
		c.id = graph->make_vertex(type);
		graph->add_edge(a.id, c.id);
		graph->add_edge(b.id, c.id);
		auto edges_in = graph->get_edges_in(c.id);
		//c.children = std::move(edges_in);
		assert(sanity_check());
		return c;
	}
	friend dag_node unary_op(vertex_type type, dag_node a) {
		dag_node c;
		c.id = graph->make_vertex(type);
		graph->add_edge(a.id, c.id);
//		c.children.clear();
//		c.children.push_back(a.id);
		return c;
	}

	friend dag_node operator+(dag_node a, dag_node b) {
		auto tmp = binary_op(ADD, a, b);
		assert(sanity_check());
		return tmp;
	}

	static bool sanity_check() {
		return graph->sanity_check();
	}

	friend dag_node operator*(dag_node a, dag_node b) {
		auto tmp = binary_op(MUL, a, b);
		assert(sanity_check());
		return tmp;
	}

	friend dag_node operator-(dag_node a) {
		auto tmp = unary_op(NEG, a);
		assert(sanity_check());
		return tmp;
	}

	friend dag_node operator-(dag_node a, dag_node b) {
		dag_node tmp;
		if (a.id < b.id) {
			tmp = binary_op(SUB, a, b);
		} else {
			tmp = binary_op(NSUB, std::move(b), std::move(a));
		}
		assert(sanity_check());
		return tmp;
	}
	dag_node operator+=(dag_node b) {
		auto c = *this + b;
		*this = c;
		assert(sanity_check());
		return *this;
	}
	dag_node operator-=(dag_node b) {
		auto c = *this - b;
		*this = c;
		assert(sanity_check());
		return *this;
	}
	dag_node operator*=(dag_node b) {
		auto c = *this * b;
		*this = c;
		assert(sanity_check());
		return *this;
	}
	dag_node() {
	}
	dag_node(double a) {
		id = graph->make_vertex(CON);
		graph->set_name(id, std::to_string(a));
	}
	static std::vector<dag_node> create_inputs(int N) {
		assert(sanity_check());
		std::vector<dag_node> rc;
		for (int n = 0; n < N; n++) {
			dag_node a;
			a.id = graph->make_vertex(IN);
			graph->set_name(a.id, std::string("x[") + std::to_string(n) + "]");
			rc.push_back(a);
		}
		assert(sanity_check());
		return rc;
	}
	static void set_outputs(std::vector<dag_node> outputs) {
		for (int i = 0; i < outputs.size(); i++) {
			graph->set_output(outputs[i].id);
			graph->set_name(outputs[i].id, std::string("x[") + std::to_string(i) + "]");
		}
	}

	static std::vector<vertex> search(vertex v, std::set<vertex>& touched) {
		std::vector<vertex> rc;
		if (touched.find(v) == touched.end()) {
			touched.insert(v);
			auto edges = graph->get_edges_in(v);
			for (auto c : edges) {
				auto this_search = search(c, touched);
				rc.insert(rc.end(), this_search.begin(), this_search.end());
			}
			rc.push_back(v);
		}
		return rc;
	}
	static std::vector<vertex> sort() {
		std::vector<vertex> L;
		std::set<vertex> touched;
		auto outputs = graph->get_outputs();
		for (auto o : outputs) {
			auto tmp = search(o, touched);
			L.insert(L.end(), tmp.begin(), tmp.end());
		}
		return L;

		/*
		 graph->sanity_check();
		 graph->backup();
		 auto inputs = graph->get_inputs();
		 auto S = std::set<vertex>(inputs.begin(), inputs.end());
		 while (S.size()) {
		 std::vector<vertex> C;
		 for (auto n : S) {
		 auto edges_in = graph->get_edges_in(n);
		 bool flag = true;
		 for (auto m : edges_in) {
		 if (graph->get_edges_in(m).size() > 0) {
		 flag = false;
		 break;
		 }
		 }
		 if (flag) {
		 C.push_back(n);
		 }
		 }
		 int best_cnt = -1;
		 int best_index;
		 for (int i = 0; i < C.size(); i++) {
		 auto edges_in = graph->get_edges_in(C[i]);
		 int cnt = 0;
		 for (auto m : edges_in) {
		 if (graph->get_edges_out(m).size() == 1) {
		 cnt++;
		 }
		 }
		 if (cnt > best_cnt) {
		 best_cnt = cnt;
		 best_index = i;
		 }
		 }
		 auto n = C[best_index];
		 auto edges_in = graph->get_edges_in(n);
		 for (auto m : edges_in) {
		 graph->remove_edge(m, n);
		 }
		 auto edges_out = graph->get_edges_out(n);
		 for (auto m : edges_out) {
		 S.insert(m);
		 }
		 S.erase(n);
		 L.push_back(n);
		 }
		 graph->restore();
		 graph->sanity_check();
		 return L;*/
	}
	static void print_list() {
		auto nodes = sort();
		int regcnt = 0;
		std::set<vertex> completed;
		std::set<std::string> free_vars;
		std::unordered_set<std::string> out_vars;
		std::unordered_map<std::string, vertex> used_vars;
		const auto outputs = graph->get_outputs();
		const auto inputs = graph->get_inputs();
		const auto genvar = [&regcnt,&free_vars, &out_vars]() {
			std::string nm;
			bool done = false;
			while(!done) {
				if (free_vars.size()) {
					nm = *free_vars.begin();
					free_vars.erase(nm);
				} else {
					nm = std::string("r") + std::to_string(regcnt++);
					printf("\tdouble %s;\n", nm.c_str());
				}
				done = (out_vars.find(nm) == out_vars.end());
			}
			return nm;
		};
		for (auto n : inputs) {
			used_vars.insert(std::make_pair(graph->get_name(n), n));
		}
		for (auto n : nodes) {
			auto type = graph->get_type(n);
			if (is_arithmetic(type)) {
				auto in = graph->get_edges_in(n);
				std::string nm;
				if (outputs.find(n) != outputs.end()) {
					nm = graph->get_name(n);
					if (free_vars.find(nm) == free_vars.end()) {
						assert(used_vars.find(nm) != used_vars.end());
						auto other = used_vars[nm];
						auto edges_out = graph->get_edges_out(other);
						bool flag = true;
						for (auto edge : edges_out) {
							if (completed.find(edge) == completed.end()) {
								flag = false;
								break;
							}
						}
						if (!flag) {
							auto other_nm = genvar();
							used_vars[other_nm] = other;
							used_vars.erase(nm);
							graph->set_name(other, other_nm);
							printf("\t%s = %s; // !!!!!\n", other_nm.c_str(), nm.c_str());
						}
					}
					out_vars.insert(nm);
				} else {
					nm = genvar();
					graph->set_name(n, nm);
				}
				used_vars[nm] = n;
				auto A = graph->get_name(n);
				auto names = graph->get_names(in);
				auto B = names.front();
				auto C = names.back();
				switch (type) {
				case ADD:
					assert(in.size() == 2);
					printf("\t%s = %s + %s;\n", A.c_str(), B.c_str(), C.c_str());
					break;
				case SUB:
					assert(in.size() == 2);
					printf("\t%s = %s - %s;\n", A.c_str(), B.c_str(), C.c_str());
					break;
				case NSUB:
					assert(in.size() == 2);
					printf("\t%s = %s - %s;\n", A.c_str(), C.c_str(), B.c_str());
					break;
				case MUL:
					assert(in.size() == 2);
					printf("\t%s = %s * %s;\n", A.c_str(), B.c_str(), C.c_str());
					break;
				case NEG:
					assert(in.size() == 1);
					printf("\t%s = -%s;\n", A.c_str(), B.c_str());
					break;
				default:
					assert(in.size() == 0);
					break;
				}
				completed.insert(n);
				auto edges_in = graph->get_edges_in(n);
				for (auto edge_in : edges_in) {
					if (is_arithmetic(graph->get_type(edge_in)) || graph->get_type(edge_in) == IN) {
						auto edges_out = graph->get_edges_out(edge_in);
						bool flag = true;
						for (auto edge : edges_out) {
							if (completed.find(edge) == completed.end()) {
								flag = false;
								break;
							}
						}
						if (flag) {
							auto nm = graph->get_name(edge_in);
							if (out_vars.find(nm) == out_vars.end()) {
								free_vars.insert(nm);
								used_vars.erase(nm);
							}
						}
					}
				}
			}
		}
	}
	static opcnt_t get_operation_count() {
		opcnt_t cnt;
		cnt.add = cnt.mul = cnt.neg = cnt.tot = 0;
		auto L = sort();
		for (auto n : L) {
			switch (graph->get_type(n)) {
			case ADD:
			case SUB:
			case NSUB:
				cnt.add++;
				cnt.tot++;
				break;
			case MUL:
				cnt.mul++;
				cnt.tot++;
				break;
			case NEG:
				cnt.neg++;
				cnt.tot++;
				break;
			default:
				break;
			}
		}
		return cnt;
	}

};

std::vector<dag_node> fft_prime_power(int R, std::vector<dag_node> xin, int N);

/*

 inline bool is_additive(vertex_type t) {
 switch (t) {
 case ADD:
 case SUB:
 case NSUB:
 return true;
 default:
 return false;
 }
 }

 inline bool is_binary(vertex_type t) {
 switch (t) {
 case ADD:
 case SUB:
 case MUL:
 case NSUB:
 return true;
 default:
 return false;
 }
 }


 class dag;
 class vertex_key;

 class vertex {
 std::shared_ptr<int> id;
 std::shared_ptr<dag> graph;
 vertex(int i) {
 id = std::make_shared<int>(i);
 }
 public:
 ~vertex();
 vertex() {
 }
 bool null() const {
 return id == nullptr;
 }
 friend class dag;
 friend class vertex_key;
 bool operator<(vertex other) const {
 if (id == nullptr && other.id == nullptr) {
 return true;
 } else if (id != nullptr && other.id != nullptr) {
 return *id < *other.id;
 }
 return false;
 }
 bool operator==(vertex other) const {
 if (id == nullptr && other.id == nullptr) {
 return true;
 } else if (id != nullptr && other.id != nullptr) {
 return *id == *other.id;
 }
 return false;
 }
 };

 struct vertex_key {
 size_t operator()(vertex v) const {
 std::hash<int> hash;
 return hash(*v.id);
 }
 };

 class dag {
 struct map_entry {
 std::set<vertex> in;
 std::set<vertex> out;
 };
 std::unordered_map<vertex, map_entry, vertex_key> map;
 std::unordered_map<vertex, vertex_type, vertex_key> type;
 std::unordered_map<vertex, std::string, vertex_key> name;
 std::unordered_map<vertex, double, vertex_key> value;
 std::set<vertex> outputs;
 std::set<vertex> inputs;
 int next_id;
 public:
 dag() {
 next_id = 0;
 }
 dag(const dag& other) = default;
 dag(dag&& other) = default;
 dag& operator=(const dag& other) = default;
 dag& operator=(dag&& other) = default;
 std::vector<vertex> get_edges_in(vertex v) {
 assert(map.find(v) != map.end());
 return std::vector<vertex>(map[v].in.begin(), map[v].in.end());
 }
 std::vector<vertex> get_edges_out(vertex v) {
 assert(map.find(v) != map.end());
 return std::vector<vertex>(map[v].out.begin(), map[v].out.end());
 }
 std::vector<vertex> get_edges(vertex v) {
 assert(map.find(v) != map.end());
 auto edges = std::vector<vertex>(map[v].out.begin(), map[v].out.end());
 edges.insert(edges.end(), map[v].in.begin(), map[v].in.end());
 return edges;
 }
 void set_value(vertex v, double val) {
 assert(map.find(v) != map.end());
 value[v] = val;
 }
 double get_value(vertex v) {
 assert(value.find(v) != value.end());
 return value[v];
 }
 void set_name(vertex v, const std::string& nm) {
 assert(map.find(v) != map.end());
 name[v] = nm;
 }
 std::vector<std::string> get_names(std::vector<vertex> indices) {
 std::vector<std::string> rc;
 for (auto i : indices) {
 rc.push_back(get_name(i));
 }
 return rc;
 }

 std::vector<vertex> search(vertex v, std::unordered_set<vertex, vertex_key>& touched) {
 assert(map.find(v) != map.end());
 std::vector<vertex> rc;
 if (touched.find(v) == touched.end()) {
 touched.insert(v);
 auto edges = get_edges_in(v);
 for (auto e : edges) {
 auto this_search = search(e, touched);
 rc.insert(rc.end(), this_search.begin(), this_search.end());
 }
 rc.push_back(v);
 }
 return rc;
 }
 vertex search_next(vertex v, std::unordered_set<vertex, vertex_key>& touched) {
 assert(map.find(v) != map.end());
 if (touched.find(v) == touched.end()) {
 auto edges = get_edges_in(v);
 for (auto e : edges) {
 auto this_rc = search_next(e, touched);
 if (!this_rc.null()) {
 return this_rc;
 }
 }
 touched.insert(v);
 return v;
 }
 return vertex();
 }
 void touch_all(vertex v, std::unordered_set<vertex, vertex_key>& touched) {
 auto edges_in = get_edges_in(v);
 touched.insert(v);
 for (auto edge : edges_in) {
 touch_all(edge, touched);
 }
 }
 std::vector<vertex> sort() {
 auto g = *this;
 std::vector<vertex> L;
 auto S = std::unordered_set<vertex, vertex_key>(g.inputs.begin(), g.inputs.end());
 while (S.size()) {
 std::vector<vertex> C;
 for (auto n : S) {
 auto edges_in = g.get_edges_in(n);
 bool flag = true;
 for (auto m : edges_in) {
 if (g.get_edges_in(m).size() > 0) {
 flag = false;
 break;
 }
 }
 if (flag) {
 C.push_back(n);
 }
 }
 int best_cnt = -1;
 int best_index;
 for (int i = 0; i < C.size(); i++) {
 auto edges_in = g.get_edges_in(C[i]);
 int cnt = 0;
 for (auto m : edges_in) {
 if (g.get_edges_out(m).size() == 1) {
 cnt++;
 }
 }
 if (cnt > best_cnt) {
 best_cnt = cnt;
 best_index = i;
 }
 }
 auto n = C[best_index];
 auto edges_in = g.get_edges_in(n);
 for (auto m : edges_in) {
 g.remove_edge(m, n);
 }
 auto edges_out = g.get_edges_out(n);
 for (auto m : edges_out) {
 S.insert(m);
 }
 S.erase(n);
 L.push_back(n);
 }
 return L;
 }
 vertex vertex_exists(vertex_type type, std::set<vertex> edges_in) {
 if (edges_in.size()) {
 std::set<vertex> to;
 for (auto i : edges_in) {
 if (!to.size()) {
 to = map[i].out;
 auto old = to;
 to.clear();
 for (auto test : old) {
 if (get_type(test) == type) {
 to.insert(test);
 }
 }
 } else {
 auto old = to;
 to.clear();
 for (auto test : map[i].out) {
 if (old.find(test) != old.end()) {
 to.insert(test);
 }
 }
 }
 if (!to.size()) {
 return vertex();
 }
 }
 return *(to.begin());
 }
 return vertex();
 }
 std::string get_name(vertex v) {
 assert(map.find(v) != map.end());
 return name[v];
 }
 vertex_type get_type(vertex v) {
 assert(map.find(v) != map.end());
 assert(type.find(v) != type.end());
 return type[v];
 }
 void add_output(vertex v) {
 assert(map.find(v) != map.end());
 outputs.insert(v);
 }
 std::set<vertex> get_outputs() {
 return outputs;
 }
 std::set<vertex> get_inputs() {
 return inputs;
 }
 vertex add_vertex(vertex_type t) {
 int id = next_id++;
 map[id] = map_entry();
 type[id] = t;
 if (t == IN) {
 inputs.insert(id);
 }
 return id;
 }
 void add_edge(vertex from, vertex to) {
 assert(map[from].out.find(to) == map[from].out.end());
 assert(map[to].in.find(to) == map[to].in.end());
 map[from].out.insert(to);
 map[to].in.insert(from);
 }
 void remove_vertex(vertex v) {
 assert(map.find(v) != map.end());
 assert(type.find(v) != type.end());
 if (type[v] == IN) {
 inputs.erase(v);
 }
 outputs.erase(v);
 for (auto i = map[v].in.begin(); i != map[v].in.end(); i++) {
 assert(map[*i].out.find(v) != map[*i].out.end());
 map[*i].out.erase(v);
 }
 for (auto i = map[v].out.begin(); i != map[v].out.end(); i++) {
 assert(map[*i].in.find(v) != map[*i].in.end());
 map[*i].in.erase(v);
 }
 map.erase(v);
 name.erase(v);
 type.erase(v);
 value.erase(v);
 }
 void remove_edge(vertex from, vertex to) {
 assert(map[from].out.find(to) != map[from].out.end());
 assert(map[to].in.find(from) != map[to].in.end());
 map[from].out.erase(to);
 map[to].in.erase(from);
 }
 };

 class dag_node {
 vertex id;
 static std::shared_ptr<dag> graph;
 static std::map<double, vertex> const_map;
 public:
 friend dag_node binary_op(vertex_type type, dag_node a, dag_node b) {
 dag_node c;
 vertex common = graph->vertex_exists(type, std::set<vertex>( { a.id, b.id }));
 if (!common.null()) {
 dag_node c;
 c.id = common;
 return c;
 } else {
 c.id = graph->add_vertex(type);
 graph->add_edge(a.id, c.id);
 graph->add_edge(b.id, c.id);
 }
 return c;
 }
 friend dag_node unary_op(vertex_type type, dag_node a) {
 dag_node c;
 vertex common = graph->vertex_exists(type, std::set<vertex>( { a.id }));
 if (!common.null()) {
 dag_node c;
 c.id = common;
 return c;
 } else {
 c.id = graph->add_vertex(type);
 graph->add_edge(a.id, c.id);
 }
 return c;
 }
 dag_node() {
 }
 dag_node(double a) {
 *this = a;
 }
 dag_node operator=(double a) {
 bool flag = false;
 auto iter = const_map.upper_bound(a);
 if (iter != const_map.end()) {
 if (close2(iter->first, a)) {
 id = iter->second;
 flag = true;
 } else if (iter != const_map.begin()) {
 iter--;
 if (close2(iter->first, a)) {
 id = iter->second;
 flag = true;
 }
 }
 }
 if (!flag) {
 id = graph->add_vertex(CON);
 graph->set_value(id, a);
 graph->set_name(id, std::to_string(a));
 const_map[a] = id;
 }
 return *this;
 }
 bool zero() {
 if (graph->get_type(id) == CON) {
 return close2(graph->get_value(id), 0.0);
 }
 return false;
 }
 bool one() {
 if (graph->get_type(id) == CON) {
 return close2(graph->get_value(id), 1.0);
 }
 return false;
 }
 bool none() {
 if (graph->get_type(id) == CON) {
 return close2(graph->get_value(id), -1.0);
 }
 return false;
 }
 static dag_node propagate_signs(dag_node n) {
 auto edges_in = graph->get_edges_in(n.id);
 dag_node next;
 for (auto m : edges_in) {
 next.id = m;
 propagate_signs(next);
 }
 auto type = graph->get_type(n.id);
 next.id = n.id;
 if (is_arithmetic(type)) {
 dag_node a, b;
 a.id = edges_in.front();
 b.id = edges_in.back();
 dag_node achild, bchild;
 bool aneg = graph->get_type(a.id) == NEG;
 bool bneg = graph->get_type(b.id) == NEG;
 if (aneg) {
 achild.id = graph->get_edges_in(a.id).front();
 }
 if (bneg) {
 bchild.id = graph->get_edges_in(b.id).front();
 }
 switch (type) {
 case ADD:
 if (aneg && bneg) {
 next = -(achild + bchild);
 } else if (aneg) {
 next = b - achild;
 } else if (bneg) {
 next = a - bchild;
 }
 break;
 case SUB:
 if (aneg && bneg) {
 next = bchild - achild;
 } else if (aneg) {
 next = -(achild + b);
 } else if (bneg) {
 next = a + bchild;
 }
 break;
 case NSUB:
 if (aneg && bneg) {
 next = achild - bchild;
 } else if (aneg) {
 next = achild + b;
 } else if (bneg) {
 next = -(a + bchild);
 }
 break;
 case MUL:
 if (aneg && bneg) {
 next = achild * bchild;
 } else if (aneg) {
 next = -(achild * b);
 } else if (bneg) {
 next = -(bchild * a);
 }
 break;
 case NEG:
 if (aneg) {
 next = achild;
 }
 break;
 }
 }
 return next;
 }
 friend dag_node operator-(dag_node a) {
 if (a.zero()) {
 return a;
 } else {
 return unary_op(NEG, a);
 }
 }
 friend dag_node operator+(dag_node a, dag_node b) {
 if (a.zero()) {
 return b;
 } else if (b.zero()) {
 return a;
 } else {
 return binary_op(ADD, a, b);
 }
 }
 friend dag_node operator-(dag_node a, dag_node b) {
 if (a.zero()) {
 return -b;
 } else if (b.zero()) {
 return a;
 } else if (a.id < b.id) {
 return binary_op(SUB, a, b);
 } else {
 return binary_op(NSUB, a, b);
 }
 }
 friend dag_node operator*(dag_node a, dag_node b) {
 if (a.zero() || b.zero()) {
 return a;
 } else if (a.one()) {
 return b;
 } else if (b.one()) {
 return a;
 } else if (a.none()) {
 return -b;
 } else if (b.none()) {
 return -a;
 } else {
 return binary_op(MUL, a, b);
 }
 }
 dag_node operator+=(dag_node b) {
 auto c = *this + b;
 *this = c;
 return c;
 }
 dag_node operator-=(dag_node b) {
 auto c = *this - b;
 *this = c;
 return c;
 }
 dag_node operator*=(dag_node b) {
 auto c = *this * b;
 *this = c;
 return c;
 }
 static dag_node create_input(int num) {
 dag_node in;
 in.id = graph->add_vertex(IN);
 graph->set_name(in.id, std::string("x[") + std::to_string(num) + "]");
 return in;
 }
 static void set_ouputs(std::vector<dag_node> outs) {
 for (int i = 0; i < outs.size(); i++) {
 graph->add_output(outs[i].id);
 graph->set_name(outs[i].id, std::string("x[") + std::to_string(i) + "]");
 }
 }
 static std::vector<dag_node> list() {
 std::vector<dag_node> nodes;
 auto indices = graph->sort();
 for (auto i : indices) {
 dag_node node;
 node.id = i;
 nodes.push_back(node);
 }
 return nodes;
 }
 static void print_list() {
 auto nodes = list();
 int regcnt = 0;
 std::unordered_set<vertex, vertex_key> completed;
 std::set<std::string> free_vars;
 std::unordered_set<std::string> out_vars;
 std::unordered_map<std::string, vertex> used_vars;
 const auto outputs = graph->get_outputs();
 const auto inputs = graph->get_inputs();
 const auto genvar = [&regcnt,&free_vars, &out_vars]() {
 std::string nm;
 bool done = false;
 while(!done) {
 if (free_vars.size()) {
 nm = *free_vars.begin();
 free_vars.erase(nm);
 } else {
 nm = std::string("r") + std::to_string(regcnt++);
 printf("\tdouble %s;\n", nm.c_str());
 }
 done = (out_vars.find(nm) == out_vars.end());
 }
 return nm;
 };
 for (auto n : inputs) {
 used_vars.insert(std::make_pair(graph->get_name(n), n));
 }
 for (auto n : nodes) {
 auto type = graph->get_type(n.id);
 if (is_arithmetic(type)) {
 auto in = graph->get_edges_in(n.id);
 std::string nm;
 if (outputs.find(n.id) != outputs.end()) {
 nm = graph->get_name(n.id);
 if (free_vars.find(nm) == free_vars.end()) {
 assert(used_vars.find(nm) != used_vars.end());
 auto other = used_vars[nm];
 auto edges_out = graph->get_edges_out(other);
 bool flag = true;
 for (auto edge : edges_out) {
 if (completed.find(edge) == completed.end()) {
 flag = false;
 break;
 }
 }
 if (!flag) {
 auto other_nm = genvar();
 used_vars[other_nm] = other;
 used_vars.erase(nm);
 graph->set_name(other, other_nm);
 printf("\t%s = %s; // !!!!!\n", other_nm.c_str(), nm.c_str());
 }
 }
 out_vars.insert(nm);
 } else {
 nm = genvar();
 graph->set_name(n.id, nm);
 }
 used_vars[nm] = n.id;
 auto A = graph->get_name(n.id);
 auto names = graph->get_names(in);
 auto B = names.front();
 auto C = names.back();
 switch (type) {
 case ADD:
 assert(in.size() == 2);
 printf("\t%s = %s + %s;\n", A.c_str(), B.c_str(), C.c_str());
 break;
 case SUB:
 assert(in.size() == 2);
 printf("\t%s = %s - %s;\n", A.c_str(), B.c_str(), C.c_str());
 break;
 case NSUB:
 assert(in.size() == 2);
 printf("\t%s = %s - %s;\n", A.c_str(), C.c_str(), B.c_str());
 break;
 case MUL:
 assert(in.size() == 2);
 printf("\t%s = %s * %s;\n", A.c_str(), B.c_str(), C.c_str());
 break;
 case NEG:
 assert(in.size() == 1);
 printf("\t%s = -%s;\n", A.c_str(), B.c_str());
 break;
 default:
 assert(in.size() == 0);
 break;
 }
 completed.insert(n.id);
 auto edges_in = graph->get_edges_in(n.id);
 for (auto edge_in : edges_in) {
 if (is_arithmetic(graph->get_type(edge_in)) || graph->get_type(edge_in) == IN) {
 auto edges_out = graph->get_edges_out(edge_in);
 bool flag = true;
 for (auto edge : edges_out) {
 if (completed.find(edge) == completed.end()) {
 flag = false;
 break;
 }
 }
 if (flag) {
 auto nm = graph->get_name(edge_in);
 if (out_vars.find(nm) == out_vars.end()) {
 free_vars.insert(nm);
 used_vars.erase(nm);
 }
 }
 }
 }
 }
 }
 }
 static opcnt_t get_operation_count() {
 opcnt_t cnt;
 cnt.add = cnt.mul = cnt.neg = cnt.tot = 0;
 auto L = list();
 for (auto n : L) {
 switch (graph->get_type(n.id)) {
 case ADD:
 case SUB:
 case NSUB:
 cnt.add++;
 cnt.tot++;
 break;
 case MUL:
 cnt.mul++;
 cnt.tot++;
 break;
 case NEG:
 cnt.neg++;
 cnt.tot++;
 break;
 default:
 break;
 }
 }
 return cnt;
 }
 static std::vector<dag_node> create_inputs(int sz) {
 std::vector<dag_node> in;
 for (int i = 0; i < sz; i++) {
 in.push_back(dag_node::create_input(i));
 }
 return in;
 }
 };
 */
#endif /* FFTDAG_HPP_ */
