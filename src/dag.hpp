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
	ADD, NEG, MUL, IN, CON, SUB, NSUB
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

class dag {
	struct map_entry {
		std::set<int> in;
		std::set<int> out;
	};
	std::unordered_map<int, map_entry> map;
	std::unordered_map<int, vertex_type> type;
	std::unordered_map<int, std::string> name;
	std::unordered_map<int, double> value;
	std::set<int> outputs;
	std::set<int> inputs;
	int next_id;
public:
	dag() {
		next_id = 0;
	}
	dag(const dag& other) = default;
	dag(dag&& other) = default;
	dag& operator=(const dag& other) = default;
	dag& operator=(dag&& other) = default;
	std::vector<int> get_edges_in(int v) {
		assert(map.find(v) != map.end());
		return std::vector<int>(map[v].in.begin(), map[v].in.end());
	}
	std::vector<int> get_edges_out(int v) {
		assert(map.find(v) != map.end());
		return std::vector<int>(map[v].out.begin(), map[v].out.end());
	}
	std::vector<int> get_edges(int v) {
		assert(map.find(v) != map.end());
		auto edges = std::vector<int>(map[v].out.begin(), map[v].out.end());
		edges.insert(edges.end(), map[v].in.begin(), map[v].in.end());
		return edges;
	}
	void set_value(int v, double val) {
		assert(map.find(v) != map.end());
		value[v] = val;
	}
	double get_value(int v) {
		assert(value.find(v) != value.end());
		return value[v];
	}
	void set_name(int v, const std::string& nm) {
		assert(map.find(v) != map.end());
		name[v] = nm;
	}
	std::vector<std::string> get_names(std::vector<int> indices) {
		std::vector<std::string> rc;
		for (auto i : indices) {
			rc.push_back(get_name(i));
		}
		return rc;
	}

	std::vector<int> search(int v, std::unordered_set<int>& touched) {
		assert(map.find(v) != map.end());
		std::vector<int> rc;
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
	int search_next(int v, std::unordered_set<int>& touched) {
		assert(map.find(v) != map.end());
		std::vector<int> rc;
		if (touched.find(v) == touched.end()) {
			auto edges = get_edges_in(v);
			for (auto e : edges) {
				auto this_rc = search_next(e, touched);
				if (this_rc != -1) {
					return this_rc;
				}
			}
			touched.insert(v);
			return v;
		}
		return -1;
	}
	std::vector<int> sort() {
		std::vector<int> rc;
		std::unordered_set<int> touched;
		for (auto v : outputs) {
			auto results = search(v, touched);
			rc.insert(rc.end(), results.begin(), results.end());
		}
		return rc;
	}
	std::vector<int> sort3() {
		auto g = *this;
		std::vector<int> L;
		auto S = std::unordered_set<int>(g.inputs.begin(), g.inputs.end());
		while (S.size()) {
			std::vector<int> C;
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
	std::vector<int> sort2() {
		auto oldgraph = *this;
		std::vector<int> L;
		auto S = std::deque<int>(inputs.begin(), inputs.end());
		int sw = 1;
		while (S.size()) {
			auto n = S.back();
			S.pop_back();
			L.push_back(n);
			auto eout = get_edges_out(n);
			for (auto m : eout) {
				remove_edge(n, m);
				auto edges_in = get_edges_in(m);
				bool flag = false;
				for (auto in : edges_in) {
					if (type[in] != CON) {
						flag = true;
						break;
					}
				}
				if (!flag) {
					if (sw) {
						S.push_back(m);
						sw = 0;
					} else {
						S.push_front(m);
						sw = 1;
					}
				}
			}
		}
		*this = oldgraph;
		return L;
	}
	int vertex_exists(vertex_type type, std::set<int> edges_in) {
		if (edges_in.size()) {
			std::set<int> to;
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
					return -1;
				}
			}
			return *(to.begin());
		}
		return -1;
	}
	std::string get_name(int v) {
		assert(map.find(v) != map.end());
		return name[v];
	}
	vertex_type get_type(int v) {
		assert(map.find(v) != map.end());
		assert(type.find(v) != type.end());
		return type[v];
	}
	void add_output(int v) {
		assert(map.find(v) != map.end());
		outputs.insert(v);
	}
	std::set<int> get_outputs() {
		return outputs;
	}
	std::set<int> get_inputs() {
		return inputs;
	}
	int add_vertex(vertex_type t) {
		int id = next_id++;
		map[id] = map_entry();
		type[id] = t;
		if (t == IN) {
			inputs.insert(id);
		}
		return id;
	}
	void add_edge(int from, int to) {
		assert(map[from].out.find(to) == map[from].out.end());
		assert(map[to].in.find(to) == map[to].in.end());
		map[from].out.insert(to);
		map[to].in.insert(from);
	}
	void remove_vertex(int v) {
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
	void remove_edge(int from, int to) {
		assert(map[from].out.find(to) != map[from].out.end());
		assert(map[to].in.find(from) != map[to].in.end());
		map[from].out.erase(to);
		map[to].in.erase(from);
	}
};

class dag_node {
	int id;
	static std::shared_ptr<dag> graph;
	static std::map<double, int> const_map;
public:
	friend dag_node binary_op(vertex_type type, dag_node a, dag_node b) {
		dag_node c;
		int common = graph->vertex_exists(type, std::set<int>( { a.id, b.id }));
		if (common != -1) {
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
		int common = graph->vertex_exists(type, std::set<int>( { a.id }));
		if (common != -1) {
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
			} else if (const_map.size() > 1) {
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
	friend dag_node operator-(dag_node a) {
		if (a.zero()) {
			return a;
		}
		return unary_op(NEG, a);
	}
	friend dag_node operator+(dag_node a, dag_node b) {
		if (a.zero()) {
			return b;
		} else if (b.zero()) {
			return a;
		}
		return binary_op(ADD, a, b);
	}
	friend dag_node operator-(dag_node a, dag_node b) {
		if (a.zero()) {
			return -b;
		} else if (b.zero()) {
			return a;
		}
		if (a.id < b.id) {
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
		}
		return binary_op(MUL, a, b);
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
	static std::vector<dag_node> list(dag& g) {
		std::vector<dag_node> nodes;
		auto indices = g.sort3();
		for (auto i : indices) {
			dag_node node;
			node.id = i;
			nodes.push_back(node);
		}
		return nodes;
	}
	static void print_list() {
		auto nodes = list(*graph);
		int regcnt = 0;
		std::unordered_set<int> completed;
		std::set<std::string> free_vars;
		std::unordered_set<std::string> out_vars;
		std::unordered_map<std::string, int> used_vars;
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
	static std::vector<dag_node> create_inputs(int sz) {
		std::vector<dag_node> in;
		for (int i = 0; i < sz; i++) {
			in.push_back(dag_node::create_input(i));
		}
		return in;
	}
};

#endif /* FFTDAG_HPP_ */
