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


enum vertex_type {
	ADD, NEG, MUL, IN, CON
};

inline bool is_arithmetic(vertex_type t) {
	switch (t) {
	case ADD:
	case MUL:
	case NEG:
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
	static int next_id;
public:
	std::vector<int> get_edges_in(int v) {
		assert(map.find(v) != map.end());
		return std::vector<int>(map[v].in.begin(), map[v].in.end());
	}
	std::vector<int> get_edges_out(int v) {
		assert(map.find(v) != map.end());
		return std::vector<int>(map[v].out.begin(), map[v].out.end());
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
	std::vector<int> sort() {
		std::vector<int> rc;
		std::unordered_set<int> touched;
		for (auto v : outputs) {
			auto results = search(v, touched);
			rc.insert(rc.end(), results.begin(), results.end());
		}
		return rc;
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
		assert(map[to].in.find(to) != map[to].in.end());
		map[from].out.erase(to);
		map[to].in.erase(from);
	}
};


class dag_node {
	int id;
	static std::shared_ptr<dag> graph;
public:
	friend dag_node binary_op(vertex_type type, dag_node a, dag_node b) {
		dag_node c;
		c.id = graph->add_vertex(type);
		graph->add_edge(a.id, c.id);
		graph->add_edge(b.id, c.id);
		return c;
	}
	friend dag_node unary_op(vertex_type type, dag_node a) {
		dag_node c;
		c.id = graph->add_vertex(type);
		graph->add_edge(a.id, c.id);
		return c;
	}
	dag_node() {
	}
	dag_node(double a) {
		*this = a;
	}
	dag_node operator=(double a) {
		id = graph->add_vertex(CON);
		graph->set_value(id, a);
		graph->set_name(id, std::to_string(a));
		return *this;
	}
	friend dag_node operator-(dag_node a) {
		return unary_op(NEG, a);
	}
	friend dag_node operator+(dag_node a, dag_node b) {
		return binary_op(ADD, a, b);
	}
	friend dag_node operator-(dag_node a, dag_node b) {
		return binary_op(ADD, a, unary_op(NEG, b));
	}
	friend dag_node operator*(dag_node a, dag_node b) {
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
		graph->set_name(in.id, std::string("xin[") + std::to_string(num) + "]");
		return in;
	}
	static void set_ouputs(std::vector<dag_node> outs) {
		for (int i = 0; i < outs.size(); i++) {
			graph->add_output(outs[i].id);
			graph->set_name(outs[i].id, std::string("xout[") + std::to_string(i) + "]");
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
		for (auto n : nodes) {
			auto type = graph->get_type(n.id);
			if (is_arithmetic(type)) {
				auto in = graph->get_edges_in(n.id);
				if (graph->get_name(n.id) == "") {
					auto regnum = regcnt++;
					printf( "\tdouble r%i;\n", regnum);
					graph->set_name(n.id, std::string("r") + std::to_string(regnum));
				}
				auto A = graph->get_name(n.id);
				auto names = graph->get_names(in);
				auto B = names.front();
				auto C = names.back();
				switch (type) {
				case ADD:
					assert(in.size() == 2);
					printf("\t%s = %s + %s;\n", A.c_str(), B.c_str(), C.c_str());
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
