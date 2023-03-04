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
#include <algorithm>
#include <numeric>
#include <deque>
#include <map>

enum vertex_type {
	INVALID, ADD, NEG, MUL, IN, CON, SUB, NSUB
};

inline bool close2(double a, double b) {
	return std::abs(a - b) < 1.0e-12;
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

inline bool is_subtraction(vertex_type t) {
	switch (t) {
	case SUB:
	case NSUB:
		return true;
	default:
		return false;
	}
}

inline bool is_antiops(vertex_type a, vertex_type b) {
	if (a == SUB && b == NSUB) {
		return true;
	}
	if (b == SUB && a == NSUB) {
		return true;
	}
	return false;
}

template<class Properties>
class dag;

struct opcnt_t {
	int add;
	int mul;
	int neg;
	int tot;
};

struct set_key {
	size_t operator()(const std::set<int>& s) const {
		int key = 1234;
		for (auto i : s) {
			key += 42;
			key *= (1234 + i);
		}
		return key;
	}
};

std::set<int> intersection(const std::set<int>& A, const std::set<int>& B);

struct math_props {
	vertex_type type;
	std::string name;
	double value;
	math_props() {
		type = INVALID;
	}
};

template<class Properties>
class dag {
	struct dag_entry {
		std::set<int> out;
		std::set<int> in;
		Properties props;
		int cnt;
		dag_entry() {
			cnt = 0;
		}
	};
	int next_id;
	std::unordered_map<int, dag_entry> map;
	std::unordered_map<int, std::set<int>> edge_map;
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
				edge_map[j].erase(i);
				map[i].out.erase(j);
				dec(i);
			}
			inputs.erase(j);
			outputs.erase(j);
			map.erase(j);
		}
	}
public:
	class vertex {
		dag* graph_ptr;
		int id;
		friend class dag;
		void destroy() {
			if (graph_ptr) {
				assert(graph_ptr->map.find(id) != graph_ptr->map.end());
				graph_ptr->dec(id);
				graph_ptr = nullptr;
			}
		}
	public:
		bool operator==(void* ptr) const {
			return (ptr == nullptr) && (id == -1);
		}
		bool operator!=(void* ptr) const {
			return !((ptr == nullptr) && (id == -1));
		}
		operator int() const {
			return id;
		}
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
		bool operator>(const vertex& other) const {
			return id > other.id;
		}
		bool operator==(const vertex& other) const {
			return id == other.id;
		}
	};
	friend class vertex;
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
	dag() {
		next_id = 1;
	}
	void backup() {
		backup_map = map;
	}
	void restore() {
		map = std::move(backup_map);
	}
	vertex make_vertex(bool in) {
		vertex rc;
		rc.id = next_id++;
		rc.graph_ptr = this;
		inc(rc.id);
		if (in) {
			inputs.insert(rc.id);
		}
		return rc;
	}
	void set_output(vertex v) {
		assert(map.find(v) != map.end());
		outputs.insert(v);
	}
	void clear_outputs() {
		outputs.clear();
	}
	void add_edge(vertex from, vertex to) {
		assert(map[from].out.find(to) == map[from].out.end());
		assert(map[to].in.find(from) == map[to].in.end());
		map[from].out.insert(to);
		edge_map[from].insert(to);
		map[to].in.insert(from);
		inc(from);
	}
	void remove_edge(vertex from, vertex to) {
		assert(map[from].out.find(to) != map[from].out.end());
		assert(map[to].in.find(from) != map[to].in.end());
		map[from].out.erase(to);
		edge_map[from].erase(to);
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
	std::vector<vertex> find_vertices(std::vector<vertex> ins) {
		std::vector<int> I(ins.begin(), ins.end());
		if (ins.size()) {
			auto x = edge_map[I[0]];
			for (int i = 1; i < ins.size(); i++) {
				x = intersection(x, edge_map[I[i]]);
			}
			std::vector<vertex> y;
			for( auto ele : x ) {
				vertex v;
				v.id = ele;
				v.graph_ptr = this;
				inc(ele);
				y.push_back(v);
			}
			return std::vector<vertex>(y.begin(), y.end());
		}
		return std::vector<vertex>();
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
	std::vector<vertex> search(vertex v, std::set<vertex>& touched) {
		std::vector<vertex> rc;
		if (touched.find(v) == touched.end()) {
			touched.insert(v);
			auto edges = get_edges_in(v);
			for (auto c : edges) {
				auto this_search = search(c, touched);
				rc.insert(rc.end(), this_search.begin(), this_search.end());
			}
			rc.push_back(v);
		}
		return rc;
	}

	std::vector<vertex> sort() {
		std::vector<vertex> L;
		std::set<vertex> touched;
		auto outputs = get_outputs();
		for (auto o : outputs) {
			auto tmp = search(o, touched);
			L.insert(L.end(), tmp.begin(), tmp.end());
		}
		return L;

	}
	Properties& props(vertex v) {
		return map[v].props;
	}
};



class math_dag : public dag<math_props> {
public:
	vertex make_vertex(vertex_type type) {
		vertex rc = dag<math_props>::make_vertex(type==IN);
		set_type(rc, type);
		return rc;
	}
	void set_name(vertex v, std::string nm) {
		props(v).name = nm;
	}
	std::string get_name(vertex v) {
		return props(v).name;
	}
	std::vector<std::string> get_names(std::vector<vertex> vs) {
		std::vector<std::string> rc;
		for (auto v : vs) {
			rc.push_back(props(v).name);
		}
		return rc;
	}
	vertex_type get_type(vertex v) {
		return props(v).type;
	}
	void set_value(vertex v, double val) {
		props(v).value = val;
	}
	void set_type(vertex v, vertex_type type) {
		props(v).type = type;
	}
	double get_value(vertex v) {
		return props(v).value;
	}
};



class dag_node {
	using vertex = dag<math_props>::vertex;
	vertex id;
	static std::shared_ptr<math_dag> graph;
	static std::map<double, vertex> const_map;
public:
	struct add_type {
		vertex v;
		double c;
	};

	static void reset() {
		const_map.clear();
		graph = std::make_shared<math_dag>();
	}
	dag_node(vertex v) {
		id = v;
	}

	bool operator<(const dag_node& other) {
		return id < other.id;
	}
	bool operator==(const dag_node& other) {
		return id == other.id;
	}
	static void optimize_adds() {
		auto all = collect_all_adds();
		for (auto a : all) {
			fprintf( stderr, "\n");
			for (auto b : a) {
				fprintf( stderr, "%i\n", (int) b);
			}
		}
	}
	static vertex find_vertex(vertex_type type, std::vector<vertex> ins) {
		auto vs = graph->find_vertices(ins);
		for (auto v : vs) {
			if (graph->get_type(v) == type) {
				return v;
			}
		}
		return vertex();
	}
	bool zero() {
		if (graph->get_type(id) == CON) {
			return close2(graph->get_value(id), 0.0);
		} else {
			return false;
		}
	}

	bool one() {
		if (graph->get_type(id) == CON) {
			return close2(graph->get_value(id), 1.0);
		} else {
			return false;
		}
	}

	bool none() {
		if (graph->get_type(id) == CON) {
			return close2(graph->get_value(id), -1.0);
		} else {
			return false;
		}
	}

	dag_node operator=(dag_node other) {
		id = other.id;
		return *this;
	}

	bool neg() const {
		return graph->get_type(id) == NEG;
	}
	dag_node neg_in() {
		dag_node in;
		in.id = graph->get_edges_in(id)[0];
		return in;
	}
	static std::vector<add_type> collect_mults(vertex v) {
		std::vector<add_type> rc;
		if (is_additive(graph->get_type(v))) {
			int i = 0;
			auto edges_in = graph->get_edges_in(v);
			for (auto e : edges_in) {
				auto tmp = collect_mults(e);
				if ((i == 0 && graph->get_type(v) == NSUB) || (i == 1 && graph->get_type(v) == SUB)) {
					for (auto& t : tmp) {
						t.c = -t.c;
					}
				}
				rc.insert(rc.end(), tmp.begin(), tmp.end());
				i++;
			}
		} else {
			add_type add;
			auto edges = graph->get_edges_in(v);
			if (graph->get_type(v) == MUL) {
				if (graph->get_type(edges[0]) == CON) {
					add.v = edges[1];
					add.c = graph->get_value(edges[0]);
				} else if (graph->get_type(edges[1]) == CON) {
					add.v = edges[0];
					add.c = graph->get_value(edges[1]);
				} else {
					add.v = v;
					add.c = 1;
				}
			} else if (graph->get_type(v) == IN) {
				add.c = 1.0;
				add.v = v;
			} else {
				add.c = -1.0;
				add.v = edges[0];
			}
			rc.push_back(add);
		}
		return rc;
	}
	static std::vector<int> collect_adds(vertex v) {
		std::vector<int> rc;
		if (is_additive(graph->get_type(v))) {
			int i = 0;
			auto edges_in = graph->get_edges_in(v);
			for (auto e : edges_in) {
				auto tmp = collect_adds(e);
				if ((i == 0 && graph->get_type(v) == NSUB) || (i == 1 && graph->get_type(v) == SUB)) {
					for (auto& t : tmp) {
						t = -t;
					}
				}
				rc.insert(rc.end(), tmp.begin(), tmp.end());
				i++;
			}
		} else {
			rc.push_back(v);
		}
		return rc;
	}
	static std::vector<std::vector<int>> collect_all_adds() {
		std::vector<std::vector<int>> rc;
		auto nodes = graph->sort();
		for (auto n : nodes) {
			auto type = graph->get_type(n);
			if (!is_additive(type)) {
				auto edges_out = graph->get_edges_out(n);
				for (auto edge : edges_out) {
					auto this_type = graph->get_type(edge);
					if (is_additive(this_type)) {
						rc.push_back(collect_adds(edge));
					}
				}
			}
		}
		return rc;
	}

	static dag_node optimize_mults(dag_node nd) {
		auto adds = collect_mults(nd.id);
		std::map<double, std::vector<add_type>> sorted_adds;
		for (int i = 0; i < adds.size(); i++) {
			auto j = sorted_adds.begin();
			while (j != sorted_adds.end()) {
				if (close2(j->first, std::abs(adds[i].c))) {
					break;
				}
				j++;
			}
			if (j == sorted_adds.end()) {
				sorted_adds[std::abs(adds[i].c)].clear();
				j = sorted_adds.find(std::abs(adds[i].c));
			}
			j->second.push_back(adds[i]);
		}
		dag_node main_sum(0.0);
		dag_node sum(0.0);
		for (auto j = sorted_adds.begin(); j != sorted_adds.end(); j++) {
			sum = dag_node(0.0);
			for (int i = 0; i < j->second.size(); i++) {
				auto v = j->second[i];
				assert(close2(std::abs(v.c), j->first));
				if (v.c > 0.0) {
					sum = sum + dag_node(v.v);
				} else {
					sum = sum - dag_node(v.v);
				}
			}
			auto coeff = j->first;
			if (close2(coeff, 1.0)) {
			} else if (close2(coeff, -1.0)) {
				sum = -sum;
			} else {
				sum = sum * dag_node(coeff);
			}
			main_sum = main_sum + sum;
		}
		sum = main_sum;
		auto ei = graph->get_edges_in(nd.id);
		for (auto e : ei) {
			graph->remove_edge(e, nd.id);
		}
		ei = graph->get_edges_in(sum.id);
		for (auto e : ei) {
			graph->remove_edge(e, sum.id);
			graph->add_edge(e, nd.id);
		}
		graph->set_type(nd.id, graph->get_type(sum.id));
		graph->set_value(nd.id, graph->get_value(sum.id));
		return nd;
	}

	static void optimize() {
		auto outputs = graph->get_outputs();
		graph->clear_outputs();
		for (auto o : outputs) {
			dag_node c;
			c.id = o;
			auto nm = graph->get_name(o);
			c = optimize_mults(c);
			c = optimize(c);
			graph->set_name(c.id, nm);
			graph->set_output(c.id);
		}
	}

	static dag_node optimize(dag_node node) {
		auto edges = graph->get_edges_in(node.id);
		for (auto e : edges) {
			dag_node other;
			other.id = e;
			other = optimize(other);
		}
		dag_node A;
		dag_node B;
		auto type = graph->get_type(node.id);
		if (is_arithmetic(type)) {
			A.id = edges.front();
			B.id = edges.back();
			switch (type) {
			case ADD:
				if (A.zero() && B.zero()) {
					node = dag_node(0.0);
				} else if (A.zero()) {
					node = B;
				} else if (B.zero()) {
					node = A;
				} else if (A.neg() && B.neg()) {
					node = -(A.neg_in() + B.neg_in());
				} else if (A.neg()) {
					node = B - A.neg_in();
				} else if (B.neg()) {
					node = A - B.neg_in();
				}
				break;
			case SUB:
				if (A.zero() && B.zero()) {
					node = dag_node(0.0);
				} else if (A.zero()) {
					node = -B;
				} else if (B.zero()) {
					node = A;
				} else if (A.neg() && B.neg()) {
					node = B.neg_in() - A.neg_in();
				} else if (A.neg()) {
					node = -(B + A.neg_in());
				} else if (B.neg()) {
					node = A + B.neg_in();
				}
				break;
			case NSUB:
				if (A.zero() && B.zero()) {
					node = dag_node(0.0);
				} else if (A.zero()) {
					node = B;
				} else if (B.zero()) {
					node = -A;
				} else if (A.neg() && B.neg()) {
					node = A.neg_in() - B.neg_in();
				} else if (B.neg()) {
					node = -(A + B.neg_in());
				} else if (A.neg()) {
					node = B + A.neg_in();
				}
				break;
			case MUL:
				if (A.zero() || B.zero()) {
					node = dag_node(0.0);
				} else if ((A.one() && B.one()) || (A.none() && B.none())) {
					node = dag_node(1.0);
				} else if ((A.one() && B.none()) || (A.none() && B.one())) {
					node = dag_node(-1.0);
				} else if (A.one()) {
					node = B;
				} else if (B.one()) {
					node = A;
				} else if (A.none()) {
					node = -B;
				} else if (B.none()) {
					node = -A;
				} else if (A.neg() && B.neg()) {
					node = B.neg_in() * A.neg_in();
				} else if (A.neg()) {
					node = -(B * A.neg_in());
				} else if (B.neg()) {
					node = -(B.neg_in() * A);
				}
				break;
			case NEG:
				if (A.neg()) {
					node = A.neg_in();
				}
				break;
			}
		}

		return node;
	}

	friend dag_node binary_op(vertex_type type, dag_node a, dag_node b) {
		auto common = find_vertex(type, std::vector<vertex>( { a.id, b.id }));
		dag_node c;
		if (common != nullptr) {
			c.id = common;
			return c;
		} else {
			if (type == SUB) {
				common = find_vertex(NSUB, std::vector<vertex>( { a.id, b.id }));
			} else if (type == NSUB) {
				common = find_vertex(SUB, std::vector<vertex>( { a.id, b.id }));
			} else {
				common = vertex();
			}
			if (common != nullptr) {
				c.id = common;
				return -c;
			}
		}
		if (a.id > b.id) {
			std::swap(a, b);
		}
		c.id = graph->make_vertex(type);
		graph->add_edge(a.id, c.id);
		graph->add_edge(b.id, c.id);
		auto edges_in = graph->get_edges_in(c.id);
		return optimize(c);
	}
	friend dag_node unary_op(vertex_type type, dag_node a) {
		dag_node c;
		auto common = find_vertex(type, std::vector<vertex>( { a.id }));
		if (common != nullptr) {
			c.id = common;
			return c;
		}
		c.id = graph->make_vertex(type);
		graph->add_edge(a.id, c.id);
		return optimize(c);
	}

	friend dag_node operator+(dag_node a, dag_node b) {
		if (a.id == b.id) {
			return dag_node(2.0) * a;
		} else {
			return binary_op(ADD, a, b);
		}
	}

	friend dag_node operator*(dag_node a, dag_node b) {
		if (a.id == b.id && graph->get_type(a.id) == CON) {
			auto v = graph->get_value(a.id);
			return dag_node(v * v);
		} else {
			return binary_op(MUL, a, b);
		}
	}

	friend dag_node operator-(dag_node a) {
		return unary_op(NEG, a);
	}

	friend dag_node operator-(dag_node a, dag_node b) {
		if (a.id < b.id) {
			return binary_op(SUB, a, b);
		} else if (a.id > b.id) {
			return binary_op(NSUB, b, a);
		} else {
			return dag_node(0.0);
		}
	}
	dag_node operator+=(dag_node b) {
		auto c = *this + b;
		*this = c;
		return *this;
	}
	dag_node operator-=(dag_node b) {
		auto c = *this - b;
		*this = c;
		return *this;
	}
	dag_node operator*=(dag_node b) {
		auto c = *this * b;
		*this = c;
		return *this;
	}
	dag_node() {
	}
	dag_node(double a) {
		if (a >= 0.0) {
			auto iter = const_map.upper_bound(a);
			if (iter != const_map.end()) {
				if (close2(a, iter->first)) {
					id = iter->second;
					return;
				} else if (iter != const_map.begin()) {
					iter--;
					if (close2(a, iter->first)) {
						id = iter->second;
						return;
					}
				}
			}
			id = graph->make_vertex(CON);
			graph->set_name(id, std::to_string(a));
			const_map[a] = id;
			graph->set_value(id, a);
		} else {
			*this = -dag_node(-a);
		}
	}
	static std::vector<dag_node> create_inputs(int N) {
		std::vector<dag_node> rc;
		for (int n = 0; n < N; n++) {
			dag_node a;
			a.id = graph->make_vertex(IN);
			graph->set_name(a.id, std::string("x[") + std::to_string(n) + "]");
			rc.push_back(a);
		}
		return rc;
	}
	static void set_outputs(std::vector<dag_node> outputs) {
		for (int i = 0; i < outputs.size(); i++) {
			graph->set_output(outputs[i].id);
			graph->set_name(outputs[i].id, std::string("x[") + std::to_string(i) + "]");
		}
	}

	static void print_list() {
		auto nodes = graph->sort();
		for (auto n : nodes) {
			switch (graph->get_type(n)) {
			case ADD:
				printf("ADD : ");
				break;
			case SUB:
				printf("SUB : ");
				break;
			case NSUB:
				printf("NSUB : ");
				break;
			case MUL:
				printf("MUL : ");
				break;
			case NEG:
				printf("NEG : ");
				break;
			case IN:
				printf("IN : ");
				break;
			case CON:
				printf("CON : ");
				break;
			}
			printf("%i | ", (int) n);
			auto edges = graph->get_edges_in(n);
			for (auto e : edges) {
				printf("%i ", (int) e);
			}
			printf("\n");
		}
	}
	static void print_code() {
		auto nodes = graph->sort();
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
							printf("\t%s = %s;\n", other_nm.c_str(), nm.c_str());
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
					printf("\t%s = -(%s);\n", A.c_str(), B.c_str());
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
		auto L = graph->sort();
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

}
;

std::vector<dag_node> fft_prime_power(int R, std::vector<dag_node> xin, int N);

#endif /* FFTDAG_HPP_ */
