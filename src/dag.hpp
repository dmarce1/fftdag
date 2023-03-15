/*
 * fftdag.hpp
 *
 *  Created on: Mar 2, 2023
 *      Author: dmarce1
 */

#ifndef FFTDAG_HPP_
#define FFTDAG_HPP_

#include <cassert>
#include <functional>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template<class Properties>
class dag_vertex {
public:
	using function_type = std::function<void(dag_vertex&, const std::vector<dag_vertex>&)>;
private:
	struct state {
		int id;
		Properties props;
		std::vector<dag_vertex> edges_in;
//		std::vector<dag_vertex> edges_out;
	};
	std::shared_ptr<state> state_ptr;
	static int next_id;
public:
	struct key {
		size_t operator()(const std::weak_ptr<state>& s) const {
			std::hash<state*> hash;
			return hash(std::shared_ptr<state>(s).get());
		}
	};
	struct equal {
		size_t operator()(const std::weak_ptr<state>& a, const std::weak_ptr<state>& b) const {
			return std::shared_ptr<state>(a).get() == std::shared_ptr<state>(b).get();
		}
	};
	struct executor {
		std::unordered_set<int> touched;
		bool free;
		executor(bool f = true) {
			free = f;
		}
	};
public:
	class weak_ref {
		std::weak_ptr<state> ptr;
		int id;
	public:
		weak_ref() = default;
		weak_ref(const dag_vertex& v) {
			ptr = v.state_ptr;
			id = v.state_ptr->id;
		}
		bool operator<(const weak_ref& other) const {
			return id < other.id;
		}
		bool operator==(const weak_ref& other) const {
			return id == other.id;
		}
		friend class dag_vertex;
	};
	dag_vertex(const weak_ref& ref) {
		state_ptr = std::shared_ptr<state>(ref.ptr);
	}
	dag_vertex() = default;
	~dag_vertex() {
		if (state_ptr) {
			if (state_ptr.use_count() == 1) {
				free_edges();
				directory.erase(state_ptr->id);
			}
		}
	}
	bool valid() const {
		return state_ptr != nullptr;
	}
	static dag_vertex generate_from_id(int id) {
		assert(directory[id].ptr.use_count());
		return directory[id];
	}
	static dag_vertex new_(Properties&& props) {
		dag_vertex v;
		auto sptr = new state;
		v.state_ptr = std::shared_ptr<state>(sptr);
		v.state_ptr->props = std::move(props);
		v.state_ptr->id = next_id;
		directory[next_id] = v;
		next_id++;
		return std::move(v);
	}
	int use_count() const {
		return state_ptr.use_count();
	}
	void free_edges() {
		auto& edges_in = state_ptr->edges_in;
		while (edges_in.size()) {
		//	edges_in.back().remove_edge_out(*this);
			edges_in.pop_back();
		}
	}
	void execute(executor& exe, const function_type& func) {
		if (exe.touched.find(state_ptr->id) == exe.touched.end()) {
			auto& edges_in = state_ptr->edges_in;
			for (auto& e : edges_in) {
				e.execute(exe, func);
			}
			func(*this, edges_in);
			if (exe.free) {
				free_edges();
			}
			exe.touched.insert(state_ptr->id);
		}
	}
	std::vector<dag_vertex> sort(executor& exe) {
		std::vector<dag_vertex> list;
		if (exe.touched.find(state_ptr->id) == exe.touched.end()) {
			std::vector<Properties> props;
			auto& edges_in = state_ptr->edges_in;
			for (auto& e : edges_in) {
				auto tmp = e.sort(exe);
				list.insert(list.end(), tmp.begin(), tmp.end());
			}
			list.push_back(*this);
			exe.touched.insert(state_ptr->id);
		}
		return std::move(list);
	}
	static std::vector<dag_vertex> sort(executor& exe, std::vector<dag_vertex>& outputs) {
		std::vector<dag_vertex> list;
		executor touched;
		for (auto o : outputs) {
			auto tmp = o.sort(exe);
			list.insert(list.end(), tmp.begin(), tmp.end());
		}
		return std::move(list);
	}
	void execute(const function_type& func) {
		executor touched;
		execute(touched);
	}
	int get_unique_id() const {
		return state_ptr->id;
	}
	bool operator==(const dag_vertex& v) const {
		if (!state_ptr && !v.state_ptr) {
			return true;
		} else if (state_ptr && v.state_ptr) {
			return state_ptr->id == v.state_ptr->id;
		} else {
			return false;
		}
	}
	bool operator!=(const dag_vertex& v) const {
		return !operator==(v);
	}
	bool operator<(const dag_vertex& v) const {
		return state_ptr->id < v.state_ptr->id;
	}
	bool operator==(void* ptr) {
		return ptr == nullptr && state_ptr == nullptr;
	}
	bool operator!=(void* ptr) {
		return ptr == nullptr && state_ptr != nullptr;
	}
//	void add_edge_out(const dag_vertex& v) {
//		state_ptr->edges_out.push_back(v);
//	}
	void add_edge_in(dag_vertex& v) {
		state_ptr->edges_in.push_back(v);
//		v.add_edge_out(*this);
	}
/*	void remove_edge_out(const dag_vertex& v) {
		auto& edges_out = state_ptr->edges_out;
		for (int i = 0; i < edges_out.size(); i++) {
			if (edges_out[i] == v) {
				edges_out[i] = edges_out.back();
				edges_out.pop_back();
				return;
			}
		}
		assert(false);
	}*/
	void remove_edge_in(const dag_vertex& v) {
		auto& edges_in = state_ptr->edges_in;
		for (int i = 0; i < edges_in.size(); i++) {
			if (edges_in[i] == v) {
				//edges_in[i].remove_edge_out(*this);
				//for (int j = i; j < edges_in.size() - 1; j++) {
				//	edges_in[j] = edges_in[j + 1];
			//	}
				edges_in.pop_back();
				return;
			}
			assert(false);
		}
	}
	dag_vertex get_edge_in(int i) const {
		return state_ptr->edges_in[i];
	}
	int get_edge_in_count() const {
		return state_ptr->edges_in.size();
	}
//	dag_vertex get_edge_out(int i) const {
//		return state_ptr->edges_out[i];
//	}
//	int get_edge_out_count() const {
//		return state_ptr->edges_out.size();
//	}
	void replace_edge_in(const dag_vertex& v, dag_vertex&& u) {
		for (auto& edge : state_ptr->edges_in) {
			if (edge == v) {
//				edge.remove_edge_out(*this);
				edge = std::move(u);
//				edge.add_edge_out(*this);
				return;
			}
		}
	}
	Properties& properties() {
		return state_ptr->props;
	}
	const Properties& properties() const {
		return state_ptr->props;
	}
	static void reset() {
		next_id = 1;
		directory.clear();
	}
private:
	static std::unordered_map<int, dag_vertex<Properties>::weak_ref> directory;
};

template<class Properties>
int dag_vertex<Properties>::next_id = 1;

template<class Properties>
std::unordered_map<int, typename dag_vertex<Properties>::weak_ref> dag_vertex<Properties>::directory;

#endif /* FFTDAG_HPP_ */
