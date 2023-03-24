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
#include <algorithm>
#include <vector>
#include <cmath>

template<class Properties>
class dag_vertex {
public:
	using function_type = std::function<void(dag_vertex&, const std::vector<dag_vertex>&)>;
private:
	struct state {
		int id;
		Properties props;
		int score;
		std::vector<dag_vertex> edges_in;
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
	int compute_labels(executor& exe, bool left_desc = false) {
		if (exe.touched.find(state_ptr->id) == exe.touched.end()) {
			auto& edges_in = state_ptr->edges_in;
			if (edges_in.size() == 2) {
				std::vector<int> cscores;
				left_desc = true;
				for (auto& c : edges_in) {
					if (exe.touched.find(c.state_ptr->id) == exe.touched.end()) {
						int ci = c.compute_labels(exe, left_desc);
						cscores.push_back(ci);
						left_desc = false;
					}
				}
				if (cscores.size() == 2) {
					if (cscores[0] != cscores[1]) {
						state_ptr->score = std::max(cscores[0], cscores[1]);
					} else {
						state_ptr->score = cscores[0] + 1;
					}
				} else if (cscores.size()) {
					state_ptr->score = cscores[0];
				} else {
					state_ptr->score = 0;
				}
			} else if (edges_in.size() == 1) {
				state_ptr->score = get_edge_in(0).compute_labels(exe, true);
			} else {
				state_ptr->score = left_desc ? 1 : 0;
			}
			exe.touched.insert(state_ptr->id);
		}
		return state_ptr->score;
	}
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
		int use_count() const {
			return ptr.use_count();
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
			}
		}
	}
	bool valid() const {
		return state_ptr != nullptr;
	}
	static dag_vertex new_(Properties&& props) {
		dag_vertex v;
		auto sptr = new state;
		v.state_ptr = std::shared_ptr<state>(sptr);
		v.state_ptr->props = std::move(props);
		v.state_ptr->id = next_id;
		next_id++;
		return std::move(v);
	}
	int use_count() const {
		return state_ptr.use_count();
	}
	void free_edges() {
		auto& edges_in = state_ptr->edges_in;
		while (edges_in.size()) {
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
	std::vector<dag_vertex> pebble_game_sort(executor& exe) {
		std::vector<dag_vertex> list;
		if (exe.touched.find(state_ptr->id) == exe.touched.end()) {
			auto edges_in = state_ptr->edges_in;
			std::sort(edges_in.begin(), edges_in.end(), [](dag_vertex a, dag_vertex b) {
				return a.state_ptr->score > b.state_ptr->score;
			});
			if (edges_in.size() == 2) {
				if (edges_in[1].state_ptr->score >= edges_in[0].state_ptr->score) {
					std::swap(edges_in[1], edges_in[0]);
				}
			}
			for (auto& e : edges_in) {
				auto tmp = e.pebble_game_sort(exe);
				list.insert(list.end(), tmp.begin(), tmp.end());
			}
			list.push_back(*this);
			exe.touched.insert(state_ptr->id);
		}
		return std::move(list);
	}
	static std::vector<dag_vertex> pebble_game_sort(executor& exe, std::vector<dag_vertex>& outputs) {
		for (auto o : outputs) {
			o.compute_labels(exe);
		}
		exe = executor();
		std::vector<dag_vertex> list;
		for (auto o : outputs) {
			auto tmp = o.pebble_game_sort(exe);
			list.insert(list.end(), tmp.begin(), tmp.end());
		}
		return std::move(list);
	}
	void execute(const function_type& func) {
		executor touched;
		execute(touched);
	}
	int get_unique_id() const {
		if (state_ptr) {
			return state_ptr->id;
		} else {
			return -1;
		}
	}
	dag_vertex(const dag_vertex& other) = default;
	dag_vertex& operator=(const dag_vertex& other) {
		if (get_unique_id() != other.get_unique_id()) {
			state_ptr = other.state_ptr;
		}
		return *this;
	}
	dag_vertex& operator=(dag_vertex&& other) {
		if (get_unique_id() != other.get_unique_id()) {
			state_ptr = std::move(other.state_ptr);
		}
		return *this;
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
	void add_edge_in(dag_vertex& v) {
		state_ptr->edges_in.push_back(v);
	}
	void remove_edge_in(const dag_vertex& v) {
		auto& edges_in = state_ptr->edges_in;
		for (int i = 0; i < edges_in.size(); i++) {
			if (edges_in[i] == v) {
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
	void replace_edge_in(const dag_vertex& v, dag_vertex&& u) {
		for (auto& edge : state_ptr->edges_in) {
			if (edge == v) {
				edge = std::move(u);
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
	}
};

template<class Properties>
int dag_vertex<Properties>::next_id = 1;

#endif /* FFTDAG_HPP_ */
