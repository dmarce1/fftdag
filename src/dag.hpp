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
#include <vector>

template<class Properties>
class dag_vertex {
public:
	using function_type = std::function<void(Properties&, const std::vector<Properties>&)>;
private:
	struct state {
		int id;
		Properties props;
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
	using executor = std::unordered_map<std::weak_ptr<state>, bool, key, equal>;
private:
	void sort(executor& touched, std::vector<dag_vertex>& vertices) const {
		const auto& edges_in = state_ptr->edges_in;
		for (const auto& e : edges_in) {
			e.sort(touched, vertices);
		}
		if (!touched[state_ptr]) {
			vertices.push_back(*this);
			touched[state_ptr] = true;
		}
	}
public:
	class weak_ref {
		std::weak_ptr<state> ptr;
	public:
		weak_ref() = default;
		weak_ref(const dag_vertex& v) {
			ptr = v.state_ptr;
		}
		friend class dag_vertex;
	};
	dag_vertex(const weak_ref& ref) {
		state_ptr = std::shared_ptr<state>(ref.ptr);
	}
	dag_vertex() = default;
	static dag_vertex new_(Properties&& props) {
		dag_vertex v;
		auto sptr = new state;
		v.state_ptr = std::shared_ptr<state>(sptr);
		v.state_ptr->props = std::move(props);
		v.state_ptr->id = next_id;
		next_id++;
		return std::move(v);
	}
	void execute(executor& touched, const function_type& func) {
		if (!touched[state_ptr]) {
			std::vector<Properties> props;
			auto& edges_in = state_ptr->edges_in;
			for (auto& e : edges_in) {
				e.execute(touched, func);
			}
			for (const auto& e : edges_in) {
				props.push_back(e.properties());
			}
			func(state_ptr->props, std::move(props));
			state_ptr->edges_in = std::vector<dag_vertex>();
			touched[state_ptr] = true;
		}
	}
	std::vector<dag_vertex> sort() {
		std::vector<dag_vertex> vertices;
		executor touched;
		sort(vertices, touched);
		return std::move(vertices);
	}
	void execute(const function_type& func) {
		executor touched;
		execute(touched);
	}
	int get_unique_id() const {
		return state_ptr->id;
	}
	bool operator==(const dag_vertex& v) const {
		return state_ptr->id == v.state_ptr->id;
	}
	bool operator!=(const dag_vertex& v) const {
		return state_ptr->id != v.state_ptr->id;
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
	void add_edge_in(const dag_vertex& v) {
		state_ptr->edges_in.push_back(v);
	}
	void remove_edge_in(const dag_vertex& v) {
		auto& edges_in = state_ptr->edges_in;
		for (int i = 0; i < edges_in.size(); i++) {
			if (edges_in[i] == v) {
				for (int j = i; j < edges_in.size() - 1; j++) {
					edges_in[j] = edges_in[j + 1];
				}
				edges_in.pop_back();
				return;
			}
			assert(false);
		}
	}
	dag_vertex get_edge_in(int i) const {
		return state_ptr->edges_in[i];
	}
	int get_edge_count() const {
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
	int use_count() const {
		return state_ptr.use_count();
	}
};

template<class Properties>
int dag_vertex<Properties>::next_id = 1;

#endif /* FFTDAG_HPP_ */
