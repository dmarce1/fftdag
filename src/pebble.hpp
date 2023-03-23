/*
 * pebble.hpp
 *
 *  Created on: Mar 19, 2023
 *      Author: dmarce1
 */

#ifndef PEBBLE_HPP_
#define PEBBLE_HPP_

#include <stack>
#include <unordered_set>
#include <vector>
#include "network.hpp"

class pebble_game {
	struct pebbles_type {
		std::vector<bool> pebbled;
		std::unordered_set<int> pebbles;
	};

	int N;
	std::vector<std::vector<int>> input_map;
	std::vector<std::vector<int>> output_map;
	std::vector<int> stack;
	std::unordered_set<int> inputs;
	std::unordered_set<int> outputs;
	std::vector<bool> retired;
	std::vector<bool> exists;
	pebbles_type board;
public:
	void swap_pebble(int sq) {
		if (board.pebbled[sq]) {
			board.pebbled[sq] = false;
			board.pebbles.erase(sq);
		} else {
			board.pebbled[sq] = true;
			board.pebbles.insert(sq);
		}
	}
	bool make_move(int mv) {
		assert(legal_move(mv));
		if (!board.pebbled[mv]) {
			retired[mv] = true;
		}
		swap_pebble(mv);
		stack.push_back(mv);
		return true;
	}
	void undo_move() {
		auto mv = stack.back();
		stack.pop_back();
		if (board.pebbled[mv]) {
			retired[mv] = false;
		}
		swap_pebble(mv);
	}
	void prepare() {
		for (int i = 0; i < N; i++) {
			if (input_map[i].size() == 0 && exists[i]) {
				inputs.insert(i);
				make_move(i);
			}
		}
		for (int i = 0; i < N; i++) {
			if (output_map[i].size() == 0 && exists[i]) {
				outputs.insert(i);
			}
		}
		stack.clear();
	}
	bool complete() {
		int cnt = 0;
		for (auto oi : outputs) {
			if (!board.pebbled[oi] && exists[oi]) {
				return false;
			}
		}
		return true;
	}
	int completeness() {
		int cnt = 0;
		for (auto oi : outputs) {
			if (board.pebbled[oi] && exists[oi]) {
				cnt++;
			}
		}
		return cnt;
	}
	bool legal_move(int sq) {
		if (board.pebbled[sq]) {
			bool f = output_map[sq].size();
			for (int i : output_map[sq]) {
				if (!board.pebbled[i]) {
					f = false;
					break;
				}
			}
			return f;
		} else {
			if (retired[sq]) {
				return false;
			} else {
				bool f = true;
				for (int i : input_map[sq]) {
					if (!board.pebbled[i]) {
						f = false;
						break;
					}
				}
				return f;
			}
		}
	}
	void add_edge(int from, int to) {
		int sz = std::max(from, to) + 1;
		if (N < sz) {
			N = sz;
			input_map.resize(N);
			output_map.resize(N);
			board.pebbled.resize(N, false);
			retired.resize(N, false);
			exists.resize(N, false);
		}
		exists[from] = exists[to] = true;
		input_map[to].push_back(from);
		output_map[from].push_back(to);
	}
	pebble_game() {
		N = 0;
	}
	int get_score() const {
		int score = board.pebbles.size();
		return score;
	}
	bool pebbled(int i) {
		return board.pebbled[i];
	}
}
;

#include "math.hpp"

class math_pebble_game: public pebble_game {
public:
	size_t operator()(const std::vector<math_vertex> in, const std::vector<math_vertex> out) {
		std::vector<dag_vertex<math_vertex::properties>> vertices;
		std::unordered_map<int, int> index_map;
		for (auto o : out) {
			vertices.push_back(o);
		}
		typename dag_vertex<math_vertex::properties>::executor exe;
		vertices = dag_vertex<math_vertex::properties>::sort(exe, vertices);
		for (int vi = 0; vi < vertices.size(); vi++) {
			auto v = vertices[vi];
			index_map[v.get_unique_id()] = vi;
		}
		for (auto v : vertices) {
			if (v.properties().op != CON) {
				int to = index_map[v.get_unique_id()];
				for (int i = 0; i < v.get_edge_in_count(); i++) {
					const auto& u = v.get_edge_in(i);
					int from = index_map[u.get_unique_id()];
					if (u.properties().op != CON) {
						add_edge(from, to);
					}
				}
			}
		}
		prepare();
		int mscore = std::numeric_limits<int>::max();
		int N = vertices.size();

		ann brain(N, 4);

		int step = 0;
		int best_score = std::numeric_limits<int>::max();
		while (1) {
			int last_mv = -1;
			int score = 0;
			auto last_brain = brain;
			brain.mutate();
			int mvcnt = 0;
			while (!complete()) {
				if (score > best_score) {
					break;
				}
				for (int i = 0; i < N; i++) {
					if (pebbled(i) && legal_move(i)) {
						make_move(i);
						mvcnt++;
					}
				}
				std::vector<bool> input(N);
				for (int n = 0; n < N; n++) {
					if (pebbled(n)) {
						input[n] = true;
					} else {
						input[n] = !legal_move(n);
					}
				}
				auto output = brain.evaluate(input);
				int besti = -1;
				double bestv = -std::numeric_limits<double>::max();
				for (int i = 0; i < N; i++) {
					if (legal_move(i) && i != last_mv) {
						if (output[i] > bestv) {
							bestv = output[i];
							besti = i;
						}
					}
				}
				int mv = besti;
				make_move(mv);
				score = std::max(score, get_score());
				last_mv = mv;
				mvcnt++;
			}
			while (mvcnt) {
				undo_move();
				mvcnt--;
			}
			step++;
			bool code = score < best_score;
			bool code1 = score < best_score;
			if (!code) {
				brain = last_brain;
			} else {
				best_score = score;
			}
			if (code1) {
				fprintf( stderr, "%i. %i %i %c\n", step, score, best_score, code ? '*' : ' ');
			}
		}
		return mscore;
	}
}
;

#endif /* PEBBLE_HPP_ */
