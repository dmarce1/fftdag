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

class pebble_game {
	struct pebbles_type {
		std::vector<bool> pebbled;
		std::unordered_set<int> pebbles;
		size_t key1;
		size_t key2;
		bool operator==(const pebbles_type& other) const {
			return (key1 == other.key1) && (key2 == other.key2);
		}
	};
	struct mem_entry_key {
		size_t operator()(const pebbles_type& entry) const {
			return entry.key1 * entry.key2;
		}
	};
	int N;
	std::vector<std::vector<int>> input_map;
	std::vector<std::vector<int>> output_map;
	std::vector<bool> exists;
	std::vector<int> stack;
	std::vector<int> keys1;
	std::vector<int> keys2;
	std::unordered_set<int> inputs;
	std::unordered_set<int> outputs;
	std::vector<bool> retired;
	std::vector<int> distance;
	std::unordered_map<pebbles_type, int, mem_entry_key> prev_boards;
	int answer;
	pebbles_type board;
	std::unordered_map<pebbles_type, std::pair<int, int>, mem_entry_key> memory;
public:
	std::vector<int> generate_moves() {
		std::vector<int> moves;
		std::vector<int> priority;
		static int max_moves = 0;
		moves.reserve(max_moves);
		for (auto pebble : board.pebbles) {
			bool flag = true;
			bool interior = false;
			for (auto o : output_map[pebble]) {
				interior = true;
				if (!board.pebbled[o]) {
					flag = false;
				}
			}
			if (flag && interior) {
				moves.push_back(pebble);
				break;
			}
		}
		if (!moves.size()) {
			//		priority.reserve(max_moves);
			std::unordered_set<int> used;
			for (auto pebble : board.pebbles) {
				for (auto o : output_map[pebble]) {
					if (used.find(o) == used.end() && !retired[o]) {
						used.insert(o);
						if (!board.pebbled[o]) {
							bool flag = true;
							bool prior = false;
							for (auto i : input_map[o]) {
								if (!board.pebbled[i]) {
									flag = false;
									break;
								}
								//	int cnt = 0;
								//	for (auto o : output_map[i]) {
								//		if (!board.pebbled[o]) {
								//				cnt++;
								//				}
								//				}
								//					if (cnt == 1) {
								//							prior = true;
								//					}
							}
							if (flag) {
								//	if (prior) {
								//		priority.push_back(o);
								//	} else {
								moves.push_back(o);
								//	}
							}
						}
					}
				}
			}
		}
//		priority.insert(priority.end(), moves.begin(), moves.end());
//		moves = std::move(priority);
		max_moves = std::max(max_moves, (int) moves.size());
		return std::move(moves);
	}

	void swap_pebble(int sq) {
		if (board.pebbled[sq]) {
			board.pebbled[sq] = false;
			board.pebbles.erase(sq);
		} else {
			board.pebbled[sq] = true;
			board.pebbles.insert(sq);
		}
		board.key1 ^= keys1[sq];
		board.key2 ^= keys2[sq];
	}
	bool make_move(int mv) {
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
		distance.resize(N, 0);
		std::unordered_set<int> touched;
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
			if (!board.pebbled[oi]) {
				return false;
			}
		}
		return true;
	}
	int search(int depth, int lev, int alpha = std::numeric_limits<int>::min(), int beta = std::numeric_limits<int>::max()) {
		beta = std::min(beta, get_score());
		if (alpha >= beta) {
			return alpha;
		}
		if (depth == 0) {
			return get_score();
		}
		int upper_bound = std::numeric_limits<int>::max();
		int lower_bound = std::numeric_limits<int>::min();
		if (memory.find(board) != memory.end()) {
			auto& entry = memory[board];
			upper_bound = entry.second;
			lower_bound = entry.first;
		} else {
			auto& entry = memory[board];
			entry.second = std::numeric_limits<int>::max();
			entry.first = std::numeric_limits<int>::min();
		}
		if (upper_bound <= alpha) {
			return upper_bound;
		} else {
			beta = std::min(beta, upper_bound);
		}
		if (lower_bound >= beta) {
			return beta;
		} else {
			alpha = std::max(alpha, lower_bound);
		}
		auto moves = generate_moves();
		if (moves.size()) {
			for (auto mv : moves) {
				make_move(mv);
				int this_score = search(depth - 1, lev + 1, alpha, beta);
				undo_move();
				if (this_score > alpha) {
					alpha = this_score;
					if (lev == 0) {
						answer = mv;
					}
				}
				if (alpha >= beta) {
					alpha = beta;
					break;
				}
			}
			auto& entry = memory[board];
			if (alpha >= beta) {
				entry.first = std::max(entry.first, alpha);
			} else {
				entry.second = std::min(entry.second, beta);
			}
			return alpha;
		} else {
			return get_score();
		}
	}
	void clear_memory() {
		//	memory.clear();
		for (auto i = memory.begin(); i != memory.end(); i++) {
			i->second.first = std::numeric_limits<int>::min();
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
			int M = keys1.size();
			keys1.resize(N);
			for (int i = M; i < N; i++) {
				keys1[i] = rand() * rand();
			}
			keys2.resize(N);
			for (int i = M; i < N; i++) {
				keys2[i] = rand() * rand();
			}
		}
		exists[from] = exists[to] = true;
		input_map[to].push_back(from);
		output_map[from].push_back(to);
	}
	pebble_game() {
		board.key1 = 0;
		board.key2 = 0;
		N = 0;
	}
	int get_score() const {
		int score = -board.pebbles.size() * 100;
		int minpeb = std::numeric_limits<int>::max();
		int maxpeb = std::numeric_limits<int>::min();
		for (auto pebble : board.pebbles) {
			minpeb = std::min(minpeb, pebble);
			maxpeb = std::max(maxpeb, pebble);
			int pcnt = 0;
			int ocnt = 0;
			for (auto o : output_map[pebble]) {
				if (board.pebbled[o]) {
					pcnt++;
				} else {
					ocnt++;
				}
			}
			if (pcnt != 0 && ocnt != 0) {
				score -= (21 + (ocnt - 1) * 11);
			} else if (pcnt == 0 && ocnt == 0) {
				score -= 3;
			}
		}
		score -= (maxpeb - minpeb) * 50 / output_map.size();
		return score;
	}
	int pebble_count() const {
		return board.pebbles.size();
	}
	int get_best_move() const {
		return answer;
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
		int score;
		int pmax = 0;
		while (!complete()) {
			int last_score;
			int depth;
			int cnt = 0;
			score = search(10, 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
			mscore = std::min(mscore, score);
			int mv = get_best_move();
			make_move(mv);
			pmax = std::max(pmax, pebble_count());
			printf("%i %i %i pcmax: %i pc: %i d: %i\n", mv, score, mscore, pmax, pebble_count(), depth);
			clear_memory();
		}
		return mscore;
	}
}
;

#endif /* PEBBLE_HPP_ */
