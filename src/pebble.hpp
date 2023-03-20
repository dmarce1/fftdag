/*
 * pebble.hpp
 *
 *  Created on: Mar 19, 2023
 *      Author: dmarce1
 */

#ifndef PEBBLE_HPP_
#define PEBBLE_HPP_

#include "dag.hpp"

template<class T>
class pebble_game {
	using pebbles_type = std::unordered_set<int>;
	pebbles_type pebbles;
	struct key_type {
		size_t operator()(const pebbles_type& pebbles) const {
			std::hash<int> hash;
			size_t key = 0;
			for (auto pebble : pebbles) {
				key ^= hash(pebble);
			}
			return key;
		}
	};
	std::vector<dag_vertex<T>> inputs;
	std::vector<dag_vertex<T>> outputs;
	std::unordered_map<int, std::vector<int>> output_map;
	std::unordered_map<int, std::vector<int>> input_map;
	std::unordered_map<pebbles_type, int, key_type> score_memory;
	std::unordered_set<int> prohibited_moves;
	std::unordered_set<int> prohibited_kills;
	std::vector<int> move_list;
	struct move_undo_t {
		std::vector<int> kills;
		int move;
	};
	std::stack<move_undo_t> move_stack;

	std::vector<int> get_moves() {
		struct entry_type {
			int move;
			int score;
		};
		std::vector<int> moves;
		std::vector<entry_type> move_entries;
		for (auto& input : input_map) {
			if (prohibited_moves.find(input.first) == prohibited_moves.end()) {
				if (pebbles.find(input.first) == pebbles.end()) {
					bool flag = false;
					for (auto i : input.second) {
						if (pebbles.find(i) == pebbles.end()) {
							flag = true;
							break;
						}
					}
					if (!flag) {
						entry_type entry;
						entry.move = input.first;
						entry.score = 0;
						for (auto i : input.second) {
							const auto& outs = output_map[i];
							bool flag = false;
							for (auto o : outs) {
								if (o != input.first) {
									if (pebbles.find(o) == pebbles.end()) {
										flag = true;
										break;
									}
								}
							}
							if (!flag) {
								entry.score++;
							}
						}
					//	fprintf(stderr, "%i\n", entry.score);
						move_entries.push_back(entry);
					}
				}
			}
		}
		std::sort(move_entries.begin(), move_entries.end(), [](const entry_type& a, const entry_type& b) {
			return a.score > b.score;
		});
		for (auto entry : move_entries) {
			moves.push_back(entry.move);
		}
		return std::move(moves);
	}

	void make_move(int mv) {
		move_undo_t undo;
		std::vector<int> kills;
		move_list.push_back(mv);
		for (auto pebble : pebbles) {
			if (prohibited_kills.find(pebble) == prohibited_kills.end()) {
				auto& outputs = output_map[pebble];
				bool flag = false;
				for (auto o : outputs) {
					if (pebbles.find(o) == pebbles.end()) {
						flag = true;
						break;
					}
				}
				if (!flag) {
					assert(pebble != mv);
					kills.push_back(pebble);
				}
			}
		}
		for (auto kill : kills) {
			assert(pebbles.find(kill) != pebbles.end());
			pebbles.erase(kill);
		}
		assert(pebbles.find(mv) == pebbles.end());
		assert(prohibited_moves.find(mv) == prohibited_moves.end());
		pebbles.insert(mv);
		prohibited_moves.insert(mv);
		undo.kills = std::move(kills);
		undo.move = mv;
		move_stack.push(std::move(undo));
	}

	void undo_move() {
		auto undo = std::move(move_stack.top());
		move_stack.pop();
		for (auto kill : undo.kills) {
			assert(pebbles.find(kill) == pebbles.end());
			pebbles.insert(kill);
		}
		assert(prohibited_moves.find(undo.move) != prohibited_moves.end());
		assert(pebbles.find(undo.move) != pebbles.end());
		prohibited_moves.erase(undo.move);
		pebbles.erase(undo.move);
		move_list.pop_back();
	}

	int get_score() {
		return pebbles.size();
	}

public:

	int search(int best_score) {
		int score;
		auto mem_iter = score_memory.find(pebbles);
		int lower_bound = std::numeric_limits<int>::min();
		if (mem_iter != score_memory.end()) {
			lower_bound = mem_iter->second;
		}
		if (lower_bound >= best_score) {
			score = lower_bound;
		} else {
			score = get_score();
			if (score < best_score) {
				auto moves = get_moves();
				if (moves.size()) {
					for (auto move : moves) {
						make_move(move);
						score = search(best_score);
						undo_move();
						if (score < best_score) {
							best_score = score;
						}
					}
					score = std::max(score, best_score);
				}
			}
		}
		score_memory[pebbles] = std::max(lower_bound, score);
		if (rand() % 1000 == 0) {
			for (auto mv : move_list) {
				fprintf( stderr, "%i ", mv);
			}
			fprintf( stderr, " : %i\n", score);
		}
		return score;
	}

	int search() {
		int score;
		int test_score = std::max(outputs.size(), inputs.size());
		do {
			test_score++;
			score = search(test_score);
			fprintf( stderr, "%i %i\n", score, test_score);
		} while (score >= test_score);
		return score;
	}

	pebble_game(const std::vector<dag_vertex<T>>& i, const std::vector<dag_vertex<T>>& j) {
		inputs = i;
		outputs = j;
		typename dag_vertex<T>::executor exe;
		const auto vertices = dag_vertex<T>::sort(exe, outputs);
		for (auto vertex : vertices) {
			const auto out = vertex.get_unique_id();
			for (int i = 0; i < vertex.get_edge_in_count(); i++) {
				const auto in = vertex.get_edge_in(i).get_unique_id();
				output_map[in].push_back(out);
				input_map[out].push_back(in);
			}
		}
		for (auto& input : inputs) {
			pebbles.insert(input.get_unique_id());
		}
		for (auto& i : inputs) {
			prohibited_moves.insert(i.get_unique_id());
		}
		for (auto& o : outputs) {
			prohibited_kills.insert(o.get_unique_id());
		}
	}
}
;

#endif /* PEBBLE_HPP_ */
