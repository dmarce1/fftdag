/*
 * util.hpp
 *
 *  Created on: Mar 5, 2023
 *      Author: dmarce1
 */

#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <vector>
#include <set>
#include <unordered_map>
#include <functional>

std::vector<std::vector<int>> nchoosek(int n, int k);

template<class T>
std::set<T> intersection(const std::set<T>& A, const std::set<T>& B) {
	std::set<T> C;
	if (B.size() <= A.size()) {
		for (auto b : B) {
			if (A.find(b) != A.end()) {
				C.insert(b);
			}
		}
	} else {
		return intersection(B, A);
	}
	return C;
}

template<class T, class F>
struct set_key {
	size_t operator()(const std::set<T>& set) const {
		F hash;
		size_t key = 1;
		for (auto i = set.begin(); i != set.end(); i++) {
			int bit = key & 1;
			key = *i ^ (key << 1);
			if (bit) {
				key ^= 1;
			}
		}
		return key;
	}
};

/*

 template<class T, class F>
 allint_t<T, F> find_all_intersections(const std::vector<std::set<T>>& sets) {
 allint_t<T, F> rc;
 for (int i = 0; i < sets.size(); i++) {
 std::set<int> s;
 s.insert(i);
 rc[s] = sets[i];
 }
 for (int k = 2; k < sets.size(); k++) {
 bool flag = false;
 std::vector<std::pair<std::set<int>, std::set<T>>>new_entries;
 for (auto entry : rc) {
 auto index = entry.first;
 auto set = entry.second;
 for (int i = 0; i < sets.size(); i++) {
 auto other = sets[i];
 if (index.size() + 1 == k) {
 auto I = intersection(set, other);
 if (I.size() > 1) {
 auto new_index = index;
 new_index.insert(i);
 new_entries.push_back(std::make_pair(new_index, std::move(I)));
 flag = true;
 }
 }
 }
 }
 rc.insert(new_entries.begin(), new_entries.end());
 fprintf(stderr, "%i\n", k);
 if (!flag) {
 break;
 }

 }
 return rc;
 }
 */
template<class T>
std::set<T> operator-(std::set<T> A, const std::set<T>& B) {
	for (auto b : B) {
		A.erase(b);
	}
	return A;
}

template<class T, class F>
using allint_t = std::unordered_map<std::set<int>, std::set<T>, set_key<T, F>>;

template<class T, class F>
allint_t<T, F> find_all_intersections(const std::vector<std::set<T>>& sets) {
	allint_t<T, F> rc;
	const std::function<void(std::set<int> past, std::set<int> future, std::set<T> inter)> find = [&find, &rc, &sets](std::set<int> past, std::set<int> future, std::set<T> inter) {
		rc[past] = inter;
			for( auto i = future.begin(); i != future.end(); i++) {
			auto next_future = future;
			auto next_past = past;
			next_future.erase(*i);
			next_past.insert(*i);
			auto next_inter = intersection(inter, sets[*i]);
			if( next_inter.size() > 1 ) {
				find(std::move(next_past), std::move(next_future), std::move(next_inter));
			}
		}
	};
	std::set<int> future;
	for (int i = 0; i < sets.size(); i++) {
		future.insert(i);
	}
	for (auto i = future.begin(); i != future.end(); i++) {
		auto next_future = future;
		std::set<int> next_past;
		next_future.erase(*i);
		next_past.insert(*i);
		auto inter = sets[*i];
		if (inter.size() > 1) {
			find(std::move(next_past), std::move(next_future), std::move(inter));
		}
	}
	return rc;
}

#endif /* UTIL_HPP_ */
