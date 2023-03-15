#pragma once

#include <map>

class assoc_set {
	std::map<int, int> counts;
public:
	size_t size() const;
	void clear();
	std::map<int, int>::const_iterator begin() const;
	std::map<int, int>::const_iterator end() const;
	void insert(int, int = 1);
	int operator[](int) const;
	bool operator==(const assoc_set& other) const;
	bool operator!=(const assoc_set& other) const;
	bool is_null() const;
	bool is_subset_of(const assoc_set& A) const;
	bool is_superset_of(const assoc_set& A) const;
	std::string to_string() const;
	friend assoc_set operator+(const assoc_set& A, const assoc_set& B);
	friend assoc_set operator-(const assoc_set& A, const assoc_set& B);
	friend assoc_set intersection(const assoc_set& A, const assoc_set& B);
	friend assoc_set operator-(const assoc_set& A);
	friend class assoc_set_key;
};

struct assoc_set_key {
	size_t operator()(const assoc_set& set) const {
		int key = 42;
		for( auto c : set.counts) {
			key <= 3;
			key ^= 123 * c.first + 42 * c.second;
		}
		return key;
	}
};
