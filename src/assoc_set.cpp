#include "assoc_set.hpp"

void assoc_set::insert(int val, int cnt) {
	if (val < 0) {
		insert(-val, -cnt);
	} else {
		auto iter = counts.find(val);
		if (iter == counts.end()) {
			counts[val] = 0;
			iter = counts.find(val);
		}
		iter->second += cnt;
		if (iter->second == 0) {
			counts.erase(iter->second);
		}
	}
}

int assoc_set::operator[](int val) const {
	auto iter = counts.find(val);
	if (iter == counts.end()) {
		return 0;
	} else if (val < 0) {
		return -operator[](-val);
	} else {
		return iter->second;
	}
}

assoc_set operator+(const assoc_set& A, const assoc_set& B) {
	assoc_set C = A;
	for (auto i : B.counts) {
		if (C.counts.find(i.first) == C.counts.end()) {
			C.counts[i.first] = 0;
		}
		C.counts[i.first] += i.second;
		if( C.counts[i.first] == 0 ) {
			C.counts.erase(i.first);
		}
	}
	return C;
}

size_t assoc_set::size() const {
	return counts.size();
}

assoc_set operator-(const assoc_set& A, const assoc_set& B) {
	assoc_set C;
	for (auto i : A.counts) {
		C.counts[i.first] = i.second;
	}
	for (auto i : B.counts) {
		if (C.counts.find(i.first) == C.counts.end()) {
			C.counts[i.first] = 0;
		}
		C.counts[i.first] -= i.second;
		if( C.counts[i.first] == 0 ) {
			C.counts.erase(i.first);
		}
	}
	return C;
}

assoc_set operator*(const assoc_set& A, const assoc_set& B) {
	assoc_set C;
	for (auto i : A.counts) {
		C.counts[i.first] = i.second;
	}
	for (auto i : B.counts) {
		if (C.counts.find(i.first) != C.counts.end()) {
			C.counts[i.first] *= i.second;
		}
	}
	return C;
}

assoc_set operator-(const assoc_set& A) {
	assoc_set B;
	for (auto i : A.counts) {
		B.counts[i.first] = -i.second;
	}
	return B;
}

bool assoc_set::operator==(const assoc_set& other) const {
	return counts == other.counts;
}

std::map<int, int>::const_iterator assoc_set::begin() const {
	return counts.begin();
}

std::map<int, int>::const_iterator assoc_set::end() const {
	return counts.end();
}

bool assoc_set::operator!=(const assoc_set& other) const {
	return !(*this == other);
}

bool assoc_set::is_null() const {
	return counts.size() == 0;
}

bool assoc_set::is_subset_of(const assoc_set& A) const {
	return counts.size() == (*this + A).counts.size();
}

bool assoc_set::is_superset_of(const assoc_set& A) const {
	return A.is_subset_of(*this);
}

