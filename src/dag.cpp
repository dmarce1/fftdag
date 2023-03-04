#include "dag.hpp"

std::shared_ptr<dag> dag_node::graph = std::make_shared<dag>();
std::map<double, dag::vertex> dag_node::const_map;

std::set<int> intersection(const std::set<int>& A, const std::set<int>& B) {
	std::set<int> C;
	for (auto b : B) {
		if (A.find(b) != A.end()) {
			C.insert(b);
		}
	}
	return C;
}
