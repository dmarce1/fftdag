#include "dag.hpp"

std::shared_ptr<math_dag> dag_node::graph = std::make_shared<math_dag>();
std::map<double, dag<math_props>::vertex> dag_node::const_map;

std::set<int> intersection(const std::set<int>& A, const std::set<int>& B) {
	std::set<int> C;
	for (auto b : B) {
		if (A.find(b) != A.end()) {
			C.insert(b);
		}
	}
	return C;
}


std::set<int> antiset(const std::set<int>& A) {
	std::set<int> B;
	for( auto a : A ) {
		B.insert(-a);
	}
	return B;
}
