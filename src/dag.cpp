#include "dag.hpp"

std::shared_ptr<dag> dag_node::graph = std::make_shared<dag>();
/*std::map<double, vertex> dag_node::const_map;

inline vertex::~vertex() {
	if (id != nullptr) {
		if (id.use_count() == 1) {
			graph->remove_vertex(*id);
		}
	}
}*/
