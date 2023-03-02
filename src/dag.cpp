#include "dag.hpp"

std::shared_ptr<dag> dag_node::graph = std::make_shared<dag>();
std::map<double, int> dag_node::const_map;
