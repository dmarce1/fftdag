#include "dag.hpp"

std::shared_ptr<dag> dag_node::graph = std::make_shared<dag>();
std::map<double, dag::vertex> dag_node::const_map;

