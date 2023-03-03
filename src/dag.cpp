#include "dag.hpp"

std::shared_ptr<dag> dag_node::graph = std::make_shared<dag>();
std::map<double, dag::vertex> dag_node::const_map;
std::unordered_map<vertex_type, std::map<dag::vertex, std::map<dag::vertex, dag::vertex>>> dag_node::common_binary;
std::unordered_map<vertex_type, std::map<dag::vertex, dag::vertex>> dag_node::common_unary;
