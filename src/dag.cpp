#include "dag.hpp"

int dag::next_id = 0;
std::shared_ptr<dag> dag_node::graph = std::make_shared<dag>();
