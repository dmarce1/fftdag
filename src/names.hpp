/*
 * names.hpp
 *
 *  Created on: Mar 8, 2023
 *      Author: dmarce1
 */

#ifndef NAMES_HPP_
#define NAMES_HPP_

#include <stack>
#include <string>
#include <memory>
#include <unordered_map>

class math_vertex;

class name_server {
public:
	using name_ptr = std::shared_ptr<std::string>;
	name_server();
	name_ptr generate_name();
	std::pair<name_server::name_ptr, std::string> reserve_name(std::string&& name);
private:
	std::shared_ptr<std::stack<std::string>> available;
	std::shared_ptr<std::unordered_map<std::string, std::weak_ptr<std::string>>> in_use;
	int next_id;
};


#endif /* NAMES_HPP_ */
