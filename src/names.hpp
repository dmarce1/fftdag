/*
 * names.hpp
 *
 *  Created on: Mar 8, 2023
 *      Author: dmarce1
 */

#ifndef NAMES_HPP_
#define NAMES_HPP_

#include <set>
#include <string>
#include <memory>
#include <unordered_map>
#include <queue>
class math_vertex;

class name_server {
public:
	using name_ptr = std::shared_ptr<std::string>;
	name_server();
	name_ptr generate_name();
	std::string get_register(std::string, std::string&, bool);
	//std::pair<name_server::name_ptr, std::string> reserve_name(std::string name);
	std::string get_declarations() const;
	std::string current_name(std::string nm) const;
	void transfer_register(std::string reg, std::string from, std::string to) {
		(*reg2mem)[reg] = to;
		(*mem2reg)[to] = reg;
		mem2reg->erase(from);
	}
	bool name_in_use(std::string nm) const {
		return in_use->find(nm) != in_use->end();
	}
	int size() const {
		return next_id;
	}
private:
	std::shared_ptr<std::set<std::string>> available;
	std::shared_ptr<std::unordered_map<std::string, std::weak_ptr<std::string>>> in_use;
	std::string declarations;
	int next_id;
	std::shared_ptr<std::unordered_map<std::string, std::string>> reg2mem;
	std::shared_ptr<std::unordered_map<std::string, std::string>> mem2reg;
	std::shared_ptr<std::set<std::string>> avail_regs;
	std::shared_ptr<std::queue<std::string>> reg_q;
};

int distance(const std::string& a, const std::string& b);


#endif /* NAMES_HPP_ */
