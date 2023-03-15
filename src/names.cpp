#include "names.hpp"
#include <cassert>

name_server::name_server() {
	next_id = 0;
	available = std::make_shared<std::set<std::string>>();
	in_use = std::make_shared<std::unordered_map<std::string, std::weak_ptr<std::string>>>();
}

std::string name_server::get_declarations() const {
	return declarations;
}

name_server::name_ptr name_server::generate_name() {
	if (available->empty()) {
		auto new_name = std::string("r") + std::to_string(next_id);
		declarations += std::string("\tdouble ") + new_name + ";\n";
		available->insert(std::move(new_name));
		next_id++;
	}
	std::string name = *(available->begin());
	available->erase(name);
	auto ptr = new std::string(name);
	auto nptr = std::shared_ptr<std::string>(ptr, [this](std::string* ptr) {
		assert(in_use->find(*ptr) != in_use->end());
		in_use->erase(*ptr);
		available->insert(*ptr);
		delete ptr;
	});
	assert(in_use->find(name) == in_use->end());
	(*in_use)[name] = nptr;
	return nptr;
}

std::pair<name_server::name_ptr, std::string> name_server::reserve_name(std::string name) {
	std::string code;
	if (in_use->find(name) != in_use->end()) {
		auto other_ptr = name_ptr((*in_use)[name]);
		auto tmp = generate_name();
		std::swap(*tmp, *other_ptr);
		(*in_use)[*tmp] = tmp;
		(*in_use)[*other_ptr] = other_ptr;
		code = std::string("\t") + *other_ptr + " = " + *tmp + ";\n";
	}
	available->erase(name);
	auto ptr = new std::string(name);
	auto nptr = std::shared_ptr<std::string>(ptr, [this](std::string* ptr) {
		assert(in_use->find(*ptr) != in_use->end());
		in_use->erase(*ptr);
		available->insert(*ptr);
		delete ptr;
	});
	assert(in_use->find(name) == in_use->end());
	(*in_use)[name] = nptr;
	return std::make_pair(nptr, code);
}

