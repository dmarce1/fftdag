#include "names.hpp"
#include <cassert>
#include <limits>

name_server::name_server() {
	next_id = 0;
	available = std::make_shared<std::set<std::string>>();
	in_use = std::make_shared<std::unordered_map<std::string, std::weak_ptr<std::string>>>();
	reg2mem = std::make_shared<std::unordered_map<std::string, std::string>>();
	mem2reg = std::make_shared<std::unordered_map<std::string, std::string>>();
	avail_regs = std::make_shared<std::set<std::string>>();
	reg_q = std::make_shared<std::queue<std::string>>();
	for (int n = 0; n < 16; n++) {
		avail_regs->insert(std::string("%ymm") + std::to_string(n));
	}
}

std::string name_server::get_declarations() const {
	return declarations;
}

std::string name_server::current_name(std::string nm) const {
	if (mem2reg->find(nm) != mem2reg->end()) {
		return mem2reg->find(nm)->second;
	} else {
		return nm;
	}
}

std::string name_server::get_register(std::string mem, std::string& code, bool noload) {
	assert(mem[0] != 'C');
	assert(mem[0] != '%');
	if (mem2reg->find(mem) != mem2reg->end()) {
		return (*mem2reg)[mem];
	} else {
		char* ptr;
		std::string reg;
		if (avail_regs->size()) {
			reg = *avail_regs->begin();
			avail_regs->erase(reg);
			(*mem2reg)[mem] = reg;
			(*reg2mem)[reg] = mem;
			std::queue<std::string> newq;
			while (reg_q->size()) {
				if (reg_q->front() != reg) {
					newq.push(reg_q->front());
				}
				reg_q->pop();
			}
			*reg_q = std::move(newq);
			reg_q->push(reg);
		} else {
			if (!reg_q->size()) {
				assert(false);
				abort();
			}
			reg = reg_q->front();
			reg_q->pop();
			if (reg2mem->find(reg) != reg2mem->end()) {
				asprintf(&ptr, "%15s%-15s%s, %s\n", "", "vmovupd", reg.c_str(), (*reg2mem)[reg].c_str());
				code += ptr;
				free(ptr);
				mem2reg->erase((*reg2mem)[reg]);
				reg2mem->erase(reg);
			}
			(*mem2reg)[mem] = reg;
			(*reg2mem)[reg] = mem;
			reg_q->push(reg);
		}
		if (!noload) {
			asprintf(&ptr, "%15s%-15s%s, %s\n", "", "vmovupd", mem.c_str(), reg.c_str());
			code += ptr;
			free(ptr);
		}
		return reg;
	}
}

name_server::name_ptr name_server::generate_name() {
	if (available->empty()) {
		auto new_name = std::to_string(-32 * (1 + next_id)) + std::string("(%rbp)");
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
		if( mem2reg->find(*ptr) != mem2reg->end()) {
			reg2mem->erase((*mem2reg)[*ptr]);
			avail_regs->insert((*mem2reg)[*ptr]);
			mem2reg->erase(*ptr);
		}
		delete ptr;
	});
	assert(in_use->find(name) == in_use->end());
	(*in_use)[name] = nptr;
	return nptr;
}

/*std::pair<name_server::name_ptr, std::string> name_server::reserve_name(std::string name) {
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
 */
int distance(const std::string& a, const std::string& b) {
	int d = std::numeric_limits<int>::max() - 1;
	if (a[0] == b[0] && a[1] == b[1]) {
		d = std::abs(atoi(a.c_str() + 2) - atoi(b.c_str() + 2));
	}
	return d;
}

