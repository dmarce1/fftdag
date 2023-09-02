#include "names.hpp"
#include <cassert>
#include <limits>
#include <algorithm>
#include "instructions.hpp"

name_server::name_server() {
	next_id = 0;
	available = std::make_shared<std::set<std::string>>();
	in_use = std::make_shared<std::unordered_map<std::string, std::weak_ptr<std::string>>>();
	reg2mem = std::make_shared<std::unordered_map<std::string, std::string>>();
	mem2reg = std::make_shared<std::unordered_map<std::string, std::string>>();
	avail_regs = std::make_shared<std::set<std::string>>();
	reg_q = std::make_shared<std::queue<std::string>>();
	for (int n = 0; n < 14; n++) {
		avail_regs->insert(std::string(simd_reg()) + std::to_string(n));
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
	//assert(mem.size());
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
				asprintf(&ptr, "%15s%-15s%s, %s\n", "", mova_op(), reg.c_str(), (*reg2mem)[reg].c_str());
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
			asprintf(&ptr, "%15s%-15s%s, %s\n", "", mova_op(), mem.c_str(), reg.c_str());
			code += ptr;
			free(ptr);
		}
		return reg;
	}
}



name_server::name_ptr name_server::reserve_name(std::string name) {
	auto ptr = new std::string(name);
	auto nptr = std::shared_ptr<std::string>(ptr, [this](std::string* ptr) {
		in_use->erase(*ptr);
		available->insert(*ptr);
		if( mem2reg->find(*ptr) != mem2reg->end()) {
			reg2mem->erase((*mem2reg)[*ptr]);
			avail_regs->insert((*mem2reg)[*ptr]);
			mem2reg->erase(*ptr);
		}
		delete ptr;
	});
	(*in_use)[name] = nptr;
	return nptr;

}

name_server::name_ptr name_server::change_name(std::string& code, std::string from, std::string to) {
	char* ptr;
	if (from == to) {
		return name_server::name_ptr((*in_use)[from]);
	}
	if (in_use->find(to) != in_use->end()) {
		std::swap(*name_server::name_ptr((*in_use)[to]), *name_server::name_ptr((*in_use)[from]));
		asprintf(&ptr, "%15s%-15s%s, %s\n", "", mova_op(), current_name(from).c_str(), simd_reg(14).c_str());
		code += ptr;
		free(ptr);
		asprintf(&ptr, "%15s%-15s%s, %s\n", "", mova_op(), current_name(to).c_str(), simd_reg(15).c_str());
		code += ptr;
		free(ptr);
		asprintf(&ptr, "%15s%-15s%s, %s\n", "", mova_op(), simd_reg(15).c_str(), current_name(from).c_str());
		code += ptr;
		free(ptr);
		asprintf(&ptr, "%15s%-15s%s, %s\n", "", mova_op(), simd_reg(14).c_str(), current_name(to).c_str());
		code += ptr;
		free(ptr);
		if( (*mem2reg)[from] == (*mem2reg)[to]) {
			std::swap((*mem2reg)[from], (*mem2reg)[to]);
			(*reg2mem)[(*mem2reg)[from]] = from;
			(*reg2mem)[(*mem2reg)[to]] = to;
		}
		std::swap(*name_server::name_ptr((*in_use)[to]), *name_server::name_ptr((*in_use)[from]));
		return name_server::name_ptr((*in_use)[to]);
	} else if (available->find(to) != available->end()) {
		available->erase(to);
		auto nptr = std::shared_ptr<std::string>(new std::string(to), [this](std::string* ptr) {
			if(in_use->find(*ptr) == in_use->end()) {
				printf( "Can't find %s to erase it\n", ptr->c_str());
				assert(false);
			}
			in_use->erase(*ptr);
			available->insert(*ptr);
			if( mem2reg->find(*ptr) != mem2reg->end()) {
				reg2mem->erase((*mem2reg)[*ptr]);
				avail_regs->insert((*mem2reg)[*ptr]);
				mem2reg->erase(*ptr);
			}
			delete ptr;
		});
		(*in_use)[to] = nptr;
		asprintf(&ptr, "%15s%-15s%s, %s\n", "", mova_op(), current_name(from).c_str(), simd_reg(15).c_str());
		code += ptr;
		free(ptr);
		asprintf(&ptr, "%15s%-15s%s, %s\n", "", mova_op(), simd_reg(15).c_str(), current_name(to).c_str());
		code += ptr;
		free(ptr);
		return nptr;
	} else {
		assert(false);
	}

}

std::string name_server::free_regs() {
	std::string code;
	char* ptr;
	for (auto i = reg2mem->begin(); i != reg2mem->end(); i++) {
		auto reg = (*mem2reg)[i->second];
		asprintf(&ptr, "%15s%-15s%s, %s\n", "", mova_op(), reg.c_str(), (*reg2mem)[reg].c_str());
		code += ptr;
		free(ptr);
	}
	return code;
}

name_server::name_ptr name_server::generate_name() {
	if (available->empty()) {
		auto new_name = std::string("MEM") + &((std::to_string(next_id + 1000))[1]);
		available->insert(std::move(new_name));
		next_id++;
	}
	std::string name = *(available->begin());
	available->erase(name);
	assert(in_use->find(name) == in_use->end());
	auto ptr = new std::string(name);
	auto nptr = std::shared_ptr<std::string>(ptr, [this](std::string* ptr) {
		if(in_use->find(*ptr) == in_use->end()) {
			printf( "Can't find %s to erase it\n", ptr->c_str());
			assert(false);
		}
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

int distance(const std::string& a, const std::string& b) {
	int d = std::numeric_limits<int>::max() - 1;
	if (a[0] == b[0] && a[1] == b[1]) {
		d = std::abs(atoi(a.c_str() + 2) - atoi(b.c_str() + 2));
	}
	return d;
}

