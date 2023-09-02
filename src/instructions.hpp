/*
 * instructions.hpp
 *
 *  Created on: Aug 30, 2023
 *      Author: dmarce1
 */

#ifndef INSTRUCTIONS_HPP_
#define INSTRUCTIONS_HPP_



int simd_width();
const char* mova_op();
const char* movu_op();
const char* add_op();
const char* sub_op();
const char* mul_op();
const char* xor_op();
const char* fma_post();
const char* simd_reg();
std::string  simd_reg(int i);
void set_simd_width(int);

#endif /* INSTRUCTIONS_HPP_ */
