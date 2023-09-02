#include <string>
//#define SIMD4
//#define SIMD2
#define SCALAR

static int width = 4;

void set_simd_width(int w) {
	width = w;
}

int simd_width() {
	return 8 * width;
}

const char* mova_op() {
	return width == 1 ? "vmovq" : "vmovapd";
}

const char* movu_op() {
	return width == 1 ? "vmovq" : "vmovupd";
}

const char* add_op() {
	return width == 1 ? "vaddsd" : "vaddpd";
}

const char* sub_op() {
	return width == 1 ? "vsubsd" : "vsubpd";
}

const char* mul_op() {
	return width == 1 ? "vmulsd" : "vmulpd";
}

const char* xor_op() {
	return width == 1 ? "vpxor" : "vxorpd";
}

const char* fma_post() {
	return width == 1 ? "sd" : "pd";
}

const char* simd_reg() {
	return width == 4 ? "%ymm" : "%xmm";
}

std::string simd_reg(int i) {
	return std::string(simd_reg()) + std::to_string(i);
}

