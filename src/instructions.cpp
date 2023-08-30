
//#define SIMD4
//#define SIMD2
#define SCALAR

#ifdef SIMD4

int simd_width() {
	return 32;
}


const char* mova_op() {
	return "vmovapd";
}

const char* movu_op() {
	return "vmovupd";
}

const char* add_op() {
	return "vaddpd";
}

const char* sub_op() {
	return "vsubpd";
}

const char* mul_op() {
	return "vmulpd";
}

const char* xor_op() {
	return "vxorpd";
}

const char* fma_post() {
	return "pd";
}

const char* simd_reg() {
	return "%ymm";
}


#endif


#ifdef SIMD2

int simd_width() {
	return 16;
}

const char* mova_op() {
	return "vmovapd";
}

const char* movu_op() {
	return "vmovupd";
}

const char* add_op() {
	return "vaddpd";
}

const char* sub_op() {
	return "vsubpd";
}

const char* mul_op() {
	return "vmulpd";
}

const char* xor_op() {
	return "vxorpd";
}

const char* fma_post() {
	return "pd";
}


const char* simd_reg() {
	return "%xmm";
}

#endif

#ifdef SCALAR

int simd_width() {
	return 8;
}

const char* mova_op() {
	return "vmovq";
}

const char* movu_op() {
	return "vmovq";
}

const char* add_op() {
	return "vaddsd";
}

const char* sub_op() {
	return "vsubsd";
}

const char* mul_op() {
	return "vmulsd";
}

const char* xor_op() {
	return "vpxor";
}

const char* fma_post() {
	return "pd";
}


const char* simd_reg() {
	return "%xmm";
}


#endif

