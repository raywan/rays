#include "bxdf.h"

// Lambertian Reflection
// Lambertian Transmission
// Specular Reflection
// Specular Transmission


// sample_f
// rho_hd
// rho_hh

LambertianReflection lr_create(Vec3 r) {
	LambertianReflection result;
	result.type = (BxDFType) (BSDF_REFLECTION | BSDF_DIFFUSE);
	result.r = r;
	return result;
}

Vec3 lr_sample_f(LambertianReflection *lr) {
	return lr->r * INV_PI;
}

Vec3 lr_rho_hd(LambertianReflection *lr) {
	return lr->r;
}

Vec3 lr_rho_hh(LambertianReflection *lr) {
	return lr->r;
}
