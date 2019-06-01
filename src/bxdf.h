#ifndef __BXDF_H__
#define __BXDF_H__

#include <rw/rw_math.h>

enum BxDFType {
	BSDF_REFLECTION = 1 << 0,
	BSDF_TRANSMISSION = 1 << 1,
  BSDF_DIFFUSE = 1 << 2,
  BSDF_GLOSSY = 1 << 3,
  BSDF_SPECULAR = 1 << 4,
  BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR | BSDF_REFLECTION | BSDF_TRANSMISSION,
};

struct LambertianReflection {
	BxDFType type;
	Vec3 r; // reflectance
};

// struct SpecularReflection {
// 	Vec3 r; // reflectance
// 	Fresnel
// };

// struct SpecularTransmission {
// 	Vec3 t;
// 	float eta_a;
// 	float lta_b;
// };

// struct FresnelSpecular {

// };

struct BxDF_SOA {
	LambertianReflection lr[10];
	// SpecularReflection sr[10];
	// SpecularTransmission st[10];
	// FresnelSpecular fs[10];
};

struct BSDF {
	int bxdf_table_idx[10];
};

#endif
