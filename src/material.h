#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "primitive.h"

// enum MaterialType {
// 	M_MATTE,
// };

struct MaterialMatte {
	Vec3 Kd;
	float sigma; // For Oren-Nayar
};

MaterialMatte m_matte_create(Vec3 Kd, float sigma);
void m_matte_compute_scattering(MaterialMatte *mm, IntersectInfo *ii);


#endif
