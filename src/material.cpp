#include "material.h"

MaterialMatte m_matte_create(Vec3 Kd, float sigma) {
	MaterialMatte result;
	result.Kd = Kd;
	result.sigma = sigma;
	return result;
}

void m_matte_compute_scattering(MaterialMatte *mm, IntersectInfo *ii) {

}

