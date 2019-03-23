#pragma once
#include "Math.h"

namespace BRDF
{
	// GGX
	inline float NDF(float NdotH, float roughness)
	{
		float pi = 3.1415927f;

		float a2 = roughness*roughness;
		return a2 / (pi * Sqr((a2 - 1)*NdotH*NdotH + 1));
	}

	// Schlick
	inline float GSF(float NdotL, float NdotV, float roughness)
	{
		float k = roughness / 2;

		float SmithL = (NdotL) / (NdotL * (1 - k) + k);
		float SmithV = (NdotV) / (NdotV * (1 - k) + k);

		return (SmithL * SmithV);
	}

	inline Vec3 F(Vec3 f0, float LdotH)
	{
		return f0 + (Vec3(1.0f, 1.0f, 1.0f) - f0) * pow(LdotH, 5);
	}

	inline Vec3 CookTorrance(Vec3 f, float d, float g, float ndl, float ndv)
	{
		return f * ((d * g) / (4 * ndl * ndv));
	}
}
