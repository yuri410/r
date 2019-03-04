#pragma once

#include "Common.h"

enum struct MaterialType
{
	Diffuse,
	Specular,
};

struct Material
{
	Vec3 m_albedo;
	Vec3 m_emissive;
	float m_roughness;

	MaterialType m_type = MaterialType::Diffuse;
};