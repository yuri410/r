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

	Material() {}
	Material(Vec3 albedo)
		: m_albedo(albedo) { }
	Material(Vec3 albedo, Vec3 emissive, float roughness, MaterialType type = MaterialType::Diffuse)
		: m_albedo(albedo), m_emissive(emissive), m_roughness(roughness), m_type(type) { }
};