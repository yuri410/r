#pragma once

#include "Common.h"

class SceneObject
{
public:
	SceneObject(Material* material, IGeometry* geometry);
	~SceneObject();

	IGeometry* GetGeometry() const { return m_geometry; }
	Material* GetMaterial() const { return m_material; }
private:
	Material* m_material = nullptr;
	IGeometry* m_geometry = nullptr;
};
