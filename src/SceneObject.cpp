#include "SceneObject.h"
#include "Geometry.h"
#include "Material.h"

SceneObject::SceneObject(Material* material, IGeometry* geometry)
	: m_material(material), m_geometry(geometry)
{

}

SceneObject::~SceneObject()
{
	delete m_material;
	m_material = nullptr;

	delete m_geometry;
	m_geometry = nullptr;
}
