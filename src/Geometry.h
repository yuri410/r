#pragma once

#include "Common.h"

class IGeometry
{
public:
	virtual ~IGeometry() { }
	
	virtual bool Intersects(const Ray& ray, Vec3& point, Vec3& normal) = 0;
};

class SphereShape : public IGeometry
{
public:
	SphereShape(Vec3 center, float radius);

	bool Intersects(const Ray& ray, Vec3& point, Vec3& normal) override;

private:
	Vec3  m_center;
	float m_radius;
};

class QuadShape : public IGeometry
{
public:
	QuadShape(Vec3 center, Vec3 normal, Vec3 tangent, float width, float height);

	bool Intersects(const Ray& ray, Vec3& point, Vec3& normal) override;

	const Vec3* GetCorners() const { return m_corners; }
private:
	Vec3  m_corners[4];
	Vec3  m_normal;
	Vec3  m_edgeNormal[4];
};