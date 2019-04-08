#pragma once

#include "Common.h"
#include "Math.h"

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

private:
	PlanarConvexShape<4> m_quad;
};

class DiskShape : public IGeometry
{
public:
	DiskShape(Vec3 center, Vec3 normal, float radius);

	bool Intersects(const Ray& ray, Vec3& point, Vec3& normal) override;
	
private:
	Vec3  m_center;
	Vec3  m_normal;
	float m_radius;
};

class CylinderShape : public IGeometry
{
public:
	CylinderShape(Vec3 center, Vec3 axis, float height, float radius);

	bool Intersects(const Ray& ray, Vec3& point, Vec3& normal) override;

private:
	bool IntersectsCap(const Ray& ray, const Vec3& capCenter, const Vec3& capNormal, Vec3& point, Vec3& normal) const;
	bool IntersectsBody(const Ray& ray, Vec3& point, Vec3& normal) const;

	Vec3  m_center;
	Vec3  m_axis;
	float m_height;
	float m_radius;
};

class TrimeshShape : public IGeometry
{
public:
	TrimeshShape(const std::vector<Triangle>& triangles);
	~TrimeshShape();

	bool Intersects(const Ray& ray, Vec3& point, Vec3& normal) override;

	static TrimeshShape* GenerateNoiseQuad(Vec3 pos, double xSpan, double ySpan, int xSegments, int ySegments, const PerlinNoise& noise);

private:
	struct BvhNode final
	{
		BoundingSphere m_boundingSphere;

		std::vector<BvhNode> m_nodes;
		std::vector<Triangle> m_triangles;

		BvhNode() { }
		BvhNode(const std::vector<Triangle>& triangles);
		
		BvhNode(BvhNode&&);
		BvhNode& operator=(BvhNode&&);

		BvhNode(const BvhNode&) = delete;
		BvhNode& operator=(const BvhNode&) = delete;

		bool Intersects(const Ray& ray, Vec3& point, Vec3& normal, float& nearestDist) const;

	};

	BvhNode m_rootNode;
};