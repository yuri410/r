#include "Geometry.h"

SphereShape::SphereShape(Vec3 center, float radius)
	: m_center(center), m_radius(radius)
{ }

bool SphereShape::Intersects(const Ray& ray, Vec3& point, Vec3& normal)
{
	Vec3 sc = m_center - ray.m_origin;
	float scLen = sc.Length();
	sc = sc * (1.0f / scLen);

	float cosTheta = ray.m_direction.Dot(sc);

	if (cosTheta > 0)
	{
		float sinTheta = sqrtf(1 - cosTheta*cosTheta);
		float dist = scLen*sinTheta;

		if (dist <= m_radius)
		{
			const float dd = sqrtf(m_radius * m_radius - dist * dist);
			Vec3 t0 = ray.m_direction * (scLen*cosTheta - dd);
			point = ray.m_origin + t0;

			normal = (point - m_center).Normal();

			return true;
		}
	}

	return false;
}


QuadShape::QuadShape(Vec3 center, Vec3 normal, Vec3 tangent, float width, float height)
	: m_normal(normal)
{
	Vec3 biNormal = normal.Cross(tangent);

	m_corners[0] = center - tangent * (width*0.5f) - biNormal * (height*0.5f);
	m_corners[1] = center + tangent * (width*0.5f) - biNormal * (height*0.5f);
	m_corners[2] = center + tangent * (width*0.5f) + biNormal * (height*0.5f);
	m_corners[3] = center - tangent * (width*0.5f) + biNormal * (height*0.5f);

	m_edgeNormal[0] = (m_corners[1] - m_corners[0]).Cross(m_normal);
	m_edgeNormal[1] = (m_corners[2] - m_corners[1]).Cross(m_normal);
	m_edgeNormal[2] = (m_corners[3] - m_corners[2]).Cross(m_normal);
	m_edgeNormal[3] = (m_corners[0] - m_corners[3]).Cross(m_normal);
}

bool QuadShape::Intersects(const Ray& ray, Vec3& point, Vec3& normal)
{
	float cosTheta = ray.m_direction.Dot(m_normal);
	Vec3 sa = m_corners[0] - ray.m_origin;
	float dist = sa.Dot(m_normal);

	if (fabs(cosTheta) > FLT_EPSILON && // ray is not parallel
		cosTheta * dist > 0) // ray is pointing plane
	{
		float rayDist = (dist / cosTheta);

		Vec3 si = ray.m_origin + ray.m_direction*rayDist;

		float r1 = (si - m_corners[0]).Dot(m_edgeNormal[0]);
		float r2 = (si - m_corners[1]).Dot(m_edgeNormal[1]);
		float r3 = (si - m_corners[2]).Dot(m_edgeNormal[2]);
		float r4 = (si - m_corners[3]).Dot(m_edgeNormal[3]);

		// on the same side of edges
		if (r1*r2 >= 0 && r1*r3 >= 0 && r1*r4 >= 0)
		{
			normal = m_normal * (dist > 0 ? -1.0f : 1.0f);
			point = si;
			return true;
		}
	}
	return false;
}