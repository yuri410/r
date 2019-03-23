#include "Geometry.h"
#include "Math.h"

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
		float sinTheta = sqrt(1 - cosTheta*cosTheta);
		float dist = scLen*sinTheta;

		if (dist <= m_radius)
		{
			const float dd = sqrt(m_radius * m_radius - dist * dist);
			Vec3 t0 = ray.m_direction * (scLen*cosTheta - dd);
			point = ray.m_origin + t0;

			normal = (point - m_center).Normal();

			return true;
		}
	}

	return false;
}

//////////////////////////////////////////////////////////////////////////

QuadShape::QuadShape(Vec3 center, Vec3 normal, Vec3 tangent, float width, float height)
	: m_normal(normal.Normal())
{
	Vec3 biNormal = normal.Cross(tangent.Normal());

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

	if (abs(cosTheta) < 1 && // ray is not parallel
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

//////////////////////////////////////////////////////////////////////////

DiskShape::DiskShape(Vec3 center, Vec3 normal, float radius)
	: m_center(center), m_normal(normal.Normal()), m_radius(radius)
{

}

bool DiskShape::Intersects(const Ray& ray, Vec3& point, Vec3& normal)
{
	float cosTheta = ray.m_direction.Dot(m_normal);
	Vec3 sa = m_center - ray.m_origin;
	float dist = sa.Dot(m_normal);

	if (abs(cosTheta) < 1 && // ray is not parallel
		cosTheta * dist > 0) // ray is pointing plane
	{
		float rayDist = (dist / cosTheta);

		Vec3 si = ray.m_origin + ray.m_direction*rayDist;

		// inside the disk?
		float d = (si - m_center).Length();
		if (d < m_radius)
		{
			normal = m_normal * (dist > 0 ? -1.0f : 1.0f);
			point = si;
			return true;
		}
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////

CylinderShape::CylinderShape(Vec3 center, Vec3 axis, float height, float radius)
	: m_center(center), m_axis(axis.Normal()), m_height(height), m_radius(radius)
{

}

bool CylinderShape::Intersects(const Ray& ray, Vec3& point, Vec3& normal)
{
	// Not in parallel with body
	if (0 && abs(ray.m_direction.Dot(m_axis) < 1))
	{
		Vec3 cn = ray.m_direction.Cross(m_axis).Normal();
		Vec3 sc = m_center - ray.m_origin;
		float dist = sc.Dot(cn);
		float scLen = sc.Length();
		sc = sc * (1.0f / scLen);
		float cosTheta = ray.m_direction.Dot(sc);

		if (abs(dist) <= m_radius)
		{
			//float a = 1;
			//float b = ray.m_direction.Dot(m_axis);
			//float c = 1;
			//float d = ray.m_direction.Dot(sc);
			//float e = m_axis.Dot(sc);

			//float s = (b*e - c*d) / (a*c - b*b);
			//float t = (a*e - b*d) / (a*c - b*b);

			//Vec3 discCenter = m_center + m_axis * t;
			//Vec3 si = ray.m_origin + ray.m_direction*s;


			float cosPhi = abs(ray.m_direction.Dot(m_axis));
			const float dd = sqrt(m_radius * m_radius - dist * dist);
			Vec3 t0 = ray.m_direction * (scLen*cosTheta - dd / cosPhi);
			point = ray.m_origin + t0;

			// Is in height range?
			if (abs((point - m_center).Dot(m_axis)) <= m_height*0.5f)
			{
				normal = m_axis.Cross(point - m_center).Cross(m_axis);
				normal = normal.Normal();

				//normal = normal.Normal() * -1;
				//normal = normal.Normal() * -1;
				//normal = (point - m_center).Normal();
				return true;
			}
		}
	}

	bool result = false;
	float minDist = FLT_MAX;

	Vec3 p, n;
	if (result |= IntersectsCap(ray, m_center + m_axis*(m_height*0.5f), m_axis, p, n))
	{
		float d = (p - ray.m_origin).LengthSqr();
		if (d < minDist)
		{
			minDist = d;
			point = p;
			normal = n;
		}
	}
	if (result |= IntersectsCap(ray, m_center - m_axis*(m_height*0.5f), m_axis*-1, p, n))
	{
		float d = (p - ray.m_origin).LengthSqr();
		if (d < minDist)
		{
			minDist = d;
			point = p;
			normal = n;
		}
	}
	if (result |= IntersectsBody(ray, p, n))
	{
		float d = (p - ray.m_origin).LengthSqr();
		if (d < minDist)
		{
			minDist = d;
			point = p;
			normal = n;
		}
	}

	return result;
}

bool CylinderShape::IntersectsCap(const Ray& ray, const Vec3& capCenter, const Vec3& capNormal, Vec3& point, Vec3& normal) const
{
	Vec3 dn = capNormal;
	float cosTheta = ray.m_direction.Dot(dn);
	Vec3 sa = capCenter - ray.m_origin;
	float dist = sa.Dot(dn);

	if (abs(cosTheta) < 1 && // ray is not parallel
		cosTheta * dist > 0) // ray is pointing plane
	{
		float rayDist = (dist / cosTheta);

		Vec3 si = ray.m_origin + ray.m_direction*rayDist;

		// inside the disk?
		float d = (si - capCenter).Length();
		if (d < m_radius)
		{
			normal = capNormal * (dist > 0 ? -1.0f : 1.0f);
			point = si;
			return true;
		}
	}
	return false;
}

bool CylinderShape::IntersectsBody(const Ray& ray, Vec3& point, Vec3& normal) const
{
	Vec3 sc = ray.m_origin - m_center;
	Vec3 m = sc.Cross(m_axis);
	Vec3 cn = ray.m_direction.Cross(m_axis);

	float a = cn.LengthSqr();
	float b = 2 * cn.Dot(m);
	float c = m.LengthSqr() - (m_radius*m_radius);
	float d = b * b - 4 * a * c;
	if (d < 0)
		return false;

	float time = (-b - sqrt(d)) / (2 * a);
	if (time < 0)
		return false;

	point = ray.m_origin + ray.m_direction * time;

	// Is in height range?
	if (abs((point - m_center).Dot(m_axis)) <= m_height*0.5f)
	{
		normal = m_axis.Cross(point - m_center).Cross(m_axis);
		return true;
	}
	return false;
}