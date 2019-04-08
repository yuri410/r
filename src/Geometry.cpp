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
{
	Vec3 biNormal = normal.Cross(tangent.Normal());
	Vec3 corners[4];
	
	corners[0] = center - tangent * (width*0.5f) - biNormal * (height*0.5f);
	corners[1] = center + tangent * (width*0.5f) - biNormal * (height*0.5f);
	corners[2] = center + tangent * (width*0.5f) + biNormal * (height*0.5f);
	corners[3] = center - tangent * (width*0.5f) + biNormal * (height*0.5f);

	m_quad = PlanarConvexShape<4>(corners, normal);
}

bool QuadShape::Intersects(const Ray& ray, Vec3& point, Vec3& normal)
{
	return m_quad.Intersects(ray, point, normal);
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

//////////////////////////////////////////////////////////////////////////

TrimeshShape::TrimeshShape(const std::vector<Triangle>& triangles)
	: m_rootNode(triangles)
{

}

TrimeshShape::~TrimeshShape()
{

}

bool TrimeshShape::Intersects(const Ray& ray, Vec3& point, Vec3& normal)
{
	float nearestDist = FLT_MAX;
	return m_rootNode.Intersects(ray, point, normal, nearestDist);
}

TrimeshShape* TrimeshShape::GenerateNoiseQuad(Vec3 pos, double xSpan, double ySpan, int xSegments, int ySegments, const PerlinNoise& noise)
{
	double noiseScale = min(xSpan, ySpan)*0.05;

	std::vector<Triangle> triangles;

	for (int y = 0; y < ySegments; y++)
	{
		double yPos = ySpan * (((double)y / ySegments) * 2 - 1);
		double yPosN = ySpan * (((double)(y + 1) / ySegments) * 2 - 1);

		for (int x = 0; x < xSegments; x++)
		{
			double xPos = xSpan * (((double)x / xSegments) * 2 - 1);
			double xPosN = xSpan * (((double)(x + 1) / xSegments) * 2 - 1);

			double noiseTL = noise.GetValue2D(xPos, yPos) * noiseScale;
			double noiseTR = noise.GetValue2D(xPosN, yPos) * noiseScale;
			double noiseBL = noise.GetValue2D(xPos, yPosN) * noiseScale;
			double noiseBR = noise.GetValue2D(xPosN, yPosN) * noiseScale;

			Triangle triA =
			{
				Vec3((float)xPos,  (float)yPosN, (float)noiseBL) + pos,
				Vec3((float)xPosN, (float)yPos,  (float)noiseTR) + pos,
				Vec3((float)xPos,  (float)yPos,  (float)noiseTL) + pos,
			};

			Triangle triB =
			{
				Vec3((float)xPosN, (float)yPos,  (float)noiseTR) + pos,
				Vec3((float)xPos,  (float)yPosN, (float)noiseBL) + pos,
				Vec3((float)xPosN, (float)yPosN, (float)noiseBR) + pos,
			};

			triangles.push_back(triA);
			triangles.push_back(triB);
		}
	}

	return new TrimeshShape(triangles);
}


TrimeshShape::BvhNode::BvhNode(const std::vector<Triangle>& triangles)
{
	m_boundingSphere = BoundingSphere(triangles);

	if (triangles.size() > 8)
	{
		std::vector<Triangle> partitions[8];

		Vec3 center = m_boundingSphere.m_center;

		for (auto& tri : triangles)
		{
			Vec3 triCenter = tri.GetCenter();

			float x = triCenter.X - center.X;
			float y = triCenter.Y - center.Y;
			float z = triCenter.Z - center.Z;

			int partitionIndex = 0;

			if (y > 0)
			{
				if (x > 0)
				{
					partitionIndex = z > 0 ? 0 : 3;
				}
				else
				{
					partitionIndex = z > 0 ? 1 : 2;
				}
			}
			else
			{
				if (x > 0)
				{
					partitionIndex = z > 0 ? 4 : 7;
				}
				else
				{
					partitionIndex = z > 0 ? 5 : 6;
				}
			}

			partitions[partitionIndex].push_back(tri);
		}

		for (auto& part : partitions)
		{
			if (part.size())
				m_nodes.push_back(BvhNode(part));
		}
	}
	else
	{
		m_triangles = triangles;
	}
}

TrimeshShape::BvhNode::BvhNode(BvhNode&& o)
	: m_boundingSphere(o.m_boundingSphere), m_nodes(std::move(o.m_nodes)), m_triangles(std::move(o.m_triangles))
{ }

TrimeshShape::BvhNode& TrimeshShape::BvhNode::operator=(BvhNode&& o)
{
	if (this == &o)
	{
		this->~BvhNode();
		new (this)BvhNode(std::move(o));
	}
	return *this;
}

bool TrimeshShape::BvhNode::Intersects(const Ray& ray, Vec3& point, Vec3& normal, float& nearestDist) const
{
	if (m_boundingSphere.IntersectsRay(ray))
	{
		bool ret = false;
		for (auto& node : m_nodes)
		{
			ret |= node.Intersects(ray, point, normal, nearestDist);
		}

		for (auto& tri : m_triangles)
		{
			float d;
			Vec3 pt, n;
			
			if (tri.Intersects(ray, pt, n) && 
				(d = (pt - ray.m_origin).LengthSqr()) < nearestDist)
			{
				nearestDist = d;
				point = pt;
				normal = n;
				ret = true;
			}
		}
		return ret;
	}
	return false;
}