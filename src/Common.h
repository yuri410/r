#pragma once

#include <GL/OOGL.hpp>
#include <cassert>

using Vec3 = GL::Vec3;

struct Ray
{
	Vec3 m_origin;
	Vec3 m_direction;
};

inline void muAutoFrameZ(const Vec3& dz, Vec3& dx, Vec3& dy)
{
	Vec3 refVector = Vec3(0, 1, 0);
	if (fabs(refVector.Dot(dz)) > 0.99f)
	{
		refVector = Vec3(-1, 0, 0);
	}
	dx = dz.Cross(refVector);
	dx = dx.Normal();

	dy = dz.Cross(dx);
}

inline Vec3 operator*(const Vec3& a, const Vec3& b) { return Vec3(a.X*b.X, a.Y*b.Y, a.Z*b.Z); }

class IGeometry;
struct Material;
class Scene;
class SceneObject;
struct Camera;

class Random;

typedef unsigned int uint;