#pragma once

#include <GL/OOGL.hpp>
#include <cassert>

using Vec3 = GL::Vec3;

struct Ray
{
	Vec3 m_origin;
	Vec3 m_direction;
};


class IGeometry;
class Material;
class Scene;
class SceneObject;
struct Camera;

typedef unsigned int uint;