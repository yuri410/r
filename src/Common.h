#pragma once

#include <GL/OOGL.hpp>
#include <cassert>

using Vec3 = GL::Vec3;

inline Vec3 operator*(const Vec3& a, const Vec3& b) { return Vec3(a.X*b.X, a.Y*b.Y, a.Z*b.Z); }

class IGeometry;
struct Material;
class Scene;
class SceneObject;
struct Camera;
struct Ray;

class Random;

typedef unsigned int uint;