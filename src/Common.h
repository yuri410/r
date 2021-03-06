#pragma once

#include <GL/OOGL.hpp>
#include <cassert>

using Vec3 = GL::Vec3;
using Mat4 = GL::Mat4;

inline Vec3 operator*(const Vec3& a, const Vec3& b) { return Vec3(a.X*b.X, a.Y*b.Y, a.Z*b.Z); }

class IGeometry;
struct Material;
class Scene;
class SceneObject;
struct Camera;
struct Ray;
struct Triangle;

class Random;
class PerlinNoise;

typedef unsigned int uint;