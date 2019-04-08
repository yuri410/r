#pragma once
#include "Common.h"

struct HemisphereCoord;

const float PI = 3.14159f;

inline float Sqr(float v) { return v*v; }
inline float Saturate(float v) { return v < 0.0f ? 0.0f : (v > 1.0f ? 1.0f : v); }
inline double LinearInterpolate(double x, double y, double a) { return x * (1 - a) + y * a; }

void ImportanceSampleUniform(float e1, float e2, HemisphereCoord& c, float& pdf);

void ImportanceSampleCosine(float e1, float e2, HemisphereCoord& c, float& pdf);

void ImportanceSampleGGX(float e1, float e2, float a, HemisphereCoord& c, float& pdf);

void AutoFrameZ(const Vec3& dz, Vec3& dx, Vec3& dy);

struct HemisphereCoord
{
	float m_cosPhi;
	float m_sinPhi;
	float m_cosTheta;
	float m_sinTheta;

	Vec3 GetCartesianCoord_AutoFrameZ(Vec3 z) const;
};

struct Ray
{
	Vec3 m_origin;
	Vec3 m_direction;
};

struct BoundingSphere
{
	Vec3  m_center;
	float m_radius;

	BoundingSphere() { }
	BoundingSphere(const Vec3& center, float radius);
	BoundingSphere(const std::vector<Triangle>& triangles);

	bool IntersectsRay(const Ray& ray) const;
};

template <int N>
struct PlanarConvexShape
{
	static_assert(N >= 3, "Needs at least 3 vertices.");

	Vec3  m_corners[N];
	Vec3  m_normal;
	Vec3  m_edgeNormal[N];

	PlanarConvexShape() { }
	PlanarConvexShape(const Vec3(&corners)[N], const Vec3& normal);

	bool Intersects(const Ray& ray, Vec3& point, Vec3& normal) const;
	Vec3 GetCenter() const;
};

struct Triangle : public PlanarConvexShape<3>
{
	Triangle() { }
	Triangle(const Vec3& a, const Vec3& b, const Vec3& c);
};

class Random
{
public:
	Random();
	Random(int seed) { SetSeed(seed, true); }
	~Random() { }

	int  GetSeed() const { return m_seed; }
	void SetSeed(int seed, bool reset);
	void Reset();

	int IntSample() { return RawSample(); }
	float FloatSample() { return RawSample() * (1.0f / 2147483647.0f); }

private:

	int RawSample();

	uint m_state[16];
	int m_index = 0;
	int m_seed = 0;
};

class PerlinNoise
{
public:
	PerlinNoise(double frequency, double persistence, int octaves, int seed);

	double GetValue2D(double x, double y) const;

private:
	double Noise2D(int x, int y) const;
	double InterpolatedNoise2D(double x, double y) const;
	double SampleNoise2D(int x, int y) const;
	
	static void ModulateLerpAmount(double& a);

	double m_frequency;
	double m_persistence;
	int    m_octaves;
	int    m_seed;
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <int N>
PlanarConvexShape<N>::PlanarConvexShape(const Vec3(&corners)[N], const Vec3& normal)
	: m_normal(normal)
{
	for (int i = 0; i < N; i++)
		m_corners[i] = corners[i];

	for (int i = 0; i < N; i++)
		m_edgeNormal[i] = (m_corners[(i + 1) % N] - m_corners[i]).Cross(m_normal);
}

template <int N>
bool PlanarConvexShape<N>::Intersects(const Ray& ray, Vec3& point, Vec3& normal) const
{
	float cosTheta = ray.m_direction.Dot(m_normal);
	Vec3 sa = m_corners[0] - ray.m_origin;
	float dist = sa.Dot(m_normal);

	if (abs(cosTheta) < 1 && // ray is not parallel
		cosTheta * dist > 0) // ray is pointing plane
	{
		float rayDist = (dist / cosTheta);

		Vec3 si = ray.m_origin + ray.m_direction*rayDist;

		float r[N];
		for (int i = 0; i < N; i++)
		{
			r[i] = (si - m_corners[i]).Dot(m_edgeNormal[i]);
		}

		float r0 = r[0];

		bool sameSide = true;
		for (int i = 1; i < N; i++)
		{
			// on the same side of edges
			sameSide &= (r0 * r[i] >= 0);
		}

		if (sameSide)
		{
			normal = m_normal * (dist > 0 ? -1.0f : 1.0f);
			point = si;
			return true;
		}
	}
	return false;
}

template <int N>
Vec3 PlanarConvexShape<N>::GetCenter() const
{
	Vec3 c = { 0,0,0 };
	for (auto& v : m_corners)
	{
		c += v;
	}
	return c*(1.0f / N);
}