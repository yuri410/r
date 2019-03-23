#pragma once
#include "Common.h"

struct HemisphereCoord;

const float PI = 3.14159f;

inline float Sqr(float v) { return v*v; }
inline float Saturate(float v) { return v < 0.0f ? 0.0f : (v > 1.0f ? 1.0f : v); }

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
