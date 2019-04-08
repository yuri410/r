#include "Math.h"
#include <ctime>


void ImportanceSampleUniform(float e1, float e2, HemisphereCoord& c, float& pdf)
{
	float phi = e1 * 2 * PI;

	c.m_cosTheta = e2;
	c.m_sinTheta = sqrt(1 - c.m_cosTheta*c.m_cosTheta);

	c.m_cosPhi = cos(phi);
	c.m_sinPhi = sin(phi);

	pdf = 1.0f / (2 * PI);
}

void ImportanceSampleCosine(float e1, float e2, HemisphereCoord& c, float& pdf)
{
	float phi = e1 * 2 * PI;

	c.m_cosTheta = sqrt(1 - e2);
	c.m_sinTheta = sqrt(e2);

	c.m_cosPhi = cos(phi);
	c.m_sinPhi = sin(phi);

	pdf = c.m_cosTheta / PI;
}

void ImportanceSampleGGX(float e1, float e2, float a, HemisphereCoord& c, float& pdf)
{
	float phi = e1 * 2 * PI;

	float a2 = a*a;
	float rc = (1 - e2) / ((a2 - 1.0f) *e2 + 1);
	rc = Saturate(rc);

	c.m_cosTheta = sqrt(rc);
	c.m_sinTheta = sqrt(1 - rc);

	c.m_cosPhi = cos(phi);
	c.m_sinPhi = sin(phi);

	//pdf = a2 * c.m_cosTheta / (PI * Sqr((a2 - 1)*c.m_cosTheta*c.m_cosTheta + 1));
	pdf = a2 * c.m_cosTheta / (PI * Sqr((a2 - 1)*rc + 1));
}

void AutoFrameZ(const Vec3& dz, Vec3& dx, Vec3& dy)
{
	Vec3 refVector = Vec3(0, 1, 0);
	if (abs(refVector.Dot(dz)) > 0.99f)
	{
		refVector = Vec3(-1, 0, 0);
	}
	dx = dz.Cross(refVector);
	dx = dx.Normal();

	dy = dz.Cross(dx);
}

//////////////////////////////////////////////////////////////////////////

Vec3 HemisphereCoord::GetCartesianCoord_AutoFrameZ(Vec3 basisZ) const
{
	Vec3 basisX, basisY;
	AutoFrameZ(basisZ, basisX, basisY);

	Vec3 xi = basisX*(m_cosPhi*m_sinTheta) + basisY*(m_sinPhi*m_sinTheta) + basisZ * m_cosTheta;
	return xi.Normal();
}

//////////////////////////////////////////////////////////////////////////

BoundingSphere::BoundingSphere(const Vec3& center, float radius)
	: m_center(center), m_radius(radius) { }

BoundingSphere::BoundingSphere(const std::vector<Triangle>& triangles)
	: m_center({ 0,0,0 }), m_radius(0.0f)
{
	Vec3 center = { 0,0,0 };
	int vertexCount = 0;
	for (auto& tri : triangles)
	{
		for (auto& vtx : tri.m_corners)
		{
			center += vtx;
			vertexCount++;
		}
	}

	if (vertexCount > 0)
	{
		center = center * (1.0f / vertexCount);

		float r2 = 0;
		for (auto& tri : triangles)
		{
			for (auto& vtx : tri.m_corners)
			{
				float distSq = (vtx - center).LengthSqr();
				if (distSq > r2)
					r2 = distSq;
			}
		}

		m_center = center;
		m_radius = sqrtf(r2);
	}
}

bool BoundingSphere::IntersectsRay(const Ray& ray) const
{
	Vec3 sc = m_center - ray.m_origin;
	float slen = ray.m_direction.Dot(sc);

	if (slen > 0)
	{
		float distSq = sc.LengthSqr() - slen * slen;

		return distSq <= m_radius*m_radius;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////

Triangle::Triangle(const Vec3& a, const Vec3& b, const Vec3& c)
	: PlanarConvexShape({ a,b,c }, (c - a).Cross(b - a).Normal())
{

}

//////////////////////////////////////////////////////////////////////////

Random::Random()
{
	const int seed = static_cast<int>(time(0));
	SetSeed(seed, true);
}

void Random::SetSeed(int seed, bool reset)
{
	if (reset)
	{
		Reset();
	}

	m_seed = seed;
	int holdRand = seed;

	for (uint& s : m_state)
	{
		s = (((holdRand = holdRand * 214013L + 2531011L) >> 16) & 0x7fff);
	}
}

void Random::Reset()
{
	memset(m_state, 0, sizeof(m_state));
	m_seed = 0;
	m_index = 0;
}

int Random::RawSample()
{
	// WELLRNG512
	uint a, b, c, d;
	a = m_state[m_index];
	c = m_state[(m_index + 13) & 15];
	b = a^c ^ (a << 16) ^ (c << 15);
	c = m_state[(m_index + 9) & 15];
	c ^= (c >> 11);
	a = m_state[m_index] = b^c;
	d = a ^ ((a << 5) & 0xDA442D20UL);
	m_index = (m_index + 15) & 15;
	a = m_state[m_index];
	m_state[m_index] = a^b^d ^ (a << 2) ^ (b << 18) ^ (c << 28);
	return static_cast<int>(m_state[m_index] & 0x7fffffffU);
}

//////////////////////////////////////////////////////////////////////////

PerlinNoise::PerlinNoise(double frequency, double persistence, int octaves, int seed)
	: m_frequency(frequency), m_persistence(persistence), m_octaves(octaves), m_seed(seed)
{

}

double PerlinNoise::GetValue2D(double x, double y) const
{
	double freq = m_frequency;
	double persistence = m_persistence;
	int octaves = m_octaves;
	int seed = m_seed;

	double t = 0;
	double amplitude = 1;

	for (int k = 0; k < octaves; k++)
	{
		t += InterpolatedNoise2D(x * freq + seed, y * freq + seed) * amplitude;
		amplitude *= persistence;
		freq *= 2;
	}

	return t;
}

double PerlinNoise::InterpolatedNoise2D(double x, double y) const
{
	int Xint = x >= 0 ? (int)x : (int)x - 1;
	int Yint = y >= 0 ? (int)y : (int)y - 1;

	double Xfrac = x - Xint;
	double Yfrac = y - Yint;

	// pre-modulation to prevent duplicated calculation
	ModulateLerpAmount(Xfrac);
	ModulateLerpAmount(Yfrac);

	// sampling original: 
	//  x0y0 = Sample(x, y)
	//  x1y0 = Sample(x+1, y) 
	//  x0y1 = Sample(x, y+1)
	//  x1y1 = Sample(x+1, y+1)
	//  
	//  where Sample() is:
	//     corners = ([x-1, y-1] + [x+1, y-1] + [x-1, y+1] + [x+1, y+1]) / 16.0;
	//     sides =   ([x-1, y] + [x+1, y] + [x, y-1] + [x, y+1]) / 8.0;
	//     center =   [x, y] / 4.0;
	//     return corners + sides + center;

	// expanded to reduce duplicated sampling
	// generates less instructions
	double n01 = Noise2D(Xint - 1, Yint - 1);
	double n02 = Noise2D(Xint + 1, Yint - 1);
	double n03 = Noise2D(Xint - 1, Yint + 1);
	double n04 = Noise2D(Xint + 1, Yint + 1);
	double n05 = Noise2D(Xint - 1, Yint);
	double n06 = Noise2D(Xint + 1, Yint);
	double n07 = Noise2D(Xint, Yint - 1);
	double n08 = Noise2D(Xint, Yint + 1);
	double n09 = Noise2D(Xint, Yint);

	double n12 = Noise2D(Xint + 2, Yint - 1);
	double n14 = Noise2D(Xint + 2, Yint + 1);
	double n16 = Noise2D(Xint + 2, Yint);

	double n23 = Noise2D(Xint - 1, Yint + 2);
	double n24 = Noise2D(Xint + 1, Yint + 2);
	double n28 = Noise2D(Xint, Yint + 2);

	double n34 = Noise2D(Xint + 2, Yint + 2);

	double x0y0 = 0.0625*(n01 + n02 + n03 + n04) + 0.125*(n05 + n06 + n07 + n08) + 0.25*(n09);
	double x1y0 = 0.0625*(n07 + n12 + n08 + n14) + 0.125*(n09 + n16 + n02 + n04) + 0.25*(n06);
	double x0y1 = 0.0625*(n05 + n06 + n23 + n24) + 0.125*(n03 + n04 + n09 + n28) + 0.25*(n08);
	double x1y1 = 0.0625*(n09 + n16 + n28 + n34) + 0.125*(n08 + n14 + n06 + n24) + 0.25*(n04);

	// bilinear
	double v1 = LinearInterpolate(x0y0, x1y0, Xfrac);
	double v2 = LinearInterpolate(x0y1, x1y1, Xfrac);
	return LinearInterpolate(v1, v2, Yfrac);
}

double PerlinNoise::SampleNoise2D(int x, int y) const
{
	double corners = Noise2D(x - 1, y - 1) + Noise2D(x + 1, y - 1) + Noise2D(x - 1, y + 1) + Noise2D(x + 1, y + 1);
	double sides = Noise2D(x - 1, y) + Noise2D(x + 1, y) + Noise2D(x, y - 1) + Noise2D(x, y + 1);
	double center = Noise2D(x, y);
	return corners * (0.25 / 4.0) + sides * (0.5 / 4.0) + center * 0.25;
}

double PerlinNoise::Noise2D(int x, int y) const
{
	int n = x + y * 3251;
	n = (n << 13) ^ n;
	int t = (n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff;
	//return ( 1.0 - ( (n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0f);    

	return 1.0 - double(t) * 0.931322574615478515625e-9;/// 1073741824.0);
}

void PerlinNoise::ModulateLerpAmount(double& a)
{
	// map to cubic S-curve
	a = a * a * (3.0 - 2.0 * a);
}