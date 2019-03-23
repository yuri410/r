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

	c.m_sinTheta = sqrt(e2);
	c.m_cosTheta = sqrt(1 - e2);

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

	pdf = a2 * c.m_cosTheta / (PI * Sqr((a2 - 1)*c.m_cosTheta*c.m_cosTheta + 1));
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