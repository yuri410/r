#include "Common.h"
#include "Random.h"
#include <ctime>

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