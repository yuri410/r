#pragma once

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
