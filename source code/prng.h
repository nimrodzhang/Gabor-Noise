#ifndef _PRNG_H_
#define _PRNG_H_

#include <climits>
#include <cmath>

#define M_PI 3.14159265358979323846

class pseudo_random_number_generator
{
public:
	void seed(unsigned s) { x_ = s; }
	unsigned operator()() { x_ *= 3039177861u; return x_; }
	float uniform_0_1() { return float(operator()()) / float(UINT_MAX); }
	float uniform(float min, float max)
	{
		return min + (uniform_0_1() * (max - min));
	}
	unsigned poisson(float mean)
	{
		float g_ = std::exp(-mean);
		unsigned em = 0;
		double t = uniform_0_1();
		while (t > g_) {
			++em;
			t *= uniform_0_1();
		}
		return em;
	}
private:
	unsigned x_;
};

#endif