#ifndef _NOISE_H_
#define _NOISE_H_

#include "prng.h"

typedef struct Noise_com {
	float noise_val;
	float noise_fre;
};


class Noise {
private:
	float K_;
	float a_;
	float F_0_;
	float omega_0_;
	bool isotropic_;
	float kernel_radius_;
	float impulse_density_;
	unsigned period_;
	unsigned random_offset_;

public:
	Noise(float K, float a, float F_0, float omega_0, bool isotropic, float number_of_impulses_per_kernel, unsigned period, unsigned random_offset)
		: K_(K), a_(a), F_0_(F_0), omega_0_(omega_0), isotropic_(isotropic), period_(period), random_offset_(random_offset)
	{
		kernel_radius_ = std::sqrt(-std::log(0.05) / M_PI) / a_;		//�˰뾶
		impulse_density_ = number_of_impulses_per_kernel / (M_PI * kernel_radius_ * kernel_radius_);
	}

	float gabor(float K, float a, float F_0, float omega_0, float x, float y);

	float gabor_fre(float K, float a, float F_0, float omega_0, float x, float y);

	unsigned morton(unsigned x, unsigned y);

	Noise_com calculate(float x, float y);

	Noise_com cell(int i, int j, float x, float y);

	float variance();
};


#endif