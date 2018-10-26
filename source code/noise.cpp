#include "noise.h"


unsigned Noise::morton(unsigned x, unsigned y) {
	unsigned z = 0;
	for (unsigned i = 0; i < (sizeof(unsigned) * CHAR_BIT); ++i) {
		z |= ((x & (1 << i)) << i) | ((y & (1 << i)) << (i + 1));
	}
	return z;
}

//Gabor Kernel
float Noise::gabor(float K, float a, float F_0, float omega_0, float x, float y) {
	float gaussian_envelop = K * std::exp(-M_PI * (a * a) * ((x * x) + (y * y)));
	float sinusoidal_carrier = std::cos(2.0 * M_PI * F_0 * ((x * std::cos(omega_0)) + (y * std::sin(omega_0))));
	return gaussian_envelop * sinusoidal_carrier;
}


//Gabor Kernel in frequency domain
float Noise::gabor_fre(float K, float a, float F_0, float omega_0, float x, float y) {
	float f_cos = F_0 * std::cos(omega_0);
	float f_sin = F_0 * std::sin(omega_0);

	float part1 = std::exp((-M_PI * (pow((x - f_cos), 2) + pow((y - f_sin), 2))) / (a * a));
	float part2 = std::exp((-M_PI * (pow((x + f_cos), 2) + pow((y + f_sin), 2))) / (a * a));
	float rst = K * (part1 + part2) / (2 * a * a);
	return rst;
}



Noise_com Noise::calculate(float x, float y) {
	x /= kernel_radius_, y /= kernel_radius_;
	float int_x = std::floor(x), int_y = std::floor(y);
	float frac_x = x - int_x, frac_y = y - int_y;
	int i = int(int_x), j = int(int_y);
	float noise_v = 0.0;
	float noise_f = 0.0;

	for (int di = -1; di <= +1; ++di) {
		for (int dj = -1; dj <= +1; ++dj) {
			Noise_com val = cell(i + di, j + dj, frac_x - di, frac_y - dj);
			noise_v += val.noise_val;
			noise_f += val.noise_fre;
		}
	}
	Noise_com rst;
	rst.noise_val = noise_v;
	rst.noise_fre = noise_f;
	return rst;
}

Noise_com Noise::cell(int i, int j, float x, float y) {
	float int_x = (float)i;
	float int_y = (float)j;
	unsigned s = morton(i, j) + random_offset_;
	if (s == 0) 
		s = 1;
	pseudo_random_number_generator prng;
	prng.seed(s);
	double number_of_impulses_per_cell = impulse_density_ * kernel_radius_ * kernel_radius_;
	unsigned number_of_impulses = prng.poisson(number_of_impulses_per_cell);
	float noise_v = 0.0;
	float noise_f = 0.0;
	for (unsigned i = 0; i < number_of_impulses; ++i) {
		float x_i = prng.uniform_0_1();
		float y_i = prng.uniform_0_1();
		float w_i = prng.uniform(-1.0, +1.0);
		float z_i = (w_i + 1.0) / 2.0;

		float x_i_x = x - x_i;
		float y_i_y = y - y_i;
		if (((x_i_x * x_i_x) + (y_i_y * y_i_y)) < 1.0) {
			if (isotropic_) {
				float omega_0_i = prng.uniform(0.0, 2.0 * M_PI);
				noise_v += w_i * gabor(K_, a_, F_0_, omega_0_i, x_i_x * kernel_radius_, y_i_y * kernel_radius_); // isotropic
				noise_f += w_i * w_i * gabor_fre(K_, a_, F_0_, omega_0_i, (int_x + x)*kernel_radius_ / 256, (int_y + y)*kernel_radius_ / 256);
			}
			else {
				noise_v += w_i * gabor(K_, a_, F_0_, omega_0_, x_i_x * kernel_radius_, y_i_y * kernel_radius_); // anisotropic
				noise_f += w_i * w_i * gabor_fre(K_, a_, F_0_, omega_0_, (int_x + x)*kernel_radius_ / 256, (int_y + y)*kernel_radius_ / 256);
			}
		}
	}
	Noise_com rst;
	rst.noise_val = noise_v;
	rst.noise_fre = noise_f;
	return rst;
}

float Noise::variance() {
	float integral_gabor_filter_squared = ((K_ * K_) / (4.0 * a_ * a_)) * (1.0 + std::exp(-(2.0 * M_PI * F_0_ * F_0_) / (a_ * a_)));
	return impulse_density_ * (1.0 / 3.0) * integral_gabor_filter_squared;
}