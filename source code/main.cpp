#define CVUI_IMPLEMENTATION
#include "cvui.h"
#include "prng.h"
#include "noise.h"
#include <ctime>
#include <fstream>
#include <iostream>
using namespace cv;


#define WINDOW_NAME "GABOR_NOISE"

unsigned char uniform(float intensity) {
	if (intensity <= 0.0) {
		return 0;
	}
	else if (intensity < 1.0) {
		return static_cast<unsigned char>(intensity * 255);
	}
	else {
		return 255;
	}
}


uchar* process(unsigned resolution, float K_, float a_, float F_0_, float omega_0_, bool isotropic, float number_of_impulses_per_kernel, unsigned period, unsigned random_offset) {
	//calculate the noise in spacial and frequency domain
	Noise noise_(K_, a_, F_0_, omega_0_, isotropic, number_of_impulses_per_kernel, period, random_offset);
	float* image_val = new float[resolution * resolution];
	float* image_fre = new float[resolution * resolution];
	float scale = 3.0 * std::sqrt(noise_.variance());
	for (unsigned row = 0; row < resolution; ++row) {
		for (unsigned col = 0; col < resolution; ++col) {
			float x = (float(col) + 0.5) - (float(resolution) / 2.0);
			float y = (float(resolution - row - 1) + 0.5) - (float(resolution) / 2.0);
			Noise_com val = noise_.calculate(x, y);
			image_val[(row * resolution) + col] = 0.5 + (0.5 * (val.noise_val / scale));
			image_fre[(row * resolution) + col] = val.noise_fre;
		}
	}

	//uniform the picture to 256 scale gray picture
	uchar* temp_val = new uchar[resolution * resolution];
	uchar* temp_fre = new uchar[resolution * resolution];
	for (unsigned row = 0; row < resolution; ++row) {
		for (unsigned col = 0; col < resolution; ++col) {
			temp_val[(row * resolution) + col] = uniform(image_val[(row * resolution) + col]);
			temp_fre[(row * resolution) + col] = uniform(image_fre[(row * resolution) + col]);
		}
	}

	//calculate a single kernel in frequency domain
	uchar* kernel_fre = new uchar[resolution * resolution];
	for (unsigned row = 0; row < resolution; ++row) {
		for (unsigned col = 0; col < resolution; ++col) {
			float x = (float(col) + 0.5) - (float(resolution) / 2.0);
			float y = (float(resolution - row - 1) + 0.5) - (float(resolution) / 2.0);
			x = x / (float)resolution;
			y = y / (float)resolution;
			float f_cos = F_0_ * std::cos(omega_0_);
			float f_sin = F_0_ * std::sin(omega_0_);

			float part1 = std::exp((-M_PI * (pow((x - f_cos), 2) + pow((y - f_sin), 2))) / (a_ * a_));
			float part2 = std::exp((-M_PI * (pow((x + f_cos), 2) + pow((y + f_sin), 2))) / (a_ * a_));
			float rst = K_ * (part1 + part2) / (2 * a_ * a_);
			kernel_fre[row * resolution + col] = uniform(rst);
		}
	}

	//fill 3 pictures in one window
	uchar* window = new uchar[resolution * resolution * 4];
	for (unsigned i = 0; i < resolution; ++i) {
		for (unsigned j = 0; j < resolution; ++j) {
			window[(i * 2 * resolution) + j] = 255;
			window[(i * 2 * resolution) + resolution + j] = temp_val[i * resolution + j];
			window[2 * resolution * (resolution + i) + j] = kernel_fre[i * resolution + j];
			window[2 * resolution * (resolution + i) + resolution + j] = temp_fre[i * resolution + j];
		}
	}

	return window;
}

int main() {

	unsigned resolution = 256;
	float K_ = 4.0;
	float a_ = 0.05;
	float F_0_ = 0.2;
	float omega_0_ = M_PI / 4.0;
	float number_of_impulses_per_kernel = 64.0;
	unsigned period = 256;
	unsigned random_offset = std::time(0);

	bool is = false;
	bool isotropic = false;

	// Init a OpenCV window and tell cvui to use it.
	namedWindow(WINDOW_NAME);
	cvui::init(WINDOW_NAME);
	
	uchar* temp = process(resolution, K_, a_, F_0_, omega_0_, isotropic, number_of_impulses_per_kernel, period, random_offset);
	Mat frame = Mat(512, 512, CV_8UC1, temp);
	
	while (true) {
				
		if (is) {
			temp = process(resolution, K_, a_, F_0_, omega_0_, isotropic, number_of_impulses_per_kernel, period, random_offset);
			frame = Mat(512, 512, CV_8UC1, temp);
			is = false;
		}


		// Render the settings window to house the UI
		cvui::window(frame, 0, 0, 256, 256, "Settings");
		
		//Button to generate
		is = cvui::button(frame, 140, 25, "Generate");
		
		// Checkbox to choose isotropic or anisotropic
		cvui::checkbox(frame, 30, 30, "isotropic", &isotropic);

		// four trackbars to control the paramters of gabor noise
		cvui::text(frame, 10, 65, "K");
		cvui::trackbar(frame, 60, 50, 165, &K_, 0.5f, 5.0f);
		cvui::text(frame, 10, 110, "a");
		cvui::trackbar(frame, 60, 95, 165, &a_, 0.005f, 0.1f, 1, "%.3Lf");
		cvui::text(frame, 10, 155, "F");
		cvui::trackbar(frame, 60, 140, 165, &F_0_, 0.01f, 0.3f, 1, "%.4Lf");
		cvui::text(frame, 10, 200, "omega");
		cvui::trackbar(frame, 60, 185, 165, &omega_0_, 0.0f, (float)M_PI, 1, "%.2Lf");

		// Update cvui internal stuff
		cvui::update();
		imshow(WINDOW_NAME, frame);

		if (waitKey(30) == 27) {
			break;
		}
	}
	
	return 0;

}

