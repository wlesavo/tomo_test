#include "CImg.h"
#include "Projector.hpp"
#include <iostream>
#include <fstream>
#include <string>

int main()
{   
    cimg_library::CImg<unsigned char> phantom("phantom.bmp"), target("sinogram.bmp");

    int size_x = phantom.width(), size_y = phantom.height();
    int detector_size = target.height(), angle_count = target.width();

    std::unique_ptr<unsigned char[]> img(new unsigned char[size_x * size_y]{});

    // copy first channel of image to array assuming grayscale bmp
    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            img[i + size_x * j] = (int)phantom(i, j, 0, 0);
        }
    }

    std::vector<float> angles;
    for (int i = 0; i < angle_count; ++i) {
        angles.push_back((float)(i) * cimg_library_suffixed::cimg::PI / 180.0f);
    }

    Geometry geom(angles, detector_size, 1.0f, 7000.0f, 500.0f);
    Projector myProj(std::move(img), size_x, size_y, geom);

    std::unique_ptr<unsigned char[]> fp = myProj.getFullProjectionImage();

    // compare the differences in image form,
    // count non-zero differences in pixels for statistics
    float delta = 0.0f, summ = 0.0f;
    int count = 0;
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            int a = (int)target(i, j, 0, 0);
            int b = (int)fp[i * geom.detectorCount + (j)];
            delta += std::abs(b - a);
            if (a!=b) {
                count += 1;
            }
        }
    }
    if (count == 0)
        std::cout << "Perfect!" << std::endl;
    else
        std::cout << "total delta: "<< delta << " count: " << count << " average: " << delta / count << std::endl;

    // build residual image
    unsigned char color[] = { 255 };
    unsigned char color2[] = { 255 , 0, 0 };
    cimg_library::CImg<unsigned char> image(angle_count, detector_size, 1, 1, 0);
    cimg_library::CImg<unsigned char> diff(angle_count, detector_size, 1, 3, 0);

    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            int b = fp[i * geom.detectorCount + j];
            color[0] = b;
            image.draw_point(i, j, color);
            int k = b > target(i, j, 0, 0) ? 0 : 1;
            color2[0] = 0;
            color2[1] = 0;
            color2[k] = abs(b - target(i, j, 0, 0));
            diff.draw_point(i, j, color2);
        }
    }

    // save normalized
    image.normalize(0, 255);
    image.save("out_sinogram.bmp");
    diff.normalize(0, 255);
    diff.save("out_residual_normed.bmp");

    // read ground truth sinogram from file
    std::ifstream input_f("sinogram.txt");
    std::unique_ptr<float[]> target_f(new float[angle_count * detector_size]);
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            float a;
            input_f >> a;
            target_f[i * detector_size + j] = a;
        }
    }
    input_f.close();

    // calculate sigma from ground truth
    delta = 0.0f;
    float total_delta = 0.0f;
    float out_delta = 1.0f, max_delta = 0.0f;
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            delta = target_f[i * detector_size + j] - myProj.fullProjection[i * detector_size + j];
            if (std::abs(delta) > max_delta) {
                max_delta = delta;
                out_delta = std::abs(delta) / myProj.fullProjection[i * detector_size + j];
            }
            float delta_sq = delta * delta;
            total_delta += delta_sq;
        }
    }
    float sigma = std::pow(total_delta/(angle_count * detector_size - 1), 0.5);
    std::cout << "sigma: " << sigma << " max_delta: " << max_delta << " delta/value: " << out_delta << std::endl;
}
