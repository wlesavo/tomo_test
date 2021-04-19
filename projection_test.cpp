#include "CImg.h"
#include "Projector.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include "Geometries.h"

std::vector<double> getAngles(double start, double end, double delta) {
    std::vector<double> angles;
    double current = start;
    while (current < end) {
        angles.push_back(current * cimg_library_suffixed::cimg::PI / 180.0);
        current += delta;
    }
    return angles;
}

void customLineProjections() {
    std::vector<double> angles = getAngles(0.0, 180.0, 5.0);
    int angle_count = angles.size();
    int detector_size = 256;
    int size_x = 256, size_y = 256;
    for (int k = 1; k <= 14; ++k) {
        std::string s = std::to_string(k);
        std::ifstream input_f("phantom\\"+s+".txt");
        std::unique_ptr<float[]> img(new float[size_x * size_y]);

        unsigned char color[] = { 255 };
        cimg_library::CImg<float> cimage(size_x, size_y, 1, 1, 0);

        for (int i = 0; i < size_x; ++i) {
            for (int j = 0; j < size_y; ++j) {
                float a;
                input_f >> a;
                img[i * size_x + j] = a;
                //color[0] = a*255;
                //cimage.draw_point(i, j, color);
            }
        }

        input_f.close();
        Projector myProj(std::move(img), size_x, size_y, std::make_unique<GeometryParallel>(GeometryParallel(angles, detector_size, 1.0, size_x * 0.5, size_y * 0.5))
            , SumAlgo::LINE);
        std::ofstream output_f("phantom\\out_sino\\" + s + "_small.txt");
        //std::ofstream output_f("phantom\\out_sino\\" + s + "_large.txt");
        output_f.setf(std::ios::scientific);
        output_f.precision(std::numeric_limits<double>::digits10 + 1);
        for (int i = 0; i < angle_count; ++i) {
            for (int j = 0; j < detector_size; ++j) {
                Line line = myProj.geometry->v_GetNextLineCenter(angles[i], j);
                float sum = myProj.sumLine(line);
                float x1 = 0.0f, x2 = 256.0f;
                float y1 = line.value(x1), y2 = line.value(x2);
                if (std::abs(line.k) > 1) {
                    y1 = 0.0f;
                    y2 = 256.0f;
                    x1 = line.coor(y1);
                    x2 = line.coor(y2);
                }
                x1 -= size_x * 0.5f;
                x2 -= size_x * 0.5f;
                y1 -= size_y * 0.5f;
                y2 -= size_y * 0.5f;
                if (line.vert) {
                    x1 = line.b - size_x * 0.5f;
                    x2 = line.b - size_x * 0.5f;
                    y1 =        - size_y * 0.5f;
                    y2 =          size_y * 0.5f;
                }
                if (line.hor) {
                    x1 =        - size_x * 0.5f;
                    x2 =          size_x * 0.5f;
                    y1 = line.b - size_y * 0.5f;
                    y2 = line.b - size_y * 0.5f;
                }
                output_f << x1 << " " << y1 << " " << x2 << " " << y2 << " " << sum << std::endl;
            }
        }
        output_f.close();
        //cimage.display();
        std::cerr << k << " finished" << std::endl;

    }
    return;
}

int main()
{   
    //customLineProjections();
    //std::cout << std::numeric_limits<double>::is_iec559 << " " << std::numeric_limits<double>::digits <<  std::endl;
    double delta_limit = 110000.1;
    double distance_obj = 500.0, distance_source = 7000.0;

    std::string model           = "";      // "", "_model";
    std::string algo_string     = "_line"; // "_area", "_line";
    std::string geometry_string = "_par";  // "_fan", "_par";

    std::string p = "assets\\phantom"  + model + geometry_string + algo_string + ".bmp";
    std::string t = "assets\\sinogram" + model + geometry_string + algo_string + ".bmp";
    std::string target_txt_path = "assets\\sinogram" + model + geometry_string + algo_string + ".txt";

    const char* const phantom_path    = p.c_str();
    const char* const target_path     = t.c_str();

    cimg_library::CImg<unsigned char> phantom(phantom_path), target(target_path);
    std::ifstream input_f(target_txt_path);

    int size_x = phantom.width(), size_y = phantom.height();
    int detector_size = target.height(), angle_count = target.width();
    std::unique_ptr<float[]> img(new float[size_x * size_y]{});

    // copy first channel of image to array assuming grayscale bmp
    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            img[i + size_x * j] = (float)phantom(j, i, 0, 0);
        }
    }
    //std::ifstream input_angles("angles.txt");
    //std::vector<double> angles;
    //input_angles >> angle_count;
    //for (int i = 0; i < angle_count; ++i) {
    //    double a;
    //    input_angles >> a;
    //    angles.push_back(a);
    //}

    std::vector<double> angles;
    for (int i = 0; i < angle_count; ++i) {
        angles.push_back((double)(i) * M_PI / 180.0);
    }
    
    auto sum_algo = SumAlgo::LINE;
    if (algo_string == "_area") {
        sum_algo = SumAlgo::AREA;
    }
    else if (algo_string == "_line") {
        sum_algo = SumAlgo::LINE;
    }
    std::unique_ptr<Geometry> input_geometry = nullptr;
    if (geometry_string == "_fan") {
        input_geometry = std::move(std::make_unique<GeometryFanBeam>(GeometryFanBeam(angles, detector_size, 1.0, distance_source, distance_obj, size_x * 0.5, size_y * 0.5)));
    }
    else {
        input_geometry = std::move(std::make_unique<GeometryParallel>(GeometryParallel(angles, detector_size, 1.0, size_x * 0.5, size_y * 0.5)));
    }

    Projector myProj(std::move(img), size_x, size_y, std::move(input_geometry), sum_algo);

    std::unique_ptr<unsigned char[]> fp = myProj.getFullProjectionImage();

    // compare the differences in image form,
    // count non-zero differences in pixels for statistics
    float delta = 0.0f, summ = 0.0f;
    int count = 0;
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            int a = (int)target(i, j, 0, 0);
            int b = (int)fp[i * detector_size + (j)];
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
            int b = fp[i * detector_size + j];
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
    int pos_delta = 0, neg_delta = 0;
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            delta = target_f[i * detector_size + j] - myProj.fullProjection[i * detector_size + j];
            if (std::abs(delta) > max_delta) {
                max_delta = std::abs(delta);
                out_delta = std::abs(delta) / myProj.fullProjection[i * detector_size + j];
            }
            float delta_sq = delta * delta;
            total_delta += delta_sq;
            if (delta == 0)
                continue;
            if (delta > 0)
                pos_delta += 1;
            else
                neg_delta += 1;
            if (abs(delta) > delta_limit) {
                std::cerr << "delta: " << delta << " target: " << target_f[i * detector_size + j] << " res: "
                    << myProj.fullProjection[i * detector_size + j] << " det: " << j << " ang: " << i << std::endl;

                myProj.getLineProjectionTest(i, j);
                /*myProj.getLineProjectionTest(i, j);
                Line line = myProj.geometry->v_GetNextLineCenter(myProj.geometry->angles[i], j);

                std::cerr << "full path ";
                if (myProj.testFlag) {
                    std::cerr << " - ";
                }
                else {
                    std::cerr << " + ";
                }
                
                std::cerr << "delta: " << delta << " target: " << target_f[i * detector_size + j] << " res: " 
                    << myProj.fullProjection[i * detector_size + j] << " det: " << j << " ang: " << i << " line_k " << line.k << std::endl;
                */
            }
        }
    }
    float sigma = std::pow(total_delta/(angle_count * detector_size - 1), 0.5);
    std::cout << "sigma: " << sigma << " max_delta: " << max_delta << " delta/value: " << out_delta << std::endl;
    std::cout << "pos: " << pos_delta << " neg: " << neg_delta << std::endl;


    // calculate residual picture from values
    unsigned char color3[] = { 255 , 0, 0 };
    cimg_library::CImg<unsigned char> diff_float(angle_count, detector_size, 1, 3, 0);
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            float value = (target_f[i * detector_size + j] - myProj.fullProjection[i * detector_size + j]) / max_delta; ;
            int k = value > 0 ? 0 : 1;
            color3[0] = 0;
            color3[1] = 0;
            color3[k] = abs(value*255);
            diff_float.draw_point(i, j, color3);
        }
    }
    diff_float.normalize(0, 255);
    diff_float.save("out_residual_float.bmp");
}
