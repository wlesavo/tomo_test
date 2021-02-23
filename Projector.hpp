#ifndef _PROJECTOR_HPP_
#define _PROJECTOR_HPP_
#include "Line.h"
#include "MyImg.h"
#include "Geometry.h"
#include "Point.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

class Projector {
public:
	const MyImg inputImg;
	const Geometry geometry;
	const int imgSize_x, imgSize_y;
	std::unique_ptr<float[]> fullProjection = nullptr;
	Projector(std::unique_ptr<unsigned char[]> inputImg, int imgSize_x, int imgSize_y, const Geometry& geometry);
	float sumLine(const Line& line) const;
	float sumNeibs(int i_min, int i_max, int j, int axis = 0) const;
	std::unique_ptr<float[]> getSingleProjection(const float& angle) const;
	void buildFullProjection();
	std::unique_ptr<unsigned char[]> getFullProjectionImage();
};

#endif