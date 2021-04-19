#ifndef _PROJECTOR_HPP_
#define _PROJECTOR_HPP_
#include "Line.h"
#include "MyImg.h"
#include "Geometry.h"
#include "Point.h"
#include "SumAlgo.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

class Projector {
public:
	// main objects
	const MyImg inputImg;
	std::unique_ptr<Geometry> geometry;
	SumAlgo sumAlgorithm;
	std::unique_ptr<float[]> fullProjection;
	const int imgSize_x, imgSize_y;

	// debug vars
	bool testFlag = false;

	// constructor
	Projector(std::unique_ptr<float[]> inputImg, int imgSize_x, int imgSize_y, std::unique_ptr<Geometry> geometry, SumAlgo sumAlgorithm);
	
	// main methods
	void buildFullProjection();
	std::unique_ptr<unsigned char[]> getFullProjectionImage();
	std::unique_ptr<float[]> getSingleProjection(int angle_i) const;

	// common methods
	std::pair<Point, Point> getIntersectionPoints(const Line& line, int pixel_i = -999, int pixel_j = -999) const;
	float sumNeibs(double i_min, double i_max, double j, int axis = 0) const;
	std::pair<Line, Line> sortLines(const Line& line1, const Line& line2) const;
	float manyPixelArea(int i_min, int i_max, int j, bool upper, const Line& line) const;
	double singlePixelArea(int i, int j, const Line& line) const;

	// sum algos
	float sumLine(const Line& line) const;
	float sumLinear(const Line& line) const;
	float sumArea(const Line& line1, const Line& line2) const;
	
	// debug methods
	float sumLineTest(const Line& line);
	float getLineProjectionTest(const int& angle, const int& detector);


};

#endif