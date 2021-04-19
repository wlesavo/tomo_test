#include "Projector.hpp"
#include <iostream>

// constructor

Projector::Projector(std::unique_ptr<float[]> i_inputImg, int i_imgSize_x, int i_imgSize_y, std::unique_ptr<Geometry> i_geometry, SumAlgo sumAlgorithm)
	: imgSize_x(i_imgSize_x), imgSize_y(i_imgSize_y), inputImg(MyImg(std::move(i_inputImg), i_imgSize_x, i_imgSize_y)), geometry(std::move(i_geometry)), sumAlgorithm(sumAlgorithm){};


// common methods

std::pair<Point, Point> Projector::getIntersectionPoints(const Line& line, int pixel_i, int pixel_j) const {
	// Method to get intersection points of line and image or pixel edges.
	// Default method will get intersection points of line and image edges.

	double coors[2][2] = {};
	int count = 0;
	int min_x, min_y, max_x, max_y;
	// get intersection points
	
	if (pixel_i == -999) {
		min_x = 0;
		min_y = 0;
		max_x = imgSize_x;
		max_y = imgSize_y;
	}
	else {
		min_x = pixel_i;
		min_y = pixel_j;
		max_x = pixel_i + 1;
		max_y = pixel_j + 1;
	}

	for (int i = 0; i < 2; ++i) {
		double x = (i == 1) ? min_x : max_x;
		double y = line.value(x);
		double y1 = (i == 1) ? min_y : max_y;
		double x1 = line.coor(y1);
		if (y > min_y && y < max_y) {
			if (count < 2) {
				if (x != coors[0][0] || y != coors[0][1]) {
					coors[count][0] = x;
					coors[count][1] = y;
					count += 1;
				}
			}
		}
		if (x1 >= min_x && x1 <= max_x) {
			if (count < 2) {
				if (x1 != coors[0][0]|| y1 != coors[0][1]) {
					coors[count][0] = x1;
					coors[count][1] = y1;
					count += 1;
				}
			}
		}
	}

	if (count < 2)
		return std::pair<Point, Point>(Point(), Point());

	Point point1, point2;
	int sign = (line.k > 0) - (line.k < 0);
	int i, j, i_max, j_max;

	// order points

	if (std::abs(line.k) >= 1.0) {
		if (coors[0][1] < coors[1][1]) {
			j = std::floor(coors[0][1]) - 1;
			j_max = std::ceil(coors[1][1]) + 1;
			if (sign > 0) {
				i = std::floor(coors[0][0]) - 1;
				i_max = std::ceil(coors[1][0]) + 1;
			}
			else {
				i = std::ceil(coors[0][0]) + 1;
				i_max = std::floor(coors[1][0]) - 1;
			}
		}
		else {
			j = std::floor(coors[1][1]) - 1;
			j_max = std::ceil(coors[0][1]) + 1;
			if (sign > 0) {
				i = std::floor(coors[1][0]) - 1;
				i_max = std::ceil(coors[0][0]) + 1;
			}
			else {
				i = std::ceil(coors[1][0]) + 1;
				i_max = std::floor(coors[0][0]) - 1;
			}
		}
	}
	else {
		if (coors[0][0] < coors[1][0]) {
			i = std::floor(coors[0][0]) - 1;
			i_max = std::ceil(coors[1][0]) + 1;
			if (sign > 0) {
				j = std::floor(coors[0][1]) - 1;
				j_max = std::ceil(coors[1][1]) + 1;
			}
			else {
				j = std::ceil(coors[0][1]) + 1;
				j_max = std::floor(coors[1][1]) - 1;
			}
		}
		else {
			i = std::floor(coors[1][0]) - 1;
			i_max = std::ceil(coors[0][0]) + 1;
			if (sign > 0) {
				j = std::floor(coors[1][1]) - 1;
				j_max = std::ceil(coors[0][1]) + 1;
			}
			else {
				j = std::ceil(coors[1][1]) + 1;
				j_max = std::floor(coors[0][1]) - 1;
			}
		}
	}

	// todo, fix: potential error for area methods
	point1 = Point(i, j);
	point2 = Point(i_max, j_max);
	return std::pair<Point, Point>(point1, point2);
	
}

//std::pair<Point, Point> Projector::getIntersectionPoints(const Line& line, int pixel_i, int pixel_j) const {
//	// Method to get intersection points of line and image or pixel edges.
//	// Default method will get intersection points of line and image edges.
//
//	int count = 0;
//	int min_x, min_y, max_x, max_y;
//	// get intersection points
//
//	if (pixel_i == -999) {
//		min_x = 0;
//		min_y = 0;
//		max_x = imgSize_x;
//		max_y = imgSize_y;
//	}
//	else {
//		min_x = pixel_i;
//		min_y = pixel_j;
//		max_x = pixel_i + 1;
//		max_y = pixel_j + 1;
//	}
//
//	double x1 = line.coor(min_y);
//	double x2 = line.coor(max_y);
//	double y1 = line.value(min_x);
//	double y2 = line.value(max_x);
//	int count = 0;
//	bool bottom = (x1 >= pixel_i && x1 <= (pixel_i + 1));
//	bool top = (x2 >= pixel_i && x2 <= (pixel_i + 1));
//	bool left = (y1 >= pixel_j && y1 <= (pixel_j + 1));
//	bool right = (y2 >= pixel_j && y2 <= (pixel_j + 1));
//	count += (int)top + (int)bottom + (int)right + (int)left;
//
//	if (count > 2) {
//		if (top && bottom) {
//			left = false;
//			right = false;
//		}
//		else if (left && right) {
//			top = false;
//			bottom = false;
//		}
//		count = 2;
//	}
//	
//	if (count < 2)
//		return std::pair<Point, Point>(Point(), Point());
//
//	Point point1, point2;
//	int sign = (line.k > 0) - (line.k < 0);
//	int i, j, i_max, j_max;
//
//	// order points
//
//	if (std::abs(line.k) >= 1.0) {
//		if (coors[0][1] < coors[1][1]) {
//			j = std::floor(coors[0][1]) - 1;
//			j_max = std::ceil(coors[1][1]) + 1;
//			if (sign > 0) {
//				i = std::floor(coors[0][0]) - 1;
//				i_max = std::ceil(coors[1][0]) + 1;
//			}
//			else {
//				i = std::ceil(coors[0][0]) + 1;
//				i_max = std::floor(coors[1][0]) - 1;
//			}
//		}
//		else {
//			j = std::floor(coors[1][1]) - 1;
//			j_max = std::ceil(coors[0][1]) + 1;
//			if (sign > 0) {
//				i = std::floor(coors[1][0]) - 1;
//				i_max = std::ceil(coors[0][0]) + 1;
//			}
//			else {
//				i = std::ceil(coors[1][0]) + 1;
//				i_max = std::floor(coors[0][0]) - 1;
//			}
//		}
//	}
//	else {
//		if (coors[0][0] < coors[1][0]) {
//			i = std::floor(coors[0][0]) - 1;
//			i_max = std::ceil(coors[1][0]) + 1;
//			if (sign > 0) {
//				j = std::floor(coors[0][1]) - 1;
//				j_max = std::ceil(coors[1][1]) + 1;
//			}
//			else {
//				j = std::ceil(coors[0][1]) + 1;
//				j_max = std::floor(coors[1][1]) - 1;
//			}
//		}
//		else {
//			i = std::floor(coors[1][0]) - 1;
//			i_max = std::ceil(coors[0][0]) + 1;
//			if (sign > 0) {
//				j = std::floor(coors[1][1]) - 1;
//				j_max = std::ceil(coors[0][1]) + 1;
//			}
//			else {
//				j = std::ceil(coors[1][1]) + 1;
//				j_max = std::floor(coors[0][1]) - 1;
//			}
//		}
//	}
//
//	// todo, fix: potential error for area methods
//	point1 = Point(i, j);
//	point2 = Point(i_max, j_max);
//	return std::pair<Point, Point>(point1, point2);
//}
std::pair<Line, Line> Projector::sortLines(const Line& line1, const Line& line2) const {
	// return pair (upper, lower) lines. for vertical lines return (right, left) lines

	if (line1.b > line2.b) {
		return std::pair<Line, Line>(line1, line2);
	}
	else if (line1.b < line2.b){
		return std::pair<Line, Line>(line2, line1);
	}
	else if (line1.b == line2.b) {
		// todo : correct upper lower check
		return std::pair<Line, Line>(line1, line2);
	}
}

float Projector::sumNeibs(double i_min, double i_max, double j, int axis) const {
	// summation along a horizontal (axis = 0) or vertical (axis = 1) direction
	
	float summ = 0;
	if (axis == 0) {
		if (i_min < 0)
			i_min = -1;
		if (i_max >= inputImg.size_x)
			i_max = inputImg.size_x;
		for (int i = i_min; i < i_max; i++) {
			summ += inputImg.get(i, j);
		}
	}
	else if (axis == 1) {
		if (i_min < 0)
			i_min = - 1;
		if (i_max >= inputImg.size_y)
			i_max = inputImg.size_y;
		for (int i = i_min; i < i_max; i++) {
			summ += inputImg.get(j, i);
		}
	}
	return summ;
}

double Projector::singlePixelArea(int pixel_i, int pixel_j, const Line& line) const {
	// return the part of the pixel above the line
	// (1 - sum) would denote the part below the line due to whole area being 1.0
	double sum = 0.0;
	double x1 = line.coor(pixel_j);
	double x2 = line.coor((pixel_j + 1));
	double y1 = line.value(pixel_i);
	double y2 = line.value((pixel_i + 1));
	int count = 0;
	bool bottom = (x1 >= pixel_i && x1 <= (pixel_i + 1));
	bool top = (x2 >= pixel_i && x2 <= (pixel_i + 1));
	bool left = (y1 >= pixel_j && y1 <= (pixel_j + 1));
	bool right = (y2 >= pixel_j && y2 <= (pixel_j + 1));
	count += (int)top + (int)bottom + (int)right + (int)left;
	if (count > 2) {
		if (top && bottom) {
			left = false;
			right = false;
		}
		else if (left && right) {
			top = false;
			bottom = false;
		}
		count = 2;
	}

	x1 -= std::floor(x1);
	x2 = 1 - (std::ceil(x2) - x2);
	y1 -= std::floor(y1);
	y2 = 1 - (std::ceil(y2) - y2);

	// todo: check for upper and lower pixel positions // fixed
	if (count < 2) {
		if ((pixel_j - line.k * pixel_i) > line.b) {
			sum = 1.0;
		}
		else {
			sum = 0.0;
		}
		return sum;
	}
	if (line.k > 0) {
		if (bottom && top) {
			sum = (x2 + x1) / 2;
		}
		else if (left && right) {
			sum = 1 - (y2 + y1) / 2;
		}
		else if (bottom) {
			assert(right);
			sum = 1 - (1 - x1) * (y2) / 2;
		}
		else if (top) {
			assert(left);
			sum = (1 - y1) * (x2) / 2;
		}
	}
	else if (line.k < 0) {
		if (bottom && top) {
			sum = 1 - (x2 + x1) / 2;
		}
		else if (left && right) {
			sum = 1 - (y2 + y1) / 2;
		}
		else if (bottom) {
			assert(left);
			sum = 1 - (x1 * y1) / 2;
		}
		else if (top) {
			assert(right);
			sum = (1 - x2) * (1 - y2) / 2;
		}
	}
	
	assert(0.0 <= sum && sum <= 1.0);
	
	return sum;
}

float Projector::manyPixelArea(int i_min, int i_max, int j, bool upper, const Line& line) const {
	float sum = 0.0f;
	if (j < 0 || j >= inputImg.size_y) {
		return sum;
	}
	if (i_min < 0)
		i_min = -1;
	if (i_max >= inputImg.size_x)
		i_max = inputImg.size_x;
	for (int i = i_min; i < i_max; i++) {
		if (upper) {
			double a = singlePixelArea(i, j, line);
			float b = inputImg.get(i, j);
			sum += a * b;
			//sum += singlePixelArea(i, j, line) * inputImg.get(i, j);
		}
		else if (!upper) {
			sum += (1.0 - singlePixelArea(i, j, line)) * inputImg.get(i, j);
		}
	}
	return sum;	
}


// sum algos

float Projector::sumLine(const Line& line) const {

	if (line.hor)
		return sumNeibs(-1, imgSize_x, line.b, 0);
	else if (line.vert)
		return sumNeibs(-1, imgSize_y, line.b, 1);

	std::pair<Point, Point> intersect = getIntersectionPoints(line);

	if (intersect.first.x == -999)
		return 0.0f;

	float sum = 0.0f;
	// choose axis to sum along based on line parameters
	// axis = 1; dy = 0, dx = 1 || axis = 0; dy = 1, dx = 0
	int axis = (std::abs(line.k) >= 1.0);
	int sign = (line.k > 0) - (line.k < 0);
	int i, j, i_max, j_max;
	i = intersect.first.x;
	j = intersect.first.y;
	i_max = intersect.second.x;
	j_max = intersect.second.y;

	// main summation cycle
	// summation is performed by sequencually finding next intersections of line and 
	// vertical (axis = 0) or horizontal (axis = 1) pixel border, 
	// summing up to this point and adding proportional values for the intersection area
	double new_x, new_y, frac;
	int new_i, new_j;

	while (true) {
		if (axis == 1) { // dy == 0
			new_i = i + sign;
			if (sign > 0) {
				new_y = line.value(new_i);
			}
			else {
				new_y = line.value(i);
			}
			new_j = std::floor(new_y);
			sum += sumNeibs(j, new_j, i, axis);
			frac = std::abs(new_y - new_j);
			float b = inputImg.get(i, new_j) * frac + inputImg.get(new_i, new_j) * (1 - frac);
			if (b > 0) {
				float a = 0;
			}
			sum += b;
			i = new_i;
			j = new_j + 1;
			assert(sum >= 0);
			if (j > j_max)
				return sum / std::abs(std::sin(line.angle));
		}
		else if (axis == 0) { // dx == 0
			new_j = j + sign;
			if (sign > 0)
				new_x = line.coor(new_j);
			else
				new_x = line.coor(j);
			new_i = std::floor(new_x);

			float a = sumNeibs(i, new_i, j, axis);
			sum += a;
			frac = std::abs(new_x - new_i);
			float b = inputImg.get(new_i, j) * frac + inputImg.get(new_i, new_j) * (1 - frac);
			sum += b;
			j = new_j;
			i = new_i + 1;
			assert(sum >= 0);
			if (i > i_max)
				return sum / std::abs(std::cos(line.angle));
		}
	}
}

float Projector::sumLinear(const Line& line) const {

	if (line.hor)
		return sumNeibs(-1, imgSize_x, line.b, 0);
	else if (line.vert)
		return sumNeibs(-1, imgSize_y, line.b, 1);

	std::pair<Point, Point> intersect = getIntersectionPoints(line);

	if (intersect.first.x == -999)
		return 0.0f;

	float sum = 0.0f;
	// choose axis to sum along based on line parameters
	// axis = 1; dy = 0, dx = 1 || axis = 0; dy = 1, dx = 0
	int axis = (std::abs(line.k) >= 1.0);
	int sign = (line.k > 0) - (line.k < 0);
	int i, j, i_max, j_max;
	i = intersect.first.x;
	j = intersect.first.y;
	i_max = intersect.second.x;
	j_max = intersect.second.y;

	// main summation cycle
	// summation is performed by sequencually finding next intersections of line and 
	// vertical (axis = 0) or horizontal (axis = 1) pixel border, 
	// summing up to this point and adding proportional values for the intersection area
	double new_x, new_y, frac;
	int new_i, new_j;
	double left, right;
	while (true) {
		if (axis == 1) { // dy == 0
			new_x = line.coor(j + 0.5);
			left = std::floor(new_x - 0.5);
			right = left + 1;
			frac = std::abs(new_x - 0.5 - left);
			sum += inputImg.get(left, j) * (1 - frac) + inputImg.get(right, j) * frac;
			j = j + 1;
			assert(sum >= 0);
			if (j > j_max + 1)
				return sum / std::abs(std::sin(line.angle));
		}
		else if (axis == 0) { // dx == 0
			new_y = line.value(i + 0.5);
			left = std::floor(new_y - 0.5);
			right = left + 1;
			frac = std::abs(new_y - 0.5 - left);
			sum += inputImg.get(i, left) * (1 - frac) + inputImg.get(i, right) * frac;
			i = i + 1;
			assert(sum >= 0);
			if (i > i_max + 1)
				return sum / std::abs(std::cos(line.angle));
		}
	}
}

float Projector::sumArea(const Line& line_1, const Line& line_2) const {
	float sum = 0.0f;
	float sum1, sum2, sum3;
	double x_left_1, x_right_1;
	double x_left_2, x_right_2;
	// todo: fan geometry hor and vert handling
	if (line_1.hor && line_2.hor)
	{
		sum1 = sumNeibs(-1, imgSize_x, std::floor(line_1.b), 0) *      (line_1.b - std::floor(line_1.b));
		sum2 = sumNeibs(-1, imgSize_x, std::floor(line_2.b), 0) * (1 - (line_2.b - std::floor(line_2.b)));
		return sum1 + sum2;
	}
	else if (line_1.vert && line_2.vert)
	{

		sum1 = sumNeibs(-1, imgSize_y, std::floor(line_1.b), 1) *      (line_1.b - std::floor(line_1.b));
		sum2 = sumNeibs(-1, imgSize_y, std::floor(line_2.b), 1) * (1 - (line_2.b - std::floor(line_2.b)));
		return sum1 + sum2;
	}

	std::pair<Point, Point> p_1, p_2;
	
	// todo: check correctness of intersection ordering
	p_1 = getIntersectionPoints(line_1);
	p_2 = getIntersectionPoints(line_2);
	bool sign_1 = line_1.k > 0, sign_2 = line_2.k > 0;

	int j, j_max;
	j = std::min({ p_1.first.y, p_1.second.y, p_2.first.y, p_2.second.y });
	j_max = std::max({ p_1.first.y, p_1.second.y, p_2.first.y, p_2.second.y });

	if (p_2.first.x == -999) {
		j = std::min(p_1.first.y, p_1.second.y);
	}
	assert(j_max >= j);
	
	//while (true) {
	//	if (sign > 0) {
	//		x_left_u = line_up.coor(j);
	//		x_right_u = line_up.coor(j + 1);
	//		x_left_l = line_low.coor(j);
	//		x_right_l = line_low.coor(j + 1);
	//		sum1 = sumNeibs(std::floor(x_left_u), std::ceil(x_right_l), j, 0);
	//	}
	//	else {
	//		x_left_u = line_up.coor(j + 1);
	//		x_right_u = line_up.coor(j);
	//		x_left_l = line_low.coor(j + 1);
	//		x_right_l = line_low.coor(j);
	//		sum1 = sumNeibs(std::floor(x_left_l), std::ceil(x_right_u), j, 0);
	//	}
	//	//if (sum1 > 0) {
	//	//	float a = 0;
	//	//}
	//	sum2 = manyPixelArea(std::floor(x_left_u), std::ceil(x_right_u), j, true, line_up);
	//	sum3 = manyPixelArea(std::floor(x_left_l), std::ceil(x_right_l), j, false, line_low);
	//	sum += std::abs(sum1 - sum2 - sum3);
	//	j += 1;
	//	if (j > j_max)
	//		return sum;
	//}
	double left, right;
	while (true) {

		x_left_1 = line_1.coor(j);
		x_right_1 = line_1.coor(j + 1);
		x_left_2 = line_2.coor(j);
		x_right_2 = line_2.coor(j + 1);

		left = std::min({ x_left_1, x_right_1, x_left_2, x_right_2 });
		right = std::max({ x_left_1, x_right_1, x_left_2, x_right_2 });

		sum2 = manyPixelArea(std::floor(left), std::ceil(right), j, sign_1, line_1);
		sum3 = manyPixelArea(std::floor(left), std::ceil(right), j, sign_2, line_2);
		sum += std::abs(sum2 - sum3);
		j += 1;
		if (j > j_max)
			return sum;
	}


}


// main methods

std::unique_ptr<float[]> Projector::getSingleProjection(int angle_i) const {
	// summation along particular angle of projection for every pixel of detector

	std::unique_ptr<float[]> outLineImg(new float[geometry->detectorCount]);
	Line line, line1, line2;
	std::pair<Line, Line> sorted_lines;
	for (int i = 0; i < geometry->detectorCount; ++i) {
		float out = 0;
		
		switch (sumAlgorithm)
		{
		case SumAlgo::LINE:
			line = geometry->v_GetNextLineCenter(angle_i, i);
			out = sumLine(line);
			break;
		case SumAlgo::LINEAR:
			line = geometry->v_GetNextLineCenter(angle_i, i);
			out = sumLinear(line);
			break;
		case SumAlgo::BINARY:
			break;
		case SumAlgo::AREA:
			line1 = geometry->v_GetNextLineEdge(angle_i, i);
			line2 = geometry->v_GetNextLineEdge(angle_i, i+1);
			sorted_lines = sortLines(line1, line2);
			out = sumArea(sorted_lines.first, sorted_lines.second);
			break;
		default:
			break;
		}

		outLineImg[i] = out;
	}

	return outLineImg;
}

void Projector::buildFullProjection(){
	// build projection for every angle
	// x axis is angle to detector as an index in input 'angles' vector
	// y axis is a detector pixel

	std::unique_ptr<float[]> outfloatImg(new float[geometry->detectorCount * geometry->angles.size()]);
	for (int i = 0; i < geometry->angles.size(); ++i) {
		std::unique_ptr<float[]> lineImg = getSingleProjection(i);
		for (int j = 0; j < geometry->detectorCount; ++j) {
			outfloatImg[i * geometry->detectorCount + j ] = lineImg[j];
		}
	}
	fullProjection = std::move(outfloatImg);
}

std::unique_ptr<unsigned char[]> Projector::getFullProjectionImage() {
	// normalize and return full projection in an image format
	// x axis is angle to detector as an index in input 'angles' vector
	// y axis is a detector pixel

	if (fullProjection == nullptr)
		buildFullProjection();

	std::unique_ptr<unsigned char[]> outCharImg(new unsigned char[geometry->detectorCount * geometry->angles.size()]);
	float max_val = *std::max_element(&fullProjection[0], &fullProjection[geometry->detectorCount * geometry->angles.size() - 1]);
	// check for too small max value; if (max < 1.0) then max = 1.0
	float reverse_max = (max_val > 1.0f) ? (1 / max_val) : 1.0f;
	for (int i = 0; i < geometry->angles.size(); ++i) {
		for (int j = 0; j < geometry->detectorCount; ++j) {
			outCharImg[i * geometry->detectorCount + j] = (unsigned char)(fullProjection[i * geometry->detectorCount + j] * 255.0f * reverse_max);
		}
	}
	return outCharImg;
};


// debug methods

float Projector::getLineProjectionTest(const int& i_angle, const int& detector) {
	// test function for checking particular angle and detector

	float out = 0.0f;
	Line line, line1, line2;
	std::pair<Line, Line> sorted_lines;

	switch (sumAlgorithm)
	{
	case SumAlgo::LINE:
		line = geometry->v_GetNextLineCenter(i_angle, detector);
		out = sumLine(line);
		break;
	case SumAlgo::LINEAR:
		line = geometry->v_GetNextLineCenter(i_angle, detector);
		out = sumLinear(line);
		break;
	case SumAlgo::BINARY:
		break;
	case SumAlgo::AREA:
		line1 = geometry->v_GetNextLineEdge(i_angle, detector);
		line2 = geometry->v_GetNextLineEdge(i_angle, detector + 1);
		sorted_lines = sortLines(line1, line2);
		out = sumArea(sorted_lines.first, sorted_lines.second);
		break;
	default:
		break;
	}
	double angle = geometry->angles[i_angle];
	out = sumLineTest(geometry->v_GetNextLineCenter(angle, detector));

	//out = sumLine(geometry->v_GetNextLine(angle, detector));
	//out = sumLineLinear(geometry->v_GetNextLine(angle, detector));

	return out;
}

float Projector::sumLineTest(const Line& line) {

	testFlag = false;

	if (line.hor)
		return sumNeibs(-1, imgSize_x, line.b, 0);
	else if (line.vert)
		return sumNeibs(-1, imgSize_y, line.b, 1);

	std::pair<Point, Point> intersect = getIntersectionPoints(line);

	if (intersect.first.x == -999)
		return 0.0f;

	float sum = 0.0f;
	// choose axis to sum along based on line parameters
	// axis = 1; dy = 0, dx = 1 || axis = 0; dy = 1, dx = 0
	int axis = (std::abs(line.k) >= 1.0);
	int sign = (line.k > 0) - (line.k < 0);
	int i, j, i_max, j_max;
	i = intersect.first.x;
	j = intersect.first.y;
	i_max = intersect.second.x;
	j_max = intersect.second.y;

	// main summation cycle
	// summation is performed by sequencually finding next intersections of line and 
	// vertical (axis = 0) or horizontal (axis = 1) pixel border, 
	// summing up to this point and adding proportional values for the intersection area
	float new_x, new_y, frac;
	int new_i, new_j;

	while (true) {
		if (axis == 1) { // dy == 0
			new_i = i + sign;
			if (sign > 0)
				new_y = line.value(new_i);
			else
				new_y = line.value(i);
			new_j = std::floor(new_y);
			sum += sumNeibs(j, new_j, i, axis);
			frac = std::abs(new_y - new_j);
			float b = inputImg.get(i, new_j) * frac + inputImg.get(new_i, new_j) * (1 - frac);
			sum += b;
			if (b > 0)
				testFlag = true;
			i = new_i;
			j = new_j + 1;
			assert(sum >= 0);
			if (j > j_max)
				return sum / std::abs(std::sin(line.angle));
		}
		else if (axis == 0) { // dx == 0
			new_j = j + sign;
			if (sign > 0)
				new_x = line.coor(new_j);
			else
				new_x = line.coor(j);
			new_i = std::floor(new_x);

			float a = sumNeibs(i, new_i, j, axis);
			sum += a;
			frac = std::abs(new_x - new_i);
			float b = inputImg.get(new_i, j) * frac + inputImg.get(new_i, new_j) * (1 - frac);
			if (b > 0)
				testFlag = true;
			sum += b;
			j = new_j;
			i = new_i + 1;
			assert(sum >= 0);
			if (i > i_max)
				return sum / std::abs(std::cos(line.angle));
		}
	}
}
