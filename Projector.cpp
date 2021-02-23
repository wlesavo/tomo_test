#include "Projector.hpp"

Projector::Projector(std::unique_ptr<unsigned char[]> i_inputImg, int i_imgSize_x, int i_imgSize_y, const Geometry& i_geometry)
	: imgSize_x(i_imgSize_x), imgSize_y(i_imgSize_y), inputImg(MyImg(std::move(i_inputImg), i_imgSize_x, i_imgSize_y)), geometry(i_geometry) {};

float Projector::sumLine(const Line& line) const {

	if (line.hor)
		return sumNeibs(-inputImg.center_x-1, inputImg.center_x+1, line.dist, 0);
	else if (line.vert)
		return sumNeibs(-inputImg.center_y-1, inputImg.center_y+1, line.dist, 1);

	float sum = 0.0f;

	float coors[2][2] = {};
	int count = 0;

	// get intersection points
	for (int i = 0; i < 2; ++i) {
		float x = (inputImg.center_x) * (2.0f * i - 1);
		float y = line.value(x);
		float y1 = (inputImg.center_y) * (2.0f * i - 1);
		float x1 = line.coor(y1);
		if (y >= -inputImg.center_y && y <= inputImg.center_y) {
			if (count < 2) {
				if (coors[0][0] != x || coors[0][1] != y) {
					coors[count][0] = x;
					coors[count][1] = y;
					count += 1;
				}
			}
		}
		if (x1 >= -inputImg.center_x && x1 <= inputImg.center_x) {
			if (count < 2) {
				if (coors[0][0] != x1 || coors[0][1] != y1) {
					coors[count][0] = x1;
					coors[count][1] = y1;
					count += 1;
				}
			}
		}
	}

	if (count < 2)
		return 0.0f;

	int axis;
	int sign = (line.k > 0) - (line.k < 0);
	int i, j, i_max, j_max;

	// order points,
	// choose axis to sum along based on line parameters
	if (std::abs(line.k)>=1.0f) {
		axis = 1; //dy = 0, dx = 1
		if (coors[0][1] < coors[1][1]) {
			j = std::floor(coors[0][1]);
			j_max = std::ceil(coors[1][1]);
			if (sign > 0) {
				i = std::floor(coors[0][0]);
				i_max = std::ceil(coors[1][0]);
			}
			else {
				i = std::ceil(coors[0][0]);
				i_max = std::floor(coors[1][0]);
			}
		}
		else {
			j = std::floor(coors[1][1]);
			j_max = std::ceil(coors[0][1]);
			if (sign > 0) {
				i = std::floor(coors[1][0]);
				i_max = std::ceil(coors[0][0]);
			}
			else {
				i = std::ceil(coors[1][0]);
				i_max = std::floor(coors[0][0]);
			}
		}
	}
	else {
		axis = 0; //dy = 1, dx = 0
		if (coors[0][0] < coors[1][0]) {
			i = std::floor(coors[0][0]);
			i_max = std::ceil(coors[1][0]);
			if (sign > 0) {
				j = std::floor(coors[0][1]);
				j_max = std::ceil(coors[1][1]);
			}
			else {
				j = std::ceil(coors[0][1]);
				j_max = std::floor(coors[1][1]);
			}
		}
		else {
			i = std::floor(coors[1][0]);
			i_max = std::ceil(coors[0][0]);
			if (sign > 0) {
				j = std::floor(coors[1][1]);
				j_max = std::ceil(coors[0][1]);
			}
			else {
				j = std::ceil(coors[1][1]);
				j_max = std::floor(coors[0][1]);
			}
		}
	}
	
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
			sum += inputImg.get_c(i, new_j) * frac + inputImg.get_c(new_i, new_j) * (1 - frac);
			i = new_i;
			j = new_j+1;
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
			sum += sumNeibs(i, new_i, j, axis);
			frac = std::abs(new_x - new_i);
			sum += inputImg.get_c(new_i, j) * frac + inputImg.get_c(new_i, new_j) * (1 - frac);
			j = new_j;
			i = new_i + 1;
			assert(sum >= 0);
			if (i > i_max)
				return sum / std::abs(std::cos(line.angle));
		}
	}
}

float Projector::sumNeibs(int i_min, int i_max, int j, int axis) const {
	// summation along a horizontal (axis = 0) or vertical (axis = 1) direction
	
	float summ = 0;
	if (axis == 0) {
		i_min += inputImg.center_x;
		i_max += inputImg.center_x;
		j += inputImg.center_y;
		if (i_min < 0)
			i_min = -1;
		else if (i_max >= inputImg.size_x)
			i_max = inputImg.size_x;
		for (int i = i_min; i < i_max; i++) {
			summ += inputImg.get(i, j);
		}
	}
	else if (axis == 1) {
		i_min += inputImg.center_y;
		i_max += inputImg.center_y;
		j += inputImg.center_x;
		if (i_min < 0)
			i_min = - 1;
		else if (i_max >= inputImg.size_y)
			i_max = inputImg.size_y;
		for (int i = i_min; i < i_max; i++) {
			summ += inputImg.get(j, i);
		}
	}
	return summ;
}

std::unique_ptr<float[]> Projector::getSingleProjection(const float& angle) const {
	// summation along particular angle of projection for every pixel of detector

	std::unique_ptr<float[]> outLineImg(new float[geometry.detectorCount]);
	float x_source = -geometry.dSourceToObj * std::cos(angle);
	float y_source = -geometry.dSourceToObj * std::sin(angle);
	for (int i = 0; i < geometry.detectorCount; ++i) {
		float x_det = geometry.distanceToDetectorPixel[i] * std::cos(angle + geometry.angleToDetectorPixel[i]);
		float y_det = geometry.distanceToDetectorPixel[i] * std::sin(angle + geometry.angleToDetectorPixel[i]);
		outLineImg[i] = sumLine(Line(x_source, y_source, x_det, y_det));
	}
	return outLineImg;
}

void Projector::buildFullProjection(){
	// build projection for every angle
	// x axis is angle to detector as an index in input 'angles' vector
	// y axis is a detector pixel

	std::unique_ptr<float[]> outfloatImg(new float[geometry.detectorCount * geometry.angles.size()]);
	for (int i = 0; i < geometry.angles.size(); ++i) {
		float angle = geometry.angles[i];
		std::unique_ptr<float[]> lineImg = getSingleProjection(angle);
		for (int j = 0; j < geometry.detectorCount; ++j) {
			outfloatImg[i * geometry.detectorCount + j ] = lineImg[j];
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

	std::unique_ptr<unsigned char[]> outCharImg(new unsigned char[geometry.detectorCount * geometry.angles.size()]);
	float max_val = *std::max_element(&fullProjection[0], &fullProjection[geometry.detectorCount * geometry.angles.size() - 1]);
	// check for too small max value; if (max < 1.0) then max = 1.0
	float reverse_max = (max_val > 1.0f) ? (1 / max_val) : 1.0f;
	for (int i = 0; i < geometry.angles.size(); ++i) {
		for (int j = 0; j < geometry.detectorCount; ++j) {
			outCharImg[i * geometry.detectorCount + j] = (unsigned char)(fullProjection[i * geometry.detectorCount + j] * 255.0f * reverse_max);
		}
	}
	return outCharImg;
};