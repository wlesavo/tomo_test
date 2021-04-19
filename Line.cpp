#include "Line.h"
#include <cmath>

Line::Line(double x1, double y1, double x2, double y2) {
	angle = std::atan2(y2 - y1, x2 - x1);
	if (std::abs(std::cos(angle)) < 1e-10) {
		vert = true;
		b = x1;
	}
	else if (std::abs(std::sin(angle)) < 1e-10) {
		hor = true;
		k = 0;
		b = y1;
	}
	else {
		k = (y2 - y1) / (x2 - x1);// std::tan(angle);
		b = y2 - k * x2;
	}
}

Line::Line(double x, double y, double i_angle) {
	angle = i_angle;
	if (std::abs(std::cos(angle)) < 1e-10) {
		vert = true;
		b = x;
	}
	else if (std::abs(std::sin(angle)) < 1e-10) {
		hor = true;
		k = 0;
		b = y;
	}
	else {
		k = std::tan(angle);
		b = y - k * x;
	}
}

double Line::value(double x) const {
	return x * k + b;
}

double Line::coor(double y) const {
	double rev = 1 / k;
	return ((y - b) * rev);
};
