#include "Line.h"
#include <cmath>

Line::Line(const float& x1, const float& y1, const float& x2, const float& y2) {
	angle = std::atan2(y2 - y1, x2 - x1);
	if (std::abs(std::cos(angle)) < 1e-9f) {
		vert = true;
		dist = x1;
	}
	else if (std::abs(std::cos(angle)) > 1.0f - 1e-9f) {
		hor = true;
		dist = y1;
	}
	else {
		k = std::tan(angle);
		b = y1 - k * x1;
	}
}

float Line::value(const float& x) const {
	return x * k + b;
}

float Line::coor(const float& y) const {
	return ((y - b) / k);
};
