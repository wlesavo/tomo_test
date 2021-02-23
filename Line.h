#ifndef _LINE_HPP_
#define _LINE_HPP_
#include <vector>

class Line {
public:
	bool vert = false, hor = false;
	float dist = 0.0f, angle = 0.0f, k = 0.0f, b = 0.0f;
	Line(const float& x1, const float& y1, const float& x2, const float& y2);
	float value(const float& x) const;
	float coor(const float& y) const;
};

#endif