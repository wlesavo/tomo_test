#ifndef _LINE_HPP_
#define _LINE_HPP_
#include <vector>

class Line {
public:
	bool vert = false, hor = false;
	double angle = 0.0, k = 0.0, b = 0.0;
	Line(double x1, double y1, double x2, double y2);
	Line(double x, double y, double angle);
	Line() {};
	double value(double x) const;
	double coor(double y) const;
};

#endif