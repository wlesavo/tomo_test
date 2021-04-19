#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
#include <vector>
#include <cmath>
class Geometry {
public:
	const double dSourceToObj, dObjToDetector, detectorSize, imgCenterX, imgCenterY;
	const int detectorCount;
	const std::vector<double> angles;
	std::vector<double> distanceToDetectorPixelCenter, angleToDetectorPixelCenter;
	std::vector<double> distanceToDetectorPixel, angleToDetectorPixel;

	Geometry(const std::vector<double>& i_angles, int i_detectorCount, double i_detectorSize, double i_dSourceToObj
		, double i_dObjToDetector, double imgCenterX, double imgCenterY)
		: dSourceToObj(i_dSourceToObj), dObjToDetector(i_dObjToDetector), detectorSize(i_detectorSize)
		, detectorCount(i_detectorCount), angles(i_angles), imgCenterX(imgCenterX), imgCenterY(imgCenterY)
	{};
	
	virtual Line v_GetNextLineCenter(int angle, int detector_i) = 0;
	virtual Line v_GetNextLineEdge(int angle, int detector_i) = 0;
	
};
#endif