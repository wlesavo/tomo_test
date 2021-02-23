#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
#include <vector>
#include <cmath>
class Geometry {
public:
	const float dSourceToObj, dObjToDetector, detectorSize;
	const int detectorCount;
	const std::vector<float> angles;
	std::vector<float> distanceToDetectorPixel, angleToDetectorPixel;

	Geometry(const std::vector<float>& i_angles, int i_detectorCount, float i_detectorSize, float i_dSourceToObj, float i_dObjToDetector) :
		dSourceToObj(i_dSourceToObj), dObjToDetector(i_dObjToDetector), detectorSize(i_detectorSize), detectorCount(i_detectorCount), angles(i_angles)
	{
		// precalculating distances and angles for pixels of detector matrix
		// in experiment pixels are rotating around origin with new distance and initial angle

		float center = (detectorCount) * 0.5f * detectorSize;
		for (int i = 0; i < detectorCount; ++i){
			float det_y = detectorSize * i + 0.5f - center;
			float dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5f);
			float angle = std::atan2(-det_y, dObjToDetector);
			distanceToDetectorPixel.push_back(dist);
			angleToDetectorPixel.push_back(angle);
		}
	};
	Geometry(const Geometry& geom) : dSourceToObj(geom.dSourceToObj), dObjToDetector(geom.dObjToDetector)
		, detectorSize(geom.detectorSize), detectorCount(geom.detectorCount), angles(geom.angles)
		, distanceToDetectorPixel(geom.distanceToDetectorPixel), angleToDetectorPixel(geom.angleToDetectorPixel)
	{};
};
#endif