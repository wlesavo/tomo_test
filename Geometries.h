#ifndef _GEOMETRIES_H_
#define _GEOMETRIES_H_
#include <cmath>
#include <vector>
#include <corecrt_math_defines.h>

class GeometryFanBeam : public Geometry {
public:
	GeometryFanBeam(const std::vector<double>& i_angles, int i_detectorCount, double i_detectorSize, double i_dSourceToObj, double i_dObjToDetector, double imgCenterX, double imgCenterY) :
		Geometry(i_angles, i_detectorCount, i_detectorSize, i_dSourceToObj, i_dObjToDetector, imgCenterX, imgCenterY)
	{
		// precalculating distances and angles for pixels of detector matrix
		// in experiment pixels are rotating around origin with new distance and initial angle

		double center = (detectorCount) * 0.5 * detectorSize;
		for (int i = 0; i < detectorCount; ++i) {
			double det_y = center - detectorSize * (i + 0.5);
			double dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5);
			double angle = std::atan2(-det_y, dObjToDetector);
			distanceToDetectorPixelCenter.push_back(dist);
			angleToDetectorPixelCenter.push_back(angle);
		}
		for (int i = 0; i < detectorCount + 1; ++i) {
			double det_y = center - detectorSize * (i);
			double dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5);
			double angle = std::atan2(-det_y, dObjToDetector);
			distanceToDetectorPixel.push_back(dist);
			angleToDetectorPixel.push_back(angle);
		}
	}

	virtual Line v_GetNextLineCenter(int i_angle, int detector_i) {
		double angle = angles[i_angle];
		double x_source = -dSourceToObj * std::sin(angle) + imgCenterX;
		double y_source = -dSourceToObj * std::cos(angle) + imgCenterY;
		double x_det = distanceToDetectorPixelCenter[detector_i] * std::sin(angle - angleToDetectorPixelCenter[detector_i]) + imgCenterX;
		double y_det = distanceToDetectorPixelCenter[detector_i] * std::cos(angle - angleToDetectorPixelCenter[detector_i]) + imgCenterY;
		return Line(x_source, y_source, x_det, y_det);
	};

	virtual Line v_GetNextLineEdge(int i_angle, int detector_i) {
		double angle = angles[i_angle];
		double x_source = -dSourceToObj * std::sin(angle) + imgCenterX;
		double y_source = -dSourceToObj * std::cos(angle) + imgCenterY;
		double x_det = distanceToDetectorPixel[detector_i] * std::sin(angle - angleToDetectorPixel[detector_i]) + imgCenterX;
		double y_det = distanceToDetectorPixel[detector_i] * std::cos(angle - angleToDetectorPixel[detector_i]) + imgCenterY;
		return Line(x_source, y_source, x_det, y_det);
	};

};

class GeometryParallel : public Geometry {
public:
	GeometryParallel(const std::vector<double>& i_angles, int i_detectorCount, double i_detectorSize, double imgCenterX, double imgCenterY) :
		Geometry(i_angles, i_detectorCount, i_detectorSize, 0, 0, imgCenterX, imgCenterY)
	{
		// precalculating distances and angles for pixels of detector matrix
		// in experiment pixels are rotating around origin with new distance and initial angle

		double center = (detectorCount) * 0.5f * detectorSize;
		for (int i = 0; i < detectorCount; ++i) {
			double det_y = detectorSize * (i + 0.5f) - center;
			distanceToDetectorPixelCenter.push_back(det_y);
		}

		for (int i = 0; i < detectorCount + 1; ++i) {
			double det_y = detectorSize * (i) - center;
			distanceToDetectorPixel.push_back(det_y);
		}
	}

	virtual Line v_GetNextLineCenter(int i_angle, int detector_i) {
		double angle = angles[i_angle];
		double x_det = distanceToDetectorPixelCenter[detector_i] * std::cos(-angle) + imgCenterX;
		double y_det = distanceToDetectorPixelCenter[detector_i] * std::sin(-angle) + imgCenterY;
		return Line(x_det, y_det, M_PI_2 - angle);
	};

	virtual Line v_GetNextLineEdge(int i_angle, int detector_i) {
		double angle = angles[i_angle];
		double x_det = distanceToDetectorPixel[detector_i] * std::cos(-angle) + imgCenterX;
		double y_det = distanceToDetectorPixel[detector_i] * std::sin(-angle) + imgCenterY;
		return Line(x_det, y_det, M_PI_2 - angle);
	};

};

#endif