#pragma once
#ifndef POS2ECF_HPP
#define POS2ECF_HPP
#define PI       3.14159265358979323846   // pi
#include "Eigen/Dense"
#include "Eigen/Core"
#include "iostream"
#include "iomanip"
#include "assert.h"
#include "string"
#include "fstream"

template<typename DataType>
class POS2ECF
{
public:
	typedef Eigen::Matrix<DataType, 3, 3> Matrix3X;
	typedef Eigen::Matrix<DataType, 3, 1> Vector3X;

	explicit POS2ECF(Vector3X  wgs84Point, Vector3X nedPos) :
		_wgs84Point(wgs84Point), _nedPos(nedPos) {}
	~POS2ECF() {}

	void Wgs84Point2EcefPoint();
	void NedAngle2EcefAngle();

	Matrix3X AngleAxisRotationMatrix(Vector3X axis, DataType angle);
	Vector3X GetEcefPoint();
	Vector3X GetEcefPos();
	void Caculate();
private:
	Vector3X _wgs84Point;						//(omega, alpha, H) means  longitude, latitude, ellipsoidal height
	Vector3X _nedPos;									//(heading, pitch, roll) in NED

	Vector3X _ecefPoint;								//(X, Y, Z)
	Vector3X _ecefPos;								//(psi, theta, phi) in ECEF
};
#endif

