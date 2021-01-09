#pragma once
#ifndef POS2ECF_H
#define POS2ECF_H
#define PI       3.14159265358979323846   // pi
#include "Eigen/Dense"
#include "Eigen/Core"
#include "iostream"
#include "iomanip"
#include "assert.h"
#include "string"
#include "fstream"
class POS2ECF
{
public:
	explicit POS2ECF(Eigen::Vector3d  wgs84Point, Eigen::Vector3d nedPos):
		_wgs84Point(wgs84Point), _nedPos(nedPos){}
	~POS2ECF() {}

	void Wgs84Point2EcefPoint();
	void NedAngle2EcefAngle();

	Eigen::Vector3d GetEcefPoint();
	Eigen::Vector3d GetEcefPos();
	void Caculate();
private:
	Eigen::Vector3d _wgs84Point;						//(omega, alpha, H) means  longitude, latitude, ellipsoidal height
	Eigen::Vector3d _nedPos;									//(heading, pitch, roll) in NED

	Eigen::Vector3d _ecefPoint;								//(X, Y, Z)
	Eigen::Vector3d _ecefPos;								//(psi, theta, phi) in ECEF
};

bool ReadPosData(std::string posAbsPath, Eigen::Vector3d& wgs84Point, Eigen::Vector3d& nedPos);
bool WriteEcefData(std::string ecefAbsPath, Eigen::Vector3d& ecefPoint, Eigen::Vector3d& ecefPos);
Eigen::Matrix3d AngleAxisRotationMatrix(Eigen::Vector3d axis, double angle);
#endif
