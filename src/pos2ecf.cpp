#include "pos2ecf.h"

template<typename T>
inline T Angle2Radians(T angle)
{
	return angle * PI / 180;
}

bool ReadPosData(std::string posAbsPath, Eigen::Vector3d& wgs84Point, Eigen::Vector3d& nedPos)
{
	std::fstream fp(posAbsPath.c_str());
	if (!fp)
	{
		std::cout << "couldn't open file!";
		return false;
	}
	std::string strLine;
	for (int i = 0; i < 3; i++)
	{
		std::getline(fp, strLine, ',');
		wgs84Point[i] = atof(strLine.c_str());
	}
	for (int i = 0; i < 3; i++)
	{
		std::getline(fp, strLine, ',');
		nedPos[i] = atof(strLine.c_str());
	}

	fp.close();
	return true;
}

bool WriteEcefData(std::string ecefAbsPath, Eigen::Vector3d& ecefPoint, Eigen::Vector3d& ecefPos)
{
	std::ofstream fp(ecefAbsPath.c_str());
	if (!fp)
	{
		std::cout << "couldn't create such file!";
		return false;
	}
	fp << std::setiosflags(std::ios::fixed) << std::setprecision(7) << ecefPoint[0] << "," << ecefPoint[1] << "," << ecefPoint[2]<<","
		<< ecefPos[0] << "," << ecefPos[1] << "," << ecefPos[2];
	
	fp.close();
	return true;
}
/*
* Summary: Get rotation matrix beyond rotate axis and angle
* Parameters:
*		@param axis:unit vector of axis.
*		@param theta: rotate angle.

*		@return rotation matrix
*/
Eigen::Matrix3d AngleAxisRotationMatrix(Eigen::Vector3d axis, double theta)
{
	Eigen::Matrix3d nX = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
	double n1 = axis[0];double n2 = axis[1];double n3 = axis[2];
	
	nX(0, 1) = -n3;
	nX(0, 2) = n2;
	nX(1, 2) = -n1;
	nX(1, 0) = -nX(0, 1);
	nX(2, 0) = -nX(0, 2);
	nX(2, 1) = -nX(1, 2);

	Eigen::Matrix3d Rn = (1 - std::cos(theta))*axis*axis.transpose() + std::cos(theta)*I3 + std::sin(theta)*nX;

	//Eigen::Matrix3d RnM = Rn.transpose()*Rn;
	//Eigen::ArrayXXf  RnA = RnM.array().pow(-0.5);//if has 0 element ,error.
	//Rn = Rn * RnA.matrix();
	return Rn;
}

/*
* Summary: This fuction is to caculate X, Y, Z in ECEF coordnate system from WGS84 coordnate system.
* Parameters:
*		Variables:omega->longitude, alpha->latitude, H->ellipsoidal height.
*Tips:
*The formula is:
*		X = (Rn + H)*cos(alpha)*cos(omega).
*		Y = (Rn + H)*cos(alpha)*sin(omega).
*		Z = [Rn(1 - e^2) + h]*sin(alpha),  this Rn(1 - e^2) will replace with Re, but result is same.
*/
void POS2ECF::Wgs84Point2EcefPoint()
{
	/* Default parameters*/
	double a = 6378137;										//Semi-major
	double b = 6356752.3142;							//Semi-minor

	_ecefPoint = Eigen::Vector3d::Zero();
	double domega = Angle2Radians<double>( _wgs84Point[0]);
	double dalpha = Angle2Radians<double>(_wgs84Point[1]);
	double H = _wgs84Point(2, 0);

	assert(a / std::sqrt(std::pow(std::cos(dalpha), 2) + std::pow(b, 2) / std::pow(a, 2) * std::pow(std::sin(dalpha), 2)));					//Rn
	assert(b / std::sqrt(std::pow(a, 2) / std::pow(b, 2)*std::pow(std::cos(dalpha), 2) + std::pow(std::sin(dalpha), 2)));					//Re

	double Rn = a / std::sqrt(std::pow(std::cos(dalpha), 2) + std::pow(b, 2) / std::pow(a, 2) * std::pow(std::sin(dalpha), 2));
	double Re = b / std::sqrt(std::pow(a, 2) / std::pow(b, 2)*std::pow(std::cos(dalpha), 2) + std::pow(std::sin(dalpha), 2));

	double X = (Rn + H)*std::cos(dalpha)*std::cos(domega);
	double Y = (Rn + H)*std::cos(dalpha)*std::sin(domega);
	double Z = (Re + H)*std::sin(dalpha);

	_ecefPoint[0] = X;
	_ecefPoint[1] = Y;
	_ecefPoint[2] = Z;
}
void POS2ECF::NedAngle2EcefAngle()
{
	double domega = Angle2Radians<double>(_wgs84Point[0]);
	double dalpha = Angle2Radians<double>(_wgs84Point[1]);

	double heading = Angle2Radians<double>(_nedPos[0]);
	double pitch = Angle2Radians<double>(_nedPos[1]);
	double roll = Angle2Radians<double>(_nedPos[2]);

	Eigen::Vector3d N0(0, 0, 1);
	Eigen::Vector3d E0(0, 1, 0);
	Eigen::Vector3d D0(-1, 0, 0);
	
	Eigen::Matrix3d RN0Omega = AngleAxisRotationMatrix(N0, domega);
	Eigen::Vector3d N1 = N0;
	Eigen::Vector3d E1 = RN0Omega * E0;
	Eigen::Vector3d D1 = RN0Omega * D0;
	
	Eigen::Matrix3d R_E1Alpha = AngleAxisRotationMatrix(-E1, dalpha);
	Eigen::Vector3d E = E1;
	Eigen::Vector3d N = R_E1Alpha * N1;
	Eigen::Vector3d D = N.cross(E);
	/* Norm the N E D*/
	//E.normalize(); N.normalize(); D.normalize();

	Eigen::Vector3d X0(1, 0, 0);
	Eigen::Vector3d Y0(0, 1, 0);
	Eigen::Vector3d Z0(0, 0, 1);

	Eigen::Matrix3d RDHeading = AngleAxisRotationMatrix(D, heading);
	Eigen::Vector3d X1 = RDHeading * N;
	Eigen::Vector3d Y1 = RDHeading * E;
	Eigen::Vector3d Z1 = RDHeading * D;

	Eigen::Matrix3d REPitch = AngleAxisRotationMatrix(Y1, pitch);
	Eigen::Vector3d X2 = REPitch * X1;
	Eigen::Vector3d Y2 = REPitch * Y1;
	Eigen::Vector3d Z2 = REPitch * Z1;

	Eigen::Matrix3d RNRoll = AngleAxisRotationMatrix(X2, roll);
	Eigen::Vector3d X3 = RNRoll * X2;
	Eigen::Vector3d Y3 = RNRoll * Y2;
	Eigen::Vector3d Z3 = RNRoll * Z2;
	
	double psi = std::atan2(X3.dot(Y0), X3.transpose()*X0);
	double theta = std::atan2(-X3.dot(Z0), std::sqrt(std::pow(X3.dot(X0), 2) + std::pow(X3.dot(Y0), 2)));

	Eigen::Matrix3d RZ0psi = AngleAxisRotationMatrix(Z0, psi);
	Y2 = RZ0psi * Y0;
	Eigen::Matrix3d RY2theta = AngleAxisRotationMatrix(Y2, theta);
	Z2 = RY2theta * Z0;

	double phi = std::atan2(Y3.dot(Z2), Y3.dot(Y2));

	_ecefPos[0] = psi * 180 / PI;
	_ecefPos[1] = theta * 180 / PI;
	_ecefPos[2] = phi * 180 / PI;
}

Eigen::Vector3d POS2ECF::GetEcefPoint()
{
	if (!_ecefPoint.isZero())
		return _ecefPoint;
	else
	{
		std::cout << "empty Vector!\n";
		return Eigen::Vector3d::Zero();
	}
}

Eigen::Vector3d POS2ECF::GetEcefPos()
{
	if (!_ecefPos.isZero())
		return _ecefPos;
	else
	{
		std::cout << "empty Vector!\n";
		return Eigen::Vector3d::Zero();
	}
}

void POS2ECF::Caculate()
{
	Wgs84Point2EcefPoint();
	NedAngle2EcefAngle();
}