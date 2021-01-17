#ifndef POS2ECF__IMPL_HPP
#define POS2ECF_IMPL_HPP
#include "pos2ecf.hpp"

template<typename T>
inline T Angle2Radians(T angle)
{
	return angle * PI / 180;
}

/*
* Summary: Get rotation matrix beyond rotate axis and angle
* Parameters:
*		@param axis:unit vector of axis.
*		@param theta: rotate angle.

*		@return rotation matrix
*/
template<typename DataType>
Eigen::Matrix<DataType, 3, 3> POS2ECF<DataType>::AngleAxisRotationMatrix(Vector3X axis, DataType angle)
{
	Matrix3X nX = Matrix3X::Zero();
	Matrix3X I3 = Matrix3X::Identity();
	DataType n1 = axis[0]; DataType n2 = axis[1]; DataType n3 = axis[2];

	nX(0, 1) = -n3;
	nX(0, 2) = n2;
	nX(1, 2) = -n1;
	nX(1, 0) = -nX(0, 1);
	nX(2, 0) = -nX(0, 2);
	nX(2, 1) = -nX(1, 2);

	Matrix3X Rn = (1 - std::cos(angle))*axis*axis.transpose() + std::cos(angle)*I3 + std::sin(angle)*nX;

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
template<typename DataType>
void POS2ECF<DataType>::Wgs84Point2EcefPoint()
{
	/* Default parameters*/
	double a = 6378137;										//Semi-major
	double b = 6356752.3142;							//Semi-minor

	_ecefPoint = Vector3X::Zero();
	DataType domega = static_cast<DataType>(Angle2Radians<double>(_wgs84Point[0]));
	DataType dalpha = static_cast<DataType>(Angle2Radians<double>(_wgs84Point[1]));
	DataType H = static_cast<DataType>(_wgs84Point(2, 0));

	assert(a / std::sqrt(std::pow(std::cos(dalpha), 2) + std::pow(b, 2) / std::pow(a, 2) * std::pow(std::sin(dalpha), 2)));					//Rn
	assert(b / std::sqrt(std::pow(a, 2) / std::pow(b, 2)*std::pow(std::cos(dalpha), 2) + std::pow(std::sin(dalpha), 2)));					//Re

	DataType Rn = static_cast<DataType>(a / std::sqrt(std::pow(std::cos(dalpha), 2) + std::pow(b, 2) / std::pow(a, 2) * std::pow(std::sin(dalpha), 2)));
	DataType Re = static_cast<DataType>(b / std::sqrt(std::pow(a, 2) / std::pow(b, 2)*std::pow(std::cos(dalpha), 2) + std::pow(std::sin(dalpha), 2)));

	DataType X = static_cast<DataType>((Rn + H)*std::cos(dalpha)*std::cos(domega));
	DataType Y = static_cast<DataType>((Rn + H)*std::cos(dalpha)*std::sin(domega));
	DataType Z = static_cast<DataType>((Re + H)*std::sin(dalpha));

	_ecefPoint[0] = X;
	_ecefPoint[1] = Y;
	_ecefPoint[2] = Z;
}

template<typename DataType>
void POS2ECF<DataType>::NedAngle2EcefAngle()
{
	DataType domega = Angle2Radians<DataType>(_wgs84Point[0]);
	DataType dalpha = Angle2Radians<DataType>(_wgs84Point[1]);

	DataType heading = Angle2Radians<DataType>(_nedPos[0]);
	DataType pitch = Angle2Radians<DataType>(_nedPos[1]);
	DataType roll = Angle2Radians<DataType>(_nedPos[2]);

	Vector3X N0(0, 0, 1);
	Vector3X E0(0, 1, 0);
	Vector3X D0(-1, 0, 0);

	Matrix3X RN0Omega = AngleAxisRotationMatrix(N0, domega);
	Vector3X N1 = N0;
	Vector3X E1 = RN0Omega * E0;
	Vector3X D1 = RN0Omega * D0;

	Matrix3X R_E1Alpha = AngleAxisRotationMatrix(-E1, dalpha);
	Vector3X E = E1;
	Vector3X N = R_E1Alpha * N1;
	Vector3X D = N.cross(E);
	/* Norm the N E D*/
	//E.normalize(); N.normalize(); D.normalize();

	Vector3X X0(1, 0, 0);
	Vector3X Y0(0, 1, 0);
	Vector3X Z0(0, 0, 1);

	Matrix3X RDHeading = AngleAxisRotationMatrix(D, heading);
	Vector3X X1 = RDHeading * N;
	Vector3X Y1 = RDHeading * E;
	Vector3X Z1 = RDHeading * D;

	Matrix3X REPitch = AngleAxisRotationMatrix(Y1, pitch);
	Vector3X X2 = REPitch * X1;
	Vector3X Y2 = REPitch * Y1;
	Vector3X Z2 = REPitch * Z1;

	Matrix3X RNRoll = AngleAxisRotationMatrix(X2, roll);
	Vector3X X3 = RNRoll * X2;
	Vector3X Y3 = RNRoll * Y2;
	Vector3X Z3 = RNRoll * Z2;

	DataType psi = static_cast<DataType>(std::atan2(X3.dot(Y0), X3.transpose()*X0));
	DataType theta = static_cast<DataType>(std::atan2(-X3.dot(Z0), std::sqrt(std::pow(X3.dot(X0), 2) + std::pow(X3.dot(Y0), 2))));

	Matrix3X RZ0psi = AngleAxisRotationMatrix(Z0, psi);
	Y2 = RZ0psi * Y0;
	Matrix3X RY2theta = AngleAxisRotationMatrix(Y2, theta);
	Z2 = RY2theta * Z0;

	DataType phi = static_cast<DataType>(std::atan2(Y3.dot(Z2), Y3.dot(Y2)));

	_ecefPos[0] = psi * 180 / PI;
	_ecefPos[1] = theta * 180 / PI;
	_ecefPos[2] = phi * 180 / PI;
}

template<typename DataType>
Eigen::Matrix<DataType, 3, 1> POS2ECF<DataType>::GetEcefPoint()
{
	if (!_ecefPoint.isZero())
		return _ecefPoint;
	else
	{
		std::cout << "empty Vector!\n";
		return Vector3X::Zero();
	}
}

template<typename DataType>
Eigen::Matrix<DataType, 3, 1> POS2ECF<DataType>::GetEcefPos()
{
	if (!_ecefPos.isZero())
		return _ecefPos;
	else
	{
		std::cout << "empty Vector!\n";
		return Vector3X::Zero();
	}
}

template<typename DataType>
void POS2ECF<DataType>::Caculate()
{
	Wgs84Point2EcefPoint();
	NedAngle2EcefAngle();
}
#endif
