#ifndef POS2X_HPP
#define POS2X_HPP
#include "Eigen/Dense"
#include "Eigen/Core"
#include <iostream>
#include <iomanip>
#include <corecrt_math_defines.h>
#include <string>
#include <fstream>

template<typename T>
inline T Deg2Rad(T deg)
{
	return deg * M_PI / 180;
}

template<typename T>
inline T Rad2Deg(T rad)
{
	return rad * 180 / M_PI;
}

template<typename T>
inline T round1(T input)
{
	T result;
	if (input > 1) result = 1;
	else if (input < -1) result = -1;
	else result = input;
	return result;
}

template<typename Vector3X>
bool ReadPosData(std::string posAbsPath, Vector3X& wgs84Point, Vector3X& nedPos)
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

template<typename Vector3X>
bool WriteEcefData(std::string ecefAbsPath, Vector3X& ecefPoint, Vector3X& ecefPos)
{
	std::ofstream fp(ecefAbsPath.c_str());
	if (!fp)
	{
		std::cout << "couldn't create such file!";
		return false;
	}
	fp << std::setiosflags(std::ios::fixed) << std::setprecision(7) << ecefPoint[0] << "," << ecefPoint[1] << "," << ecefPoint[2] << ","
		<< ecefPos[0] << "," << ecefPos[1] << "," << ecefPos[2];

	fp.close();
	return true;
}

template<typename DataType>
class POSTRANSTOOL
{
public:
	typedef Eigen::Matrix<DataType, 3, 3> Matrix3X;
	typedef Eigen::Matrix<DataType, 3, 1> Vector3X;

	POSTRANSTOOL() {}
	~POSTRANSTOOL() {}

	Matrix3X AngleAxisRotationMatrix(Vector3X axis, DataType angle);

	Matrix3X POK2Matrix(DataType phi, DataType omega, DataType kappa);

	Matrix3X YPR2Matrix(DataType yaw, DataType pitch, DataType roll);

	void R2YPR(Matrix3X R, DataType& yaw, DataType& pitch, DataType& roll);

	void R2OPK(Matrix3X R, DataType& omega, DataType& phi, DataType& kappa);

	void R2OPK2(Matrix3X R, DataType& omega, DataType& phi, DataType& kappa);

	DataType ArcLengthOfMeridian(DataType sm_a, DataType sm_b, DataType phi);

	void MapLatLonToXY(DataType sm_a, DataType sm_b, DataType phi, DataType lambda, DataType lambda0, DataType& x, DataType& y);

	void LONLATALLT2ECEF(DataType longitude, DataType latitude, DataType altitude, DataType& X, DataType& Y, DataType& Z);

	void ECEF2LONLATALLT(DataType X, DataType Y, DataType Z, DataType& longitude, DataType& latitude, DataType& altitude);

	void WGS84_HPR2ECEF_HPR(DataType wgs84_heading, DataType wgs84_pitch, DataType wgs84_roll, DataType wgs84_lon, DataType wgs84_lat, DataType& ecef_heading, DataType& ecef_pitch, DataType& ecef_roll);

	void WGS84_HPR2ECEF_YPR(DataType wgs84_heading, DataType wgs84_pitch, DataType wgs84_roll, DataType wgs84_lon, DataType wgs84_lat, DataType& ecef_yaw, DataType& ecef_pitch, DataType& ecef_roll);

	void LatLonDegreeToGauXY(DataType latitude, DataType longitude, int nType, int ParaDegree, DataType& x, DataType& y);

	void LatLonDegreeToUTMXY(DataType lat, DataType lon, int nType, int ParaDegree, DataType& x, DataType& y);

	void WGS84OPK_2_ECEFOPK(DataType wgs84_omega, DataType wgs84_phi, DataType wgs84_kappa, DataType wgs84_lon, DataType wgs84_lat, DataType& ecef_omega, DataType& ecef_phi, DataType& ecef_kappa);

	void WGS84YPR_2_CGCS2000GAUSSOPK(
		const DataType lon, const DataType lat, const DataType alt, const DataType yaw, const DataType pitch, const DataType roll,
		DataType& X, DataType& Y, DataType& Z, DataType& omega, DataType& phi, DataType& kappa);// 经纬度 yaw pitch roll 转为cgcs 2000, 高斯投影的omega phi kappa

	void WGS84YPR_2_CGCS2000UTMOPK(
		const DataType lon, const DataType lat, const DataType alt, const DataType yaw, const DataType pitch, const DataType roll,
		DataType& X, DataType& Y, DataType& Z, DataType& omega, DataType& phi, DataType& kappa);// 经纬度 yaw pitch roll 转为cgcs 2000, 墨卡托投影的omega phi kappa

	void POK_SMART3D_2_POK_PHOTOGRAMMETRY(DataType phi1_deg, DataType omega1_deg, DataType kappa1_deg, DataType& phi_rad, DataType& omega_rad, DataType& kappa_rad);
};


template<typename DataType>
Eigen::Matrix<DataType, 3, 3> POSTRANSTOOL<DataType>::AngleAxisRotationMatrix(Vector3X axis, DataType angle)
{
	Eigen::Matrix<DataType, 3, 3> nX = Eigen::Matrix<DataType, 3, 3>::Zero();
	Eigen::Matrix<DataType, 3, 3> I3 = Eigen::Matrix<DataType, 3, 3>::Identity();
	DataType n1 = axis[0]; DataType n2 = axis[1]; DataType n3 = axis[2];

	nX(0, 1) = -n3;
	nX(0, 2) = n2;
	nX(1, 2) = -n1;
	nX(1, 0) = -nX(0, 1);
	nX(2, 0) = -nX(0, 2);
	nX(2, 1) = -nX(1, 2);

	Eigen::Matrix<DataType, 3, 3> R = (1 - std::cos(angle)) * axis * axis.transpose() + std::cos(angle) * I3 + std::sin(angle) * nX;

	return R;
}

template<typename DataType>
Eigen::Matrix<DataType, 3, 3> POSTRANSTOOL<DataType>::POK2Matrix(DataType phi, DataType omega, DataType kappa)
{
	DataType radPhi = Deg2Rad<DataType>(phi);
	DataType radOmega = Deg2Rad<DataType>(omega);
	DataType radKappa = Deg2Rad<DataType>(kappa);

	DataType cp = std::cos(radPhi);
	DataType sp = std::sin(radPhi);
	DataType co = std::cos(radOmega);
	DataType so = std::sin(radOmega);
	DataType ck = std::cos(radKappa);
	DataType sk = std::sin(radKappa);

	Eigen::Matrix<DataType, 3, 3> R = Eigen::Matrix<DataType, 3, 3>::Zero();
	R(0, 0) = cp * ck;
	R(0, 1) = co * sk + so * sp * ck;
	R(0, 2) = so * sk - co * sp * ck;

	R(1, 0) = -cp * sk;
	R(1, 1) = co * ck - so * sp * sk;
	R(1, 2) = so * ck + co * sp * sk;

	R(2, 0) = sp;
	R(2, 1) = -so * cp;
	R(2, 2) = co * cp;

	return R;
}

template<typename DataType>
Eigen::Matrix<DataType, 3, 3> POSTRANSTOOL<DataType>::YPR2Matrix(DataType yaw, DataType pitch, DataType roll)
{
	DataType radYaw = Deg2Rad<DataType>(yaw);
	DataType radPitch = Deg2Rad<DataType>(pitch);
	DataType radRoll = Deg2Rad<DataType>(roll);

	DataType cr = cos(radRoll);
	DataType cy = cos(radYaw);
	DataType sr = sin(radRoll);
	DataType sp = sin(radPitch);
	DataType sy = sin(radYaw);
	DataType cp = cos(radPitch);

	Eigen::Matrix<DataType, 3, 3> R = Eigen::Matrix<DataType, 3, 3>::Zero();

	R(0, 0) = cr * cy - sr * sp * sy;
	R(0, 1) = -cr * sy - cy * sr * sp;
	R(0, 2) = cp * sr;

	R(1, 0) = cy * sr + cr * sp * sy;
	R(1, 1) = cr * cy * sp - sr * sy;
	R(1, 2) = -cr * cp;

	R(2, 0) = cp * sy;
	R(2, 1) = cp * cy;
	R(2, 2) = sp;

	return R;
}

template<typename DataType>
void POSTRANSTOOL<DataType>::R2YPR(Matrix3X R, DataType& yaw, DataType& pitch, DataType& roll)
{
	DataType phi_rad, theta_rad, psi_rad;

	if (fabs(R(2,2)) != 1)
	{
		theta_rad = asin(R(2,2));
		double cos_theta = cos(theta_rad);
		psi_rad = atan2(R(2,0) / cos_theta, R(2,1) / cos_theta);
		phi_rad = atan2(R(0,2) / cos_theta, -R(1,2) / cos_theta);
	}
	else
	{
		phi_rad = 0;
		if (R(2, 2) == -1)
		{
			theta_rad = -PI / 2;
			psi_rad = atan2(R(0,1), R(0,0));
		}
		else
		{
			theta_rad = PI / 2;
			psi_rad = atan2(R(0, 1), R(0, 0));
		}
	}

	yaw = psi_rad;
	pitch = theta_rad;
	roll = phi_rad;
}

template<typename DataType>
void POSTRANSTOOL<DataType>::R2OPK(Matrix3X R, DataType& omega, DataType& phi, DataType& kappa)
{
	phi = std::asin(R(2,0));
	omega = std::asin(-R(2,1) / cos(phi));
	kappa = std::asin(-R(1,0) / cos(phi));

	DataType kappa2 = std::acos(R(0,0) / cos(phi));

	kappa *= 180 / PI;
	kappa2 *= 180 / PI;

	//计算kappa角的范围，应该是-180~180
	if (kappa2 < 90) { ; }
	else if (kappa > 0) { kappa = kappa2; }
	else { kappa = -kappa2; }

	kappa *= PI / 180;

	DataType omega2 = std::acos(R(2,2) / cos(phi));
	omega *= 180 / PI;
	omega2 *= 180 / PI;

	//计算omega角的范围，应该是-180~180
	if (omega2 < 90) { ; }
	else if (omega > 0) { omega = omega2; }
	else { omega = -omega2; }

	omega *= PI / 180;
}

template<typename DataType>
void POSTRANSTOOL<DataType>::R2OPK2(Matrix3X R, DataType& omega, DataType& phi, DataType& kappa)
{
	omega = asin(-round1(R(1,2)));
	phi = asin(-round1(R(0,2) / cos(omega)));
	kappa = asin(round1(R(1,0) / cos(omega)));
	DataType kappa2 = acos(round1(R(1,1) / cos(omega)));

	phi = Rad2Deg<DataType>(phi);
	omega = Rad2Deg<DataType>(omega);
	kappa = Rad2Deg<DataType>(kappa);
	kappa2 = Rad2Deg<DataType>(kappa2);

	//计算kappa角的范围，应该是-180~180
	if (kappa2 < 90) { ; }
	else if (kappa > 0) { kappa = kappa2; }
	else { kappa = -kappa2; }

	phi = Deg2Rad<DataType>(phi);
	omega = Deg2Rad<DataType>(omega);
	kappa = Deg2Rad<DataType>(kappa);
}

template<typename DataType>
DataType POSTRANSTOOL<DataType>::ArcLengthOfMeridian(DataType sm_a, DataType sm_b, DataType phi)
{
	DataType alpha, beta, gamma, delta, epsilon, n;
	DataType result;
	n = (sm_a - sm_b) / (sm_a + sm_b);
	alpha = ((sm_a + sm_b) / 2.0) * (1.0 + (pow(n, 2.0) / 4.0) + (pow(n, 4.0) / 64.0));
	beta = (-3.0 * n / 2.0) + (9.0 * pow(n, 3.0) / 16.0) + (-3.0 * pow(n, 5.0) / 32.0);
	gamma = (15.0 * pow(n, 2.0) / 16.0) + (-15.0 * pow(n, 4.0) / 32.0);
	delta = (-35.0 * pow(n, 3.0) / 48.0) + (105.0 * pow(n, 5.0) / 256.0);
	epsilon = (315.0 * pow(n, 4.0) / 512.0);
	result = alpha * (phi + (beta * sin(2.0 * phi)) + (gamma * sin(4.0 * phi)) + (delta * sin(6.0 * phi)) + (epsilon * sin(8.0 * phi)));

	return result;
}

template<typename DataType>
void POSTRANSTOOL<DataType>::MapLatLonToXY(DataType sm_a, DataType sm_b, DataType phi, DataType lambda, DataType lambda0, DataType& x, DataType& y)
{
	DataType N, nu2, ep2, t, t2, l;
	DataType l3coef, l4coef, l5coef, l6coef, l7coef, l8coef;
	DataType tmp;

	ep2 = (pow(sm_a, 2.0) - pow(sm_b, 2.0)) / pow(sm_b, 2.0);
	nu2 = ep2 * pow(cos(phi), 2.0);
	N = pow(sm_a, 2.0) / (sm_b * sqrt(1 + nu2));
	t = tan(phi);
	t2 = t * t;
	tmp = (t2 * t2 * t2) - pow(t, 6.0);
	l = lambda - lambda0;
	l3coef = 1.0 - t2 + nu2;
	l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
	l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2 - 58.0 * t2 * nu2;
	l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2 - 330.0 * t2 * nu2;
	l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);
	l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);
	x = N * cos(phi) * l + (N / 6.0 * pow(cos(phi), 3.0) * l3coef * pow(l, 3.0))
		+ (N / 120.0 * pow(cos(phi), 5.0) * l5coef * pow(l, 5.0))
		+ (N / 5040.0 * pow(cos(phi), 7.0) * l7coef * pow(l, 7.0));
	y = ArcLengthOfMeridian(sm_a, sm_b, phi)
		+ (t / 2.0 * N * pow(cos(phi), 2.0) * pow(l, 2.0))
		+ (t / 24.0 * N * pow(cos(phi), 4.0) * l4coef * pow(l, 4.0))
		+ (t / 720.0 * N * pow(cos(phi), 6.0) * l6coef * pow(l, 6.0))
		+ (t / 40320.0 * N * pow(cos(phi), 8.0) * l8coef * pow(l, 8.0));
}

template<typename DataType>
void POSTRANSTOOL<DataType>::LONLATALLT2ECEF(DataType longitude, DataType latitude, DataType altitude, DataType& X, DataType& Y, DataType& Z)
{
	/* Default parameters*/
	double a = 6378137;										//Semi-major
	double b = 6356752.3142;							//Semi-minor

	DataType domega = static_cast<DataType>(Deg2Rad<double>(longitude));
	DataType dalpha = static_cast<DataType>(Deg2Rad<double>(latitude));
	DataType H = static_cast<DataType>(altitude);

	DataType Rn = static_cast<DataType>(a / std::sqrt(std::pow(std::cos(dalpha), 2) + std::pow(b, 2) / std::pow(a, 2) * std::pow(std::sin(dalpha), 2)));
	DataType Re = static_cast<DataType>(b / std::sqrt(std::pow(a, 2) / std::pow(b, 2) * std::pow(std::cos(dalpha), 2) + std::pow(std::sin(dalpha), 2)));

	X = static_cast<DataType>((Rn + H) * std::cos(dalpha) * std::cos(domega));
	Y = static_cast<DataType>((Rn + H) * std::cos(dalpha) * std::sin(domega));
	Z = static_cast<DataType>((Re + H) * std::sin(dalpha));
}

template<typename DataType>
void POSTRANSTOOL<DataType>::ECEF2LONLATALLT(DataType X, DataType Y, DataType Z, DataType& longitude, DataType& latitude, DataType& altitude)
{
	/* Default parameters*/
	double a = 6378137;									//Semi-major
	double b = 6356752.31424518;						//Semi-minor

	double f = 1.0 / 298.257223563;

	double e = sqrtl(((a * a) - (b * b)) / (a * a));

	double e1 = sqrtl(((a * a) - (b * b)) / (b * b));

	double p = sqrtl(static_cast<double>(X) * static_cast<double>(X) + static_cast<double>(Y) * static_cast<double>(Y));

	double q = atan2l((static_cast<double>(Z) * a), (p * b));

	longitude = atan2l(static_cast<double>(Y), static_cast<double>(X));

	latitude = atan2l((static_cast<double>(Z) + (e1 * e1) * b * powl(sinl(q), 3)), (p - (e * e) * a * powl(cosl(q), 3)));

	double N = a / sqrtl(1 - ((e * e) * powl(sinl(latitude), 2)));

	altitude = (p / cosl(latitude)) - N;

	longitude = Rad2Deg<DataType>(longitude);

	latitude = Rad2Deg<DataType>(latitude);
}

template<typename DataType>
void POSTRANSTOOL<DataType>::WGS84_HPR2ECEF_HPR(DataType wgs84_heading, DataType wgs84_pitch, DataType wgs84_roll, DataType wgs84_lon, DataType wgs84_lat, DataType& ecef_heading, DataType& ecef_pitch, DataType& ecef_roll)
{
	DataType domega = Deg2Rad<DataType>(wgs84_lon);
	DataType dalpha = Deg2Rad<DataType>(wgs84_lat);

	DataType heading = Deg2Rad<DataType>(wgs84_heading);
	DataType pitch = Deg2Rad<DataType>(wgs84_pitch);
	DataType roll = Deg2Rad<DataType>(wgs84_roll);

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

	DataType psi = static_cast<DataType>(std::atan2(X3.dot(Y0), X3.transpose() * X0));
	DataType theta = static_cast<DataType>(std::atan2(-X3.dot(Z0), std::sqrt(std::pow(X3.dot(X0), 2) + std::pow(X3.dot(Y0), 2))));

	Matrix3X RZ0psi = AngleAxisRotationMatrix(Z0, psi);
	Y2 = RZ0psi * Y0;
	Matrix3X RY2theta = AngleAxisRotationMatrix(Y2, theta);
	Z2 = RY2theta * Z0;

	DataType phi = static_cast<DataType>(std::atan2(Y3.dot(Z2), Y3.dot(Y2)));

	ecef_heading = psi * 180 / M_PI;
	ecef_pitch = theta * 180 / M_PI;
	ecef_roll = phi * 180 / M_PI;
}

template<typename DataType>
void POSTRANSTOOL<DataType>::WGS84_HPR2ECEF_YPR(DataType wgs84_heading, DataType wgs84_pitch, DataType wgs84_roll, DataType wgs84_lon, DataType wgs84_lat, DataType& ecef_yaw, DataType& ecef_pitch, DataType& ecef_roll)
{
	double domega = Angle2Radians<double>(wgs84_lon);
	double dalpha = Angle2Radians<double>(wgs84_lat);

	double heading = Angle2Radians<double>(wgs84_heading);
	double pitch = Angle2Radians<double>(wgs84_pitch);
	double roll = -Angle2Radians<double>(wgs84_roll);

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

	double psi = std::atan2(X3.dot(Y0), X3.transpose() * X0);
	double theta = std::atan2(-X3.dot(Z0), std::sqrt(std::pow(X3.dot(X0), 2) + std::pow(X3.dot(Y0), 2)));

	Eigen::Matrix3d RZ0psi = AngleAxisRotationMatrix(Z0, psi);
	Y2 = RZ0psi * Y0;
	Eigen::Matrix3d RY2theta = AngleAxisRotationMatrix(Y2, theta);
	Z2 = RY2theta * Z0;

	double phi = std::atan2(Y3.dot(Z2), Y3.dot(Y2));

	ecef_yaw = 90.f - psi * 180 / M_PI;
	ecef_pitch = -theta * 180 / M_PI;
	ecef_roll = -180.f - phi * 180 / M_PI;
}

template<typename DataType>
void POSTRANSTOOL<DataType>::LatLonDegreeToGauXY(DataType latitude, DataType longitude, int nType, int ParaDegree, DataType& x, DataType& y)
{
	double a = 6378137.0;
	double f = 1.0 / 298.257223;
	if (0 == nType)  //北京54
	{
		a = 6378245.0;
		f = 1.0 / 298.3;
	}
	else if (1 == nType) //西安80
	{
		a = 6378137.0;
		f = 1.0 / 298.257222101;
	}
	else if (2 == nType) //wgs84
	{
		a = 6378137.0;
		f = 1.0 / 298.257223;
	}
	else if (3 == nType) //CGCS2000
	{
		a = 6378137.0;
		f = 1 / 298.257222101;
	}

	//	longitude = Dms2Rad(longitude);
	//	latitude = Dms2Rad(latitude);
	int nLo = 1, nLa = 1;
	if (longitude < 0)
	{
		// 		longitude = fabs(longitude);
		// 		nLo = -1;
		longitude += 360;
	}
	if (latitude < 0)
	{
		latitude = fabs(latitude);
		nLa = -1;
	}

	//	latitude = 0;
	int ProjNo = 0;

	if (6 == ParaDegree)
	{
		ProjNo = (int)(longitude / ParaDegree) + 1;
	}
	else if (3 == ParaDegree)
	{
		//	longitude = 1.5;
		ProjNo = (int)((longitude - 1.5) / ParaDegree) + 1;
		if ((longitude >= 358.5 && longitude <= 360) || (longitude >= 0 && longitude <= 1.5))
		{
			ProjNo = 120;
		}
	}

	DataType longitude1, latitude1, longitude0 = 0.0, latitude0, x0, y0, xval, yval;
	DataType e2, ee, NN, T, C, A, M, iPI;
	iPI = PI / 180.0;
	if (6 == ParaDegree)
	{
		longitude0 = ProjNo * ParaDegree - ParaDegree / 2;
	}
	else if (3 == ParaDegree)
	{
		longitude0 = ProjNo * 3;
		if ((longitude >= 0 && longitude <= 1.5))
		{
			longitude0 = 0;
		}
	}

	longitude0 = longitude0 * iPI;
	latitude0 = 0;
	longitude1 = longitude * iPI;
	latitude1 = latitude * iPI;
	e2 = 2 * f - f * f;
	ee = e2 * (1.0 - e2);
	NN = a / sqrt(1.0 - e2 * sin(latitude1) * sin(latitude1));
	T = tan(latitude1) * tan(latitude1);
	C = ee * cos(latitude1) * cos(latitude1);
	A = (longitude1 - longitude0) * cos(latitude1);

	M = a * ((1 - e2 / 4 - 3 * e2 * e2 / 64 - 5 * e2 * e2 * e2 / 256) * latitude1 - (3 * e2 / 8 + 3 * e2 * e2 / 32 + 45 * e2 * e2
		* e2 / 1024) * sin(2 * latitude1)
		+ (15 * e2 * e2 / 256 + 45 * e2 * e2 * e2 / 1024) * sin(4 * latitude1) - (35 * e2 * e2 * e2 / 3072) * sin(6 * latitude1));
	yval = NN * (A + (1 - T + C) * A * A * A / 6 + (5 - 18 * T + T * T + 72 * C - 58 * ee) * A * A * A * A * A / 120);
	xval = M + NN * tan(latitude1) * (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * A * A * A * A / 24
		+ (61 - 58 * T + T * T + 600 * C - 330 * ee) * A * A * A * A * A * A / 720);
	//y0 = 1000000L*(ProjNo)+500000L; //带号+平移
	y0 = 500000L; //只平移
	x0 = 0;
	xval = xval + x0;
	yval = yval + y0;

	x = yval * nLo;
	y = xval * nLa;
}

template<typename DataType>
void POSTRANSTOOL<DataType>::LatLonDegreeToUTMXY(DataType lat, DataType lon, int nType, int ParaDegree, DataType& x, DataType& y)
{
	//ParaDegree = 6; //UTM投影一般只有6度分带

	double UTMScaleFactor = 0.9996;
	double sm_a = 6378137.0;
	double sm_b = 6356752.31414;
	double sm_EccSquared = 6.69437999013e-03;

	if (0 == nType)  //北京54
	{
		sm_a = 6378245.0;
		sm_b = 6356863.0188;
		sm_EccSquared = 0.006693421614543;
	}
	else if (1 == nType) //西安80
	{
		sm_a = 6378137.0;
		sm_b = 6356755.2882;
		sm_EccSquared = 0.00669438499959;
	}
	else if (2 == nType) //wgs84
	{
		sm_a = 6378137.0;
		sm_b = 6356752.31414;
		sm_EccSquared = 6.69437999013e-03;
	}
	else if (3 == nType) //CGCS2000
	{
		sm_a = 6378137.0;
		sm_b = 6356752.31414;
		sm_EccSquared = 0.006694380022898;
	}

	//计算投影带中心
	int zone = 0;
	DataType m = 0;

	if (ParaDegree == 6)
	{
		zone = static_cast<int>((lon + 180.0) / 6) + 1;
		lon = lon * M_PI / 180;
		lat = lat * M_PI / 180;
		m = (-183.0 + (zone * 6.0)) * M_PI / 180;
	}

	DataType _x, _y;
	MapLatLonToXY(sm_a, sm_b, lat, lon, m, _x, _y);

	_x = _x * UTMScaleFactor + 500000.0;
	_y = _y * UTMScaleFactor;
	if (_y < 0.0) //南半球
		_y += 10000000.0;

	x = _x;
	y = _y;
}

template<typename DataType>
void POSTRANSTOOL<DataType>::WGS84OPK_2_ECEFOPK(DataType wgs84_omega, DataType wgs84_phi, DataType wgs84_kappa, DataType wgs84_lon, DataType wgs84_lat, DataType& ecef_omega, DataType& ecef_phi, DataType& ecef_kappa)
{
	DataType wgs84_omega_rad = Deg2Rad<DataType>(wgs84_omega);
	DataType wgs84_phi_rad = Deg2Rad<DataType>(wgs84_phi);
	DataType wgs84_kappa_rad = -Deg2Rad<DataType>(wgs84_kappa);

	Matrix3X WGS_R = POK2Matrix(wgs84_phi, wgs84_omega, wgs84_kappa);
	DataType YPR_yaw, YPR_pitch, YPR_roll;

	R2YPR(WGS_R, YPR_yaw, YPR_pitch, YPR_roll);

	//std::cout << YPR_yaw * 180 / PI << " " << YPR_pitch * 180 / PI << " " << YPR_roll * 180 / PI << " " << std::endl;
	DataType domega = Deg2Rad<DataType>(wgs84_lon);
	DataType dalpha = Deg2Rad<DataType>(wgs84_lat);

	DataType heading = YPR_yaw;
	DataType pitch = YPR_pitch;
	DataType roll = -YPR_roll;

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

	DataType psi = static_cast<DataType>(std::atan2(X3.dot(Y0), X3.transpose() * X0));
	DataType theta = static_cast<DataType>(std::atan2(-X3.dot(Z0), std::sqrt(std::pow(X3.dot(X0), 2) + std::pow(X3.dot(Y0), 2))));

	Matrix3X RZ0psi = AngleAxisRotationMatrix(Z0, psi);
	Y2 = RZ0psi * Y0;
	Matrix3X RY2theta = AngleAxisRotationMatrix(Y2, theta);
	Z2 = RY2theta * Z0;

	DataType phi = static_cast<DataType>(std::atan2(Y3.dot(Z2), Y3.dot(Y2)));

	DataType ecef_yaw = 90 - Rad2Deg<DataType>(psi);
	DataType ecef_pitch = -Rad2Deg<DataType>(theta);
	DataType ecef_roll = -180 - Rad2Deg<DataType>(phi);

	Matrix3X ecef_r = YPR2Matrix(ecef_yaw, ecef_pitch, ecef_roll);
	R2OPK(ecef_r, ecef_omega, ecef_phi, ecef_kappa);

	ecef_omega = Rad2Deg<DataType>(ecef_omega);
	ecef_phi = Rad2Deg<DataType>(ecef_phi);
	ecef_kappa = Rad2Deg<DataType>(ecef_kappa);
}


template<typename DataType>
void POSTRANSTOOL<DataType>::WGS84YPR_2_CGCS2000GAUSSOPK(
	const DataType lon, const DataType lat, const DataType alt, const DataType yaw, const DataType pitch, const DataType roll,
	DataType& X, DataType& Y, DataType& Z, DataType& omega, DataType& phi, DataType& kappa)
{
	// lla to ecef
	DataType ecef_x, ecef_y, ecef_z;
	LONLATALLT2ECEF(lon, lat, alt, ecef_x, ecef_y, ecef_z);

	// ENU-N
	DataType n_x, n_y, n_z;

	DataType theta_geodetic = Deg2Rad<DataType>(lon);
	DataType phi_geoetic = Deg2Rad<DataType>(lat);
	n_x = -cos(theta_geodetic) * sin(phi_geoetic);
	n_y = -sin(theta_geodetic) * sin(phi_geoetic);
	n_z = cos(phi_geoetic);

	// ENU-N +1
	DataType ecef_2_x = ecef_x + n_x;
	DataType ecef_2_y = ecef_y + n_y;
	DataType ecef_2_z = ecef_z + n_z;

	// ecef to lla
	DataType lon_2, lat_2, alt_2;
	ECEF2LONLATALLT(ecef_2_x, ecef_2_y, ecef_2_z, lon_2, lat_2, alt_2);

	// lla to gause
	double gau_1_x, gau_1_y;
	LatLonDegreeToGauXY(lat, lon, 3, 6, gau_1_x, gau_1_y);

	X = gau_1_x;
	Y = gau_1_y;
	Z = alt;

	DataType gau_2_x, gau_2_y;
	LatLonDegreeToGauXY(lat_2, lon_2, 3, 6, gau_2_x, gau_2_y);

	DataType v_x = gau_2_x - gau_1_x;
	DataType v_y = gau_2_y - gau_1_y;
	DataType xxyy = sqrt(v_x * v_x + v_y * v_y);
	v_x = v_x / xxyy;
	v_y = v_y / xxyy;
	DataType angle = acos(v_y);

	DataType yaw_2 = yaw + angle;
	Matrix3X R = YPR2Matrix(yaw_2 * 180 / M_PI, pitch * 180 / M_PI, roll * 180 / M_PI);

	Matrix3X coef;
	coef << 1, 0, 0,
		0, -1, 0,
		0, 0, 1;

	Matrix3X R_ = coef * R;

	R = R_.transpose();

	R2OPK2(R, omega, phi, kappa);
}

template<typename DataType>
void POSTRANSTOOL<DataType>::WGS84YPR_2_CGCS2000UTMOPK(
	const DataType lon, const DataType lat, const DataType alt, const DataType yaw, const DataType pitch, const DataType roll,
	DataType& X, DataType& Y, DataType& Z, DataType& omega, DataType& phi, DataType& kappa)
{
	// lla to ecef
	DataType ecef_x, ecef_y, ecef_z;
	LONLATALLT2ECEF(lon, lat, alt, ecef_x, ecef_y, ecef_z);

	// ENU-N
	DataType n_x, n_y, n_z;

	DataType theta_geodetic = Deg2Rad<DataType>(lon);
	DataType phi_geoetic = Deg2Rad<DataType>(lat);
	n_x = -cos(theta_geodetic) * sin(phi_geoetic);
	n_y = -sin(theta_geodetic) * sin(phi_geoetic);
	n_z = cos(phi_geoetic);

	// ENU-N +1
	DataType ecef_2_x = ecef_x + n_x;
	DataType ecef_2_y = ecef_y + n_y;
	DataType ecef_2_z = ecef_z + n_z;

	// ecef to lla
	DataType lon_2, lat_2, alt_2;
	ECEF2LONLATALLT(ecef_2_x, ecef_2_y, ecef_2_z, lon_2, lat_2, alt_2);

	// lla to UTM
	DataType gau_1_x, gau_1_y;
	LatLonDegreeToUTMXY(lat, lon, 3, 6, gau_1_x, gau_1_y);

	X = gau_1_x;
	Y = gau_1_y;
	Z = alt;

	DataType gau_2_x, gau_2_y;
	LatLonDegreeToUTMXY(lat_2, lon_2, 3, 6, gau_2_x, gau_2_y);

	DataType v_x = gau_2_x - gau_1_x;
	DataType v_y = gau_2_y - gau_1_y;
	DataType xxyy = sqrt(v_x * v_x + v_y * v_y);
	v_x = v_x / xxyy;
	v_y = v_y / xxyy;
	DataType angle = acos(v_y);

	DataType yaw_2 = yaw + angle;

	Matrix3X R = YPR2Matrix(yaw_2 * 180 / M_PI, pitch * 180 / M_PI, roll * 180 / M_PI);

	Matrix3X coef;
	coef << 1, 0, 0,
		0, -1, 0,
		0, 0, 1;

	Matrix3X R_ = coef * R;

	R = R_.transpose();

	R2OPK2(R, omega, phi, kappa);
}

template<typename DataType>
void  POSTRANSTOOL<DataType>::POK_SMART3D_2_POK_PHOTOGRAMMETRY(DataType phi1_deg, DataType omega1_deg, DataType kappa1_deg, DataType& phi_rad, DataType& omega_rad, DataType& kappa_rad)
{
	Matrix3X R = POK2Matrix(phi1_deg, omega1_deg, kappa1_deg);

	Matrix3X coef;
	coef << 1, 0, 0,
		0, -1, 0,
		0, 0, 1;

	Matrix3X R_ = coef * R;

	R = R_.transpose();

	R2OPK2(R, omega_rad, phi_rad, kappa_rad);
}
#endif

