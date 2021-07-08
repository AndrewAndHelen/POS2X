#ifdef  _WIN32
#include<io.h>
#include<direct.h>
#include <cstdio>
#else defined linux
#include <sys/io.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

#include "pos2x.hpp"

static void Usage(const char* pszErrorMsg = NULL)
{
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "pos2ecf -srcPosPath *** -objPosPath ***\n");
	fprintf(stderr, "options:\n");
	fprintf(stderr, "[-help,-h]											[produce help message]\n");
	fprintf(stderr, "[-srcPosPath]									[the path of origin pos data]\n");
	fprintf(stderr, "[-objPosPath]									[the path of object pos data]\n");

	if (pszErrorMsg != NULL)
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);

	exit(1);
}

int main(int argc, char* argv[])
{
	std::string srcPosPath = "";
	std::string objPosPath = "";

	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "-h") == 0)
		{
			Usage();
		}
		else if (strcmp(argv[i], "-srcPosPath") == 0)
		{
			i++; if (i >= argc) continue;
			srcPosPath = argv[i];
		}
		else if (strcmp(argv[i], "-objPosPath") == 0) {
			i++; if (i >= argc) continue;
			objPosPath = argv[i];
		}
		else
		{
			Usage("Too many command options.");
		}
	}

#ifdef _WIN32
	if (0 != _access(srcPosPath.c_str(), 0))
	{
		// if this file not exist, error and quit.
		std::cout << "No such input pos path";
		exit(1);
	}
#else defined linux
	if (0 != eaccess(srcPosPath.c_str(), F_OK))
	{
		// if this file not exist, error and quit.
		std::cout << "No such input pos path";
		exit(1);
	}
#endif

	//case 1: WGS84 longitude, latitude, ellipsoidal height, heading pitch roll(NED) --->WGS84 X, Y, Z, heading pitch roll(ECEF)  
	Eigen::Vector3d wgs84Point = Eigen::Vector3d::Zero();
	Eigen::Vector3d nedPos = Eigen::Vector3d::Zero();
	double ECEF_X, ECEF_Y, ECEF_Z;
	double ecef_heading, ecef_pitch, ecef_roll;

	std::cout << "--- read wgs84 pos data---" << std::endl;
	ReadPosData<Eigen::Vector3d>(srcPosPath, wgs84Point, nedPos);

	std::cout << "--- pos2ecef ---" << std::endl;
	POSTRANSTOOL<double> postrans1;

	/*double _X, _Y, _Z, _phi, _omega, _kappa;
	postrans1.WGS84YPR_2_CGCS2000UTMOPK(121.468669454296, 38.8316363124972, 1434.30909345299,
		-68.8018892701904 * M_PI / 180.0, -88.3451156024876 * M_PI / 180.0, -156.862875188552 * M_PI / 180.0,
		_X, _Y, _Z, _omega, _phi, _kappa);*/

	postrans1.LONLATALLT2ECEF(wgs84Point[0], wgs84Point[1], wgs84Point[2], ECEF_X, ECEF_Y, ECEF_Z);

	postrans1.WGS84_HPR2ECEF_HPR(nedPos[0], nedPos[1], nedPos[2],wgs84Point[0], wgs84Point[1], ecef_heading, ecef_pitch, ecef_roll);

	std::cout << "--- write ecef pos data ---" << std::endl;
	Eigen::Vector3d ecefPoint(ECEF_X, ECEF_Y, ECEF_Z);
	Eigen::Vector3d ecefPos(ecef_heading, ecef_pitch, ecef_roll);
	WriteEcefData<Eigen::Vector3d>(objPosPath, ecefPoint, ecefPos);

	return 0;
}