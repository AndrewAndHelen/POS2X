#include "pos2ecf.h"
#ifdef  _WIN32
#include<io.h>
#include<direct.h>
#else defined linux
#include <sys/io.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

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

	Eigen::Vector3d wgs84Point = Eigen::Vector3d::Zero();
	Eigen::Vector3d nedPos = Eigen::Vector3d::Zero();

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

	ReadPosData(srcPosPath, wgs84Point, nedPos);

	POS2ECF pos2ecf(wgs84Point, nedPos);
	pos2ecf.Caculate();

	Eigen::Vector3d ecefPoint = pos2ecf.GetEcefPoint();
	Eigen::Vector3d ecefPos = pos2ecf.GetEcefPos();
	WriteEcefData(objPosPath, ecefPoint, ecefPos);

	return 0;
}