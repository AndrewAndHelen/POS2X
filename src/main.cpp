#include "pos2ecf.hpp"
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

	std::cout << "--- read pos data---" << std::endl;
	bool flag = ReadPosData<Eigen::Vector3d>(srcPosPath, wgs84Point, nedPos);

	std::cout << "--- pos2ecf ---" << std::endl;
	POS2ECF<double> pos2ecf(wgs84Point, nedPos);
	pos2ecf.Caculate();
	Eigen::Vector3d ecefPoint = pos2ecf.GetEcefPoint();
	Eigen::Vector3d ecefPos = pos2ecf.GetEcefPos();

	std::cout << "--- write pos data ---" << std::endl;
	flag = WriteEcefData<Eigen::Vector3d>(objPosPath, ecefPoint, ecefPos);

	return 0;
}