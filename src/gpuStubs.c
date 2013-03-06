#include<R.h>

void throwNoGPUError() {
	error("GPU-enabled features not available.\n\tCUDA_HOME was not set at install time.\n\tPlease check that CUDA is properly installed,\n\tand set CUDA_HOME to the location of the CUDA toolkit\n\tprior to installing the optmatch package from source.\n");
}

void gpuMahalanobisHelper(const int * vectorSetSize, const int * vectorLength,
	const float * vectorSet1, const float * vectorSet2, const float * mat,
	float * result)
{
	throwNoGPUError();
}

void gpuDoubleMaha(const int * vectorSetSize, const int * vectorLength,
	const double * vectorSet1, const double * vectorSet2, const double * mat,
	double * result)
{
	throwNoGPUError();
}
