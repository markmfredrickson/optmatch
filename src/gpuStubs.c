#include<R.h>

void printNoGPUMessage() {
	warning("\nDuring install of source package, CUDA_HOME was not defined.\nPlease check that the CUDA toolkit is properly installed and\nreinstall the package from source. In the meantime,\ncomputation with the GPU will be unavailable.\n");
}

void gpuMahalanobisHelper(const int * vectorsSize, const int * vectorLength,
	double * vectors1, double * vectors2, const double * mat,
	double * result)
{
	printNoGPUMessage();
}
