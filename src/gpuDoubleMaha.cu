#include<R.h>

#define NTHREADS 64

void safeCudaMallocDouble(double ** vect, int ndoubles) {
	cudaError_t err = cudaSuccess;
	err = cudaMalloc((void **) vect, ndoubles * sizeof(double));
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

size_t safeCudaMallocPitchDouble(double ** vect, int width, int height) {
	cudaError_t err = cudaSuccess;
	size_t pitch;
	err = cudaMallocPitch(vect, &pitch, width * sizeof(double), height);
	if(err != cudaSuccess) error(cudaGetErrorString(err));
	return pitch;
}

void safeToDeviceDouble(double * a, const double * b, int ndoubles) {
	cudaError_t err = cudaSuccess;
	err = cudaMemcpy(a, b, ndoubles * sizeof(double), cudaMemcpyHostToDevice);
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

void safeToDevice2DDouble(double * a, size_t aPitch,
	const double * b, size_t nrows, size_t ncols)
{
	cudaError_t err = cudaSuccess;
	err = cudaMemcpy2D(a, aPitch, b, ncols * sizeof(double),
		ncols * sizeof(double), nrows, cudaMemcpyHostToDevice);
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

void safeFromDeviceDouble(double * a, const double * b, int ndoubles) {
	cudaError_t err = cudaSuccess;
	err = cudaMemcpy(a, b, ndoubles * sizeof(double), cudaMemcpyDeviceToHost);
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

void checkCudaErrorDouble() {
	cudaError_t err = cudaSuccess;
	cudaGetLastError();
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

void safeCudaFreeDouble(double * a) {
	cudaError_t err = cudaSuccess;
	err = cudaFree(a);
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

__global__ void dDoubleMaha(int nvectors, int n,
	const double * vectors1, size_t v1Pitch,
	const double * vectors2, size_t v2Pitch,
	const double * mat, size_t matPitch,
	double * result)
{
	int k, i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i >= nvectors) return;

	double
		sum = 0, innerSum;
	const double
		* matCol,
		* v1i = vectors1 + i * v1Pitch,
		* v2i = vectors2 + i * v2Pitch;

	for(int j = 0; j < n; j++) {
		innerSum = 0;
		matCol = mat + j * matPitch;
		for(k = 0; k < n; k++)
			innerSum += (v1i[k] - v2i[k]) * matCol[k];
		sum += innerSum * (v1i[j] - v2i[j]);
	}
	result[i] = sqrt(sum);
}

extern "C"
void gpuDoubleMaha(const int * vectorSetSize, const int * vectorLength,
	const double * vectorSet1, const double * vectorSet2, const double * mat,
	double * result)
{
	int
		nv = *vectorSetSize, n = *vectorLength;

	double
		* dMat = NULL, * dVectorSet1 = NULL, * dVectorSet2 = NULL,
		* dResult = NULL;
	size_t
		v1Pitch, v2Pitch, matPitch;

	safeCudaMallocDouble(&dResult, nv);

	v1Pitch = safeCudaMallocPitchDouble(&dVectorSet1, n, nv);
	v2Pitch = safeCudaMallocPitchDouble(&dVectorSet2, n, nv);
	matPitch = safeCudaMallocPitchDouble(&dMat, n, nv);

	safeToDevice2DDouble(dVectorSet1, v1Pitch, vectorSet1, nv, n);
	safeToDevice2DDouble(dVectorSet2, v2Pitch, vectorSet2, nv, n);
	safeToDevice2DDouble(dMat, matPitch, mat, n, n);

	size_t nblocks = ceil((double) nv /(double) NTHREADS);

	dDoubleMaha<<<nblocks, NTHREADS>>>(nv, n,
		dVectorSet1, v1Pitch / sizeof(double),
		dVectorSet2, v2Pitch / sizeof(double),
		dMat, matPitch / sizeof(double), dResult);
	checkCudaErrorDouble();
	safeFromDeviceDouble(result, dResult, nv);

	safeCudaFreeDouble(dVectorSet1);
	safeCudaFreeDouble(dVectorSet2);
	safeCudaFreeDouble(dMat);
	safeCudaFreeDouble(dResult);
}
