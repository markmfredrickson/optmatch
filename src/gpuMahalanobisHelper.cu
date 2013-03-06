#include<R.h>

#define NTHREADS 256

texture<float, cudaTextureType2D, cudaReadModeElementType> dMat;

void safeCudaMallocFloat(float ** vect, int nfloats) {
	cudaError_t err = cudaSuccess;
	err = cudaMalloc((void **) vect, nfloats * sizeof(float));
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

size_t safeCudaMallocPitchFloat(float ** vect, int width, int height) {
	cudaError_t err = cudaSuccess;
	size_t pitch;
	err = cudaMallocPitch(vect, &pitch, width * sizeof(float), height);
	if(err != cudaSuccess) error(cudaGetErrorString(err));
	return pitch;
}

void safeToDeviceFloat(float * a, const float * b, int nfloats) {
	cudaError_t err = cudaSuccess;
	err = cudaMemcpy(a, b, nfloats * sizeof(float), cudaMemcpyHostToDevice);
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

void safeToDevice2DFloat(float * a, size_t aPitch,
	const float * b, size_t nrows, size_t ncols)
{
	cudaError_t err = cudaSuccess;
	err = cudaMemcpy2D(a, aPitch, b, ncols * sizeof(float),
		ncols * sizeof(float), nrows, cudaMemcpyHostToDevice);
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

void safeFromDeviceFloat(float * a, const float * b, int nfloats) {
	cudaError_t err = cudaSuccess;
	err = cudaMemcpy(a, b, nfloats * sizeof(float), cudaMemcpyDeviceToHost);
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

void checkCudaError() {
	cudaError_t err = cudaSuccess;
	cudaGetLastError();
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

void safeCudaFree(float * a) {
	cudaError_t err = cudaSuccess;
	err = cudaFree(a);
	if(err != cudaSuccess) error(cudaGetErrorString(err));
}

__global__ void dMahalanobisHelper(int nvectors, int n,
	const float * vectors1, size_t v1Pitch,
	const float * vectors2, size_t v2Pitch,
	float * result)
{
	int k, i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i >= nvectors) return;

	float
		sum = 0, innerSum;
	const float
		* v1i = vectors1 + i * v1Pitch,
		* v2i = vectors2 + i * v2Pitch;

	for(int j = 0; j < n; j++) {
		innerSum = 0;
		for(k = 0; k < n; k++)
			innerSum += (v1i[k] - v2i[k]) * tex2D(dMat, j, k);
		sum += innerSum * (v1i[j] - v2i[j]);
	}
	result[i] = sqrtf(sum);
}

extern "C"
void gpuMahalanobisHelper(const int * vectorSetSize, const int * vectorLength,
	const float * vectorSet1, const float * vectorSet2, const float * mat,
	float * result)
{
	int
		nv = *vectorSetSize, n = *vectorLength;

	float
		* dVectorSet1 = NULL, * dVectorSet2 = NULL, * dResult = NULL;
	size_t v1Pitch, v2Pitch;

	// safeCudaMallocFloat(&dVectorSet1, nv * n);
	// safeCudaMallocFloat(&dVectorSet2, nv * n);

	safeCudaMallocFloat(&dResult, nv);

	v1Pitch = safeCudaMallocPitchFloat(&dVectorSet1, n, nv);
	v2Pitch = safeCudaMallocPitchFloat(&dVectorSet2, n, nv);

	safeToDevice2DFloat(dVectorSet1, v1Pitch, vectorSet1, nv, n);
	safeToDevice2DFloat(dVectorSet2, v2Pitch, vectorSet2, nv, n);

	// create mat in device memory as a read only 2D texture
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0,
		cudaChannelFormatKindFloat);
	cudaArray * cuArray;
	cudaMallocArray(&cuArray, &channelDesc, n, n);
	cudaMemcpyToArray(cuArray, 0, 0, mat, n * n * sizeof(float),
		cudaMemcpyHostToDevice);
	
	dMat.addressMode[0] = cudaAddressModeClamp;
	dMat.addressMode[1] = cudaAddressModeClamp;
	dMat.filterMode = cudaFilterModePoint;
	dMat.normalized = false;

	cudaBindTextureToArray(dMat, cuArray, channelDesc);
	// end create mat

	size_t nblocks = ceil((double) nv /(double) NTHREADS);

	dMahalanobisHelper<<<nblocks, NTHREADS>>>(nv, n,
		dVectorSet1, v1Pitch / sizeof(float),
		dVectorSet2, v2Pitch / sizeof(float),
		dResult);
	checkCudaError();
	safeFromDeviceFloat(result, dResult, nv);

	cudaFreeArray(cuArray);
	safeCudaFree(dVectorSet1);
	safeCudaFree(dVectorSet2);
	// safeCudaFree(dMat);
	safeCudaFree(dResult);
}
