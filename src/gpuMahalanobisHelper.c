#include<R.h>

#include<cuda_runtime.h>
#include<cublas_v2.h>

void gpuMahalanobisHelper(const int * vectorsSize, const int * vectorLength,
	double * vectors1, double * vectors2, const double * mat,
	double * result)
{
	int
		nv = *vectorsSize, n = *vectorLength;

	cudaError_t cudaStat;
	double
		* dVectors, * dMat, * dPartialResult, * dResult;

	cudaStat = cudaMalloc((void**)&dVectors, nv * n * sizeof(*vectors1));
	if(cudaStat != cudaSuccess)
		error("gpuMahalanobisHelper:\n\tcannot allocate GPU memory\n");

	cudaStat = cudaMalloc((void**)&dPartialResult, nv * n * sizeof(*vectors2));
	if(cudaStat != cudaSuccess) {
		cudaFree(dVectors);
		error("gpuMahalanobisHelper:\n\tcannot allocate GPU memory\n");
	}

	cudaStat = cudaMalloc((void**)&dMat, n * n * sizeof(*dMat));
	if(cudaStat != cudaSuccess) {
		cudaFree(dVectors);
		cudaFree(dPartialResult);
		error("gpuMahalanobisHelper:\n\tcannot allocate GPU memory\n");
	}

	cublasStatus_t stat;
	cublasHandle_t handle;

	stat = cublasCreate(&handle);
	if(stat != CUBLAS_STATUS_SUCCESS) {
		cudaFree(dVectors);
		cudaFree(dPartialResult);
		cudaFree(dMat);
		error("gpuMahalanobisHelper:\n\tGPU API CUBLAS init failed\n");
	}

	stat = cublasSetMatrix(nv, n, sizeof(*vectors2), vectors1, nv,
		dPartialResult, nv);
	if(stat != CUBLAS_STATUS_SUCCESS) {
		cudaFree(dVectors);
		cudaFree(dPartialResult);
		cudaFree(dMat);
		cublasDestroy(handle);
		error("gpuMahalanobisHelper:\n\tcannot copy to GPU memory\n");
	}

	stat = cublasSetMatrix(nv, n, sizeof(*vectors1), vectors2, nv,
		dVectors, nv);
	if(stat != CUBLAS_STATUS_SUCCESS) {
		cudaFree(dVectors);
		cudaFree(dPartialResult);
		cudaFree(dMat);
		cublasDestroy(handle);
		error("gpuMahalanobisHelper:\n\tcannot copy to GPU memory\n");
	}

	stat = cublasSetMatrix(n, n, sizeof(*mat), mat, n, dMat, n);
	if(stat != CUBLAS_STATUS_SUCCESS) {
		cudaFree(dVectors);
		cudaFree(dPartialResult);
		cudaFree(dMat);
		cublasDestroy(handle);
		error("gpuMahalanobisHelper:\n\tcannot copy to GPU memory\n");
	}

	double
		alpha = -1.0, beta = 0.0;

	stat = cublasDaxpy(handle, nv * n, &alpha, dPartialResult, 1, dVectors, 1);
	if(stat != CUBLAS_STATUS_SUCCESS) {
		cudaFree(dVectors);
		cudaFree(dPartialResult);
		cudaFree(dMat);
		cublasDestroy(handle);
		error("gpuMahalanobisHelper:\n\tcannot perform GPU computation\n");
	}

	alpha *= -1.0;
	stat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nv, n, n, &alpha,
		dVectors, nv, dMat, n, &beta, dPartialResult, nv);
	if(stat != CUBLAS_STATUS_SUCCESS) {
		cudaFree(dVectors);
		cudaFree(dPartialResult);
		cudaFree(dMat);
		cublasDestroy(handle);
		error("gpuMahalanobisHelper:\n\tcannot perform GPU computation\n");
	}

	stat = cublasGetMatrix(nv, n, sizeof(*dVectors), dVectors, nv,
		vectors1, nv);
	if(stat != CUBLAS_STATUS_SUCCESS) {
		cudaFree(dPartialResult);
		cublasDestroy(handle);
		error("gpuMahalanobisHelper:\n\tcannot perform GPU computation\n");
	}

	stat = cublasGetMatrix(nv, n, sizeof(*dPartialResult), dPartialResult, nv,
		vectors2, nv);
	if(stat != CUBLAS_STATUS_SUCCESS) {
		cudaFree(dPartialResult);
		cublasDestroy(handle);
		error("gpuMahalanobisHelper:\n\tcannot perform GPU computation\n");
	}

	cudaFree(dMat);
	cudaFree(dVectors);
	cudaFree(dPartialResult);
	cublasDestroy(handle);

	int i, j;
	for(i = 0; i < nv * n; i++)
		vectors1[i] *= vectors2[i];

	double sum;
	for(i = 0; i < nv; i++) {
		sum = 0.0;
		for(j = 0; j < n; j++)
			sum += vectors1[i + j * nv];
		result[i] = sqrt(sum);
	}
}
