#define NTHREADS 256

__global__ void dMahalanobisHelper(int nvectors, int n,
	const float * vectors1, const float * vectors2, const float * mat,
	float * result)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i >= nvectors) return;

	int
		k, col, innerCol;
	float
		sum = 0, innerSum;

	for(int j = 0; j < n; j++) {
		innerSum = 0;
		col = j * n;
		for(k = 0; k < n; k++) {
			innerCol = i + k * nvectors;
			innerSum += (vectors1[innerCol] - vectors2[innerCol])
				* mat[k + col];
		}
		col = i + j * nvectors;
		sum += innerSum * (vectors1[col] - vectors2[col]);
	}
	result[i] = sqrtf(sum);
}

// computes sqrt((u_i - v_i) * m * (u_i - v_i)) for two sets of vectors
// u and v and a matrix m. The matrix m should have dim n x n where each
// u_i and v_i are of length n. This is intended to fill in the last part of a
// mahabalonis distance calculation for a set of vectors u and v,
// one set a treatment and the other a control.

// The matrix m is expected to be stored column major in memory
// the ith column's elements are stored consecutively at m + i * n
// the jth row will be at m + 0 * n + j, m + n + j, m + 2 * n + j, ...

// each vector is expected to be a row vector of length n, with the whole
// set stored as an nv x n matrix w in column major order :(
// to get at the jth vector in the set, you must index the elements as
// w + 0 * nv + j, w + nv + j, w + 2 * nv + j, ..., w + (n - 1) * nv + j

extern "C"
void gpuMahalanobisHelper(const int * vectorSetSize, const int * vectorLength,
	const float * vectorSet1, const float * vectorSet2, const float * mat,
	float * result)
{
	int
		nv = *vectorSetSize, n = *vectorLength;

	float
		* dVectorSet1, * dVectorSet2, * dMat, * dResult;

	size_t
		resultBytes = nv * sizeof(float),
		matBytes = n * n * sizeof(float),
		vectorSetBytes = nv * n * sizeof(float);

	cudaMalloc((void **) & dVectorSet1, vectorSetBytes);
	cudaMalloc((void **) & dVectorSet2, vectorSetBytes);
	cudaMalloc((void **) & dMat, matBytes);
	cudaMalloc((void **) & dResult, resultBytes);

	cudaMemcpy(dVectorSet1, vectorSet1, vectorSetBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(dVectorSet2, vectorSet2, vectorSetBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(dMat, mat, matBytes, cudaMemcpyHostToDevice);

	dim3 dimBlock(NTHREADS, 1, 1);
	int gx = ceil((double) nv /(double) dimBlock.x);
	dim3 dimGrid(gx, 1, 1);

	cudaDeviceSynchronize();
	dMahalanobisHelper<<<dimBlock, dimGrid>>>(nv, n, dVectorSet1, dVectorSet2,
		dMat, dResult);
	cudaDeviceSynchronize();

	cudaMemcpy(result, dResult, resultBytes, cudaMemcpyDeviceToHost);

	cudaFree(dVectorSet1);
	cudaFree(dVectorSet2);
	cudaFree(dMat);
	cudaFree(dResult);
}

__device__ double ddDotProduct(int n, const double * u, const double * v) {
	double sum = 0.0;
	for(int i = 0; i < n; i++)
		sum += u[i] * v[i];
	return sum;
}

__global__ void ddMahalanobisHelper(int vectorSetSize, int vectorLength,
	const double * vectorSet1, const double * vectorSet2, const double * mat,
	double * temp1, double * temp2, double * result)
{
	int
		i = blockDim.x * blockIdx.x + threadIdx.x,
		nv = vectorSetSize, n = vectorLength;

	if(i >= vectorSetSize) return;

	double
		* temp1i = temp1 + i * vectorLength,
		* temp2i = temp2 + i * vectorLength;

	for(int j = 0; j < vectorLength; j++)
		temp1i[j] = vectorSet1[i + j * nv] - vectorSet2[i + j * nv];
	for(int j = 0; j < vectorLength; j++)
		temp2i[j] = ddDotProduct(n, temp1i, mat + j * n);

	result[i] = sqrtf( ddDotProduct(n, temp1i, temp2i) );
}

// computes sqrt((u_i - v_i) * m * (u_i - v_i)) for two sets of vectors
// u and v and a matrix m. The matrix m should have dim n x n where each
// u_i and v_i are of length n. This is intended to fill in the last part of a
// mahabalonis distance calculation for a set of vectors u and v,
// one set a treatment and the other a control.

// The matrix m is expected to be stored column major in memory
// the ith column's elements are stored consecutively at m + i * n
// the jth row will be at m + 0 * n + j, m + n + j, m + 2 * n + j, ...

// each vector is expected to be a row vector of length n, with the whole
// set stored as an nv x n matrix w in column major order :(
// to get at the jth vector in the set, you must index the elements as
// w + 0 * nv + j, w + nv + j, w + 2 * nv + j, ..., w + (n - 1) * nv + j

extern "C"
void gpu2MahalanobisHelper(const int * vectorSetSize, const int * vectorLength,
	const double * vectorSet1, const double * vectorSet2, const double * mat,
	double * result)
{
	int
		nv = *vectorSetSize, n = *vectorLength;

	double
		* dVectorSet1, * dVectorSet2, * dMat,
		* temp1, * temp2, * dResult;

	size_t
		resultBytes = nv * sizeof(double),
		matBytes = n * n * sizeof(double),
		vectorSetBytes = nv * n * sizeof(double);

	cudaMalloc((void **) & dVectorSet1, vectorSetBytes);
	cudaMalloc((void **) & dVectorSet2, vectorSetBytes);
	cudaMalloc((void **) & dMat, matBytes);
	cudaMalloc((void **) & temp1, vectorSetBytes);
	cudaMalloc((void **) & temp2, vectorSetBytes);
	cudaMalloc((void **) & dResult, resultBytes);

	cudaMemcpy(dVectorSet1, vectorSet1, vectorSetBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(dVectorSet2, vectorSet2, vectorSetBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(dMat, mat, matBytes, cudaMemcpyHostToDevice);

	dim3 dimBlock(NTHREADS, 1, 1);
	int gx = ceil((double) nv /(double) dimBlock.x);
	dim3 dimGrid(gx, 1, 1);

	cudaDeviceSynchronize();
	ddMahalanobisHelper<<<dimBlock, dimGrid>>>(nv, n, dVectorSet1, dVectorSet2,
		dMat, temp1, temp2, dResult);
	cudaDeviceSynchronize();

	cudaFree(temp1);
	cudaFree(temp2);
	cudaFree(dVectorSet1);
	cudaFree(dVectorSet2);
	cudaFree(dMat);

	cudaMemcpy(result, dResult, resultBytes, cudaMemcpyDeviceToHost);
	cudaFree(dResult);
}

__global__ void dMahaDouble(int nvectors, int n,
	const double * vectors1, const double * vectors2, const double * mat,
	double * result)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i >= nvectors) return;

	int
		k, col, innerCol;
	double
		sum = 0, innerSum;

	for(int j = 0; j < n; j++) {
		innerSum = 0;
		col = j * n;
		for(k = 0; k < n; k++) {
			innerCol = i + k * nvectors;
			innerSum += (vectors1[innerCol] - vectors2[innerCol])
				* mat[k + col];
		}
		col = i + j * nvectors;
		sum += innerSum * (vectors1[col] - vectors2[col]);
	}
	result[i] = sqrtf(sum);
}

extern "C"
void gpuMahaDouble(const int * vectorSetSize, const int * vectorLength,
	const double * vectorSet1, const double * vectorSet2, const double * mat,
	double * result)
{
	int
		nv = *vectorSetSize, n = *vectorLength;

	double
		* dVectorSet1, * dVectorSet2, * dMat, * dResult;

	size_t
		resultBytes = nv * sizeof(double),
		matBytes = n * n * sizeof(double),
		vectorSetBytes = nv * n * sizeof(double);

	cudaMalloc((void **) & dVectorSet1, vectorSetBytes);
	cudaMalloc((void **) & dVectorSet2, vectorSetBytes);
	cudaMalloc((void **) & dMat, matBytes);
	cudaMalloc((void **) & dResult, resultBytes);

	cudaMemcpy(dVectorSet1, vectorSet1, vectorSetBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(dVectorSet2, vectorSet2, vectorSetBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(dMat, mat, matBytes, cudaMemcpyHostToDevice);

	dim3 dimBlock(NTHREADS, 1, 1);
	int gx = ceil((double) nv /(double) dimBlock.x);
	dim3 dimGrid(gx, 1, 1);

	cudaDeviceSynchronize();
	dMahaDouble<<<dimBlock, dimGrid>>>(nv, n, dVectorSet1, dVectorSet2, dMat,
		dResult);
	cudaDeviceSynchronize();

	cudaMemcpy(result, dResult, resultBytes, cudaMemcpyDeviceToHost);

	cudaFree(dVectorSet1);
	cudaFree(dVectorSet2);
	cudaFree(dMat);
	cudaFree(dResult);
}
