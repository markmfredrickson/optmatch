#include<R.h>
// computes the dot product of two vectors u and v of length n

double dotProduct(int n, const double * u, const double * v) {
	double sum = 0.0;
	for(int i = 0; i < n; i++)
		sum += u[i] * v[i];
	return sum;
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

void mahalanobisHelper(const int * vectorSetSize, const int * vectorLength,
	const double * vectorSet1, const double * vectorSet2, const double * mat,
	double * result)
{
	int
		nv = *vectorSetSize, n = *vectorLength;

	double * temp1 = Calloc(n, double);
	double * temp2 = Calloc(n, double);

	for(int i = 0; i < nv; i++) {
		for(int j = 0; j < n; j++)
			temp1[j] = vectorSet1[i + j * nv] - vectorSet2[i + j * nv];
		for(int j = 0; j < n; j++)
			temp2[j] = dotProduct(n, temp1, mat + j * n);
		result[i] = sqrt( dotProduct(n, temp1, temp2) );
	}
	Free(temp1);
	Free(temp2);
}

void unrolledMahalanobisHelper(
	const int * vectorSetSize, const int * vectorLength,
	const double * vectorSet1, const double * vectorSet2, const double * mat,
	double * result)
{
	int
		j, k,
		nv = *vectorSetSize, n = *vectorLength;
	double
		sum, innerSum;
	const double
		* v1i, * v2i, * matCol;
	
	for(int i = 0; i < nv; i++) {
		sum = 0;
		v1i = vectorSet1 + i * n;
		v2i = vectorSet2 + i * n;
		for(j = 0; j < n; j++) {
			innerSum = 0;
			matCol = mat + j * n;
			for(k = 0; k < n; k++)
				innerSum += (v1i[k] - v2i[k]) * matCol[k];
			sum += innerSum * (v1i[j] - v2i[j]);
		}
		result[i] = sqrt(sum);
	}
}
