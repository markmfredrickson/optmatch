#include<R.h>

// computes sqrt((u_i - v_i) * m * (u_i - v_i)) for two sets of vectors
// u and v and a matrix m. The matrix m should have dim n x n where each
// u_i and v_i are of length n. This is intended to fill in the last part of a
// mahabalonis distance calculation for a set of vectors u and v,
// one set a treatment and the other a control.

// The matrix m is expected to be stored column major in memory
// the ith column's elements are stored consecutively at m + i * n
// the jth row will be at m + 0 * n + j, m + n + j, m + 2 * n + j, ...

// each vector is expected to be a row vector of length n, with the whole
// set stored as an nv x n matrix w in row major order
// to get at the jth vector in the set, you must index the elements as
// w + j * n + 0, w + j * n + 1, w + j * n + 2, ..., w + j * n + (n - 1)

void mahalanobisHelper(const int * vectorSetSize, const int * vectorLength,
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
