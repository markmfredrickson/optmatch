#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

#define SMAHAL_CUTOFF 1e-10

void rank(int n, const double * data, double * ranks) {
  double
    * sorted = Calloc(n, double);
  int
    * order = Calloc(n, int);

  memcpy(sorted, data, n * sizeof(double));
  for(int i = 0; i < n; i++)
    order[i] = i;

  rsort_with_index(sorted, order, n);

  for(int i = 0; i < n; i++)
    ranks[order[i]] = i + 1.0;

  Free(order);
  Free(sorted);
}

Rboolean rerank_dups(int n, const double * data, double * ranks) {
    int
      j, ndups,
      * dup_idx = Calloc(n, int),
      * visited = Calloc(n, int);
    double mean_rank;
    Rboolean any_ties = FALSE;

    for(int i = 0; i < n ; i++) {
        if(visited[i] == 1) continue;
        ndups = 1;
        dup_idx[0] = i;
        mean_rank = ranks[i];
        for(j = i + 1; j < n; j++) {
            if(data[j] == data[i]) {
                dup_idx[ndups] = j;
                ndups++;
		any_ties = TRUE;
                mean_rank += ranks[j];
                visited[j] = 1;
            }
        }
        if(ndups == 1) continue;
        mean_rank /= (double) ndups;
        for(j = 0; j < ndups; j++)
            ranks[dup_idx[j]] = mean_rank;
    }
    Free(dup_idx);
    Free(visited);
    return any_ties;
}

double mean(const double * x, int n) {
    double sum = x[0];
    for(int i = 1; i < n; i++)
        sum += x[i];
    return sum / (double) n;
}

double cov(int n, const double * x, const double * y) {
    double
        sum = 0,
        mean_x = mean(x, n),
        mean_y = mean(y, n);

    for(int i = 0; i < n; i++)
        sum += (x[i] -  mean_x) * (y[i] - mean_y);
    return sum / (double) (n - 1);
}

double var(int n, const double * x) {
  return cov(n, x, x);    
}

/*  implements the following R code in C
    vuntied <- var(1:n)
    rat <- sqrt(vuntied/diag(cv))
    cv <- diag(rat) %*% cv %*% diag(rat)
*/
void adjust_ties(int nr, int nc, double * covs) {
  double * seq = Calloc(nr, double);
  for(int i = 0; i < nr; i++)
    seq[i] = i + 1.0;

  double vuntied = var(nr, seq);
  Free(seq);

  double * rat = Calloc(nc, double);
  for(int i = 0; i < nc; i++)
    rat[i] = sqrt(vuntied / covs[i * nc + i]);

  for(int i = 0; i < nc; i++) {
    for(int j = i; j < nc; j++) {
      double elt = rat[i] * covs[i * nc + j] * rat[j];
      covs[i * nc + j] = covs[j * nc + i] = elt;
    }
  }
  Free(rat);
}

void transpose_sq(int n, double * mat) {
  double elt;
  for(int i = 0; i < n; i++) {
    for(int j = i + 1; j < n; j++) {
      elt = mat[i * n + j];
      mat[i * n + j] = mat[j * n + i];
      mat[j * n + i] = elt;
    }
  }
}

void mult_sq_diag(int n, double * mat, const double * diag) {
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      mat[i * n + j] *= diag[i];
    }
  }
}

double dmax(int n, const double * data) {
  double max = data[0];
  for(int i = 1; i < n; i++) {
    if(data[i] > max)
      max = data[i];
  }
  return(max);
}

void ginv_square(double * square_mat, int n) {
    char jobz = 'A';
    int
      info,
      work_size = 4 * n * n + 7 * n,
      * iwork = Calloc(8 * n, int);
    double
      * s = Calloc(n, double),
      * u = Calloc(n * n, double),
      * vt = Calloc(n * n, double),
      * work = Calloc(work_size, double);

    F77_CALL(dgesdd)(&jobz, &n, &n, square_mat, &n, s, u, &n, vt, &n,
		     work, &work_size, iwork, &info);
    Free(work);
    Free(iwork);

    if(info != 0) {
      Free(u);
      Free(vt);
      Free(s);
    }

    if(info < 0)
      error("dgesdd: problem with one of the arguments");
    else if(info > 0)
      error("dgesdd: dbdsdc did not converge, updating process failed");

    double
      smax = dmax(n, s),
      cutoff = 1e-10 * smax;

    for(int i = 0; i < n; i++) {
      if(s[i] > cutoff)
	s[i] = 1.0 / s[i];
      else
	s[i] = 0.0;
    }

    transpose_sq(n, vt);
    mult_sq_diag(n, vt, s);

    char
      transa = 'N', transb = 'T';
    double
      alpha = 1.0, beta = 0.0;

    F77_CALL(dgemm)(&transa, &transb, &n, &n, &n,
		    &alpha, vt, &n, u, &n,
		    &beta, square_mat, &n);
    Free(u);
    Free(vt);
    Free(s);
}

void mahalanobis(int nr, int nc, const double * x, const double * center,
		 const double * icov, double * distances)
{
  double * recentered = Calloc(nr * nc, double);
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++)
      recentered[j * nr + i] = x[j * nr + i] - center[j];
  }

  char
    transa = 'N', transb = 'N';
  double
    alpha = 1.0, beta = 0.0,
    * mat_mult = Calloc(nr * nc, double);

  F77_CALL(dgemm)(&transa, &transb, &nr, &nc, &nc,
		  &alpha, recentered, &nr, icov, &nc,
		  &beta, mat_mult, &nr);

  for(int i = 0; i < nr * nc; i++)
    mat_mult[i] *= recentered[i];

  Free(recentered);

  // compute row sums of mat_mult for answer
  for(int i = 0; i < nr; i++) {
    double sum  = 0.0;
    for(int j = 0; j < nc; j++) {
      sum += mat_mult[j * nr + i];
    }
    distances[i] = sum;
  }
  Free(mat_mult);
}

typedef struct dmat {
  int nr, nc;
  double * data;
} DMAT;

DMAT * smahal_nosexp(int nr, int nc, double * data, int * z) {
    double
      * col_i,
      * ranks = Calloc(nr * nc, double),
      * ranks_i = Calloc(nr, double);

    Rboolean
      tie,
      any_ties = FALSE;

    memcpy(ranks, data, nr * nc * sizeof(double));
    for(int i = 0; i < nc; i++) {
      col_i = ranks + i * nr;
      rank(nr, col_i, ranks_i);
      tie = rerank_dups(nr, col_i, ranks_i);
      any_ties = tie || any_ties;
      memcpy(col_i, ranks_i, nr * sizeof(double));
    }
    Free(ranks_i);

    double * covs_inv = Calloc(nc * nc, double);
    for(int i = 0; i < nc; i++) {
      for(int j = i; j < nc; j++) {
	  covs_inv[i * nc + j] = covs_inv[j * nc + i] =
	    cov(nr, ranks + i * nr, ranks + j * nr);
      }
    }
    if(any_ties == TRUE) adjust_ties(nr, nc, covs_inv);    
    ginv_square(covs_inv, nc);

    int ncontrol, ntreat = 0;
    for(int i = 0; i < nr; i++) {
      if(z[i] == TRUE)
	ntreat++;
    }
    ncontrol = nr - ntreat;

    double
      * ranks_control = Calloc(ncontrol * nc, double),
      * ranks_treat = Calloc(ntreat * nc, double);

    int
      c_treat_row = 0, c_control_row = 0;
    for(int i = 0; i < nr; i++) {
      if(z[i] == TRUE) {
	for(int j = 0; j < nc; j++)
	  ranks_treat[c_treat_row + j * ntreat] = ranks[i + j * nr];
	c_treat_row++;
      } else {
	for(int j = 0; j < nc; j++)
	  ranks_control[c_control_row + j * ncontrol] = ranks[i + j * nr];
	c_control_row++;
      }
    }

    DMAT * out_distances = Calloc(1, DMAT);
    if(out_distances == NULL)
      error("smahal_nosexp:391:NULL Calloc\n");

    out_distances->data = Calloc(ntreat * ncontrol, double);
    out_distances->nr = ntreat;
    out_distances->nc = ncontrol;

    double * out_row_i = Calloc(ncontrol, double);
    double * treat_row_i = Calloc(nc, double);

    for(int i = 0; i < ntreat; i++) {
      for(int j = 0; j < nc; j++)
	treat_row_i[j] = ranks_treat[i + j * ntreat];

      mahalanobis(ncontrol, nc, ranks_control, treat_row_i, covs_inv,
		  out_row_i);
      for(int j = 0; j < ncontrol; j++)
	out_distances->data[i + j * ntreat] = out_row_i[j];
    }

    Free(treat_row_i);
    Free(ranks_control);
    Free(ranks_treat);
    Free(out_row_i);
    Free(covs_inv);
    Free(ranks);

    return out_distances;
}

SEXP smahal(SEXP index, SEXP data, SEXP z) {
  DMAT * ans = smahal_nosexp(nrows(data), ncols(data), REAL(data), LOGICAL(z));
  if(ans == NULL || ans->nr < 1 || ans->nc <1)
    error("smahal_nosexp returned an invalid answer");

  SEXP out;
  PROTECT(out = allocMatrix(REALSXP, ans->nr, ans->nc));
  memcpy(REAL(out), ans->data, ans->nr * ans->nc * sizeof(double));

  Free(ans->data);
  Free(ans);
  UNPROTECT(1);

  return out;
}
