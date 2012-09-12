#include <R.h>
#include "Rdefines.h"
#include "R_ext/Rdynload.h"
#ifdef WINDOWS
#  include <windows.h>
#endif

#ifndef max
 	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

double find_closest_pair(int n, double** distmatrix, int* ip, int* jp)
/*
This function searches the distance matrix to find the pair with the shortest
distance between them. The indices of the pair are returned in ip and jp; the
distance itself is returned by the function.

n          (input) int
The number of elements in the distance matrix.

distmatrix (input) double**
A ragged array containing the distance matrix. The number of columns in each
row is one less than the row index.

ip         (output) int*
A pointer to the integer that is to receive the first index of the pair with
the shortest distance.

jp         (output) int*
A pointer to the integer that is to receive the second index of the pair with
the shortest distance.
*/
{ int i, j;
  double distance = distmatrix[1][0];
  for (i = 0; i < n; i++)
  { for (j = 0; j < i; j++)
    { if (distmatrix[i][j]<distance)
      { distance = distmatrix[i][j];
        *ip = i;
        *jp = j;
      }
    }
  }
  return distance;
}





SEXP palcluster (SEXP rdistmatrix)
/*
Purpose
=======

The palcluster routine performs clustering using pairwise average
linking on the given distance matrix.

Arguments
=========

nelements     (input) int
The number of elements to be clustered.

distmatrix (input) double**
The distance matrix, with nelements rows, each row being filled up to the
diagonal. The elements on the diagonal are not used, as they are assumed to be
zero. The distance matrix will be modified by this routine.

Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, palcluster returns NULL.
========================================================================
*/
{ double *distvec;
  int nelements;
  if (isMatrix(rdistmatrix) && isReal(rdistmatrix)) {
      distvec = REAL(rdistmatrix);
      nelements = INTEGER(GET_DIM(rdistmatrix))[1];
  }
  else {
      //printf("invalid matrix pal.\n");
      return R_NilValue;
  }

  double sil[nelements-2];

  int j=0, i, n, w=0, sum, *clusterid, *map, *number;
  double within=0, between=0, **copy, **distmatrix;

  clusterid = malloc(nelements*sizeof(int));
  number = malloc(nelements*sizeof(int));

  distmatrix = malloc(nelements*sizeof(double*));
  copy = malloc(nelements*sizeof(double*));
  map = malloc(nelements*sizeof(int));
  if (copy && clusterid && number && map)
	  for (j = 0; j < nelements; j++) {
		number[j] = 1;							/* the number of elements in each cluster	*/
		clusterid[j] = j;						/* maps cluster -> node number				*/
		map[j] = j;								/* maps gene    -> cluster number			*/
		distmatrix[j] = malloc(j*sizeof(double));
		copy[j] = malloc(j*sizeof(double));		/* copy of distmatrix						*/
		if (!copy[j]) break;
		for (i=0; i<j; ++i) {
		  distmatrix[j][i] = distvec[i*nelements+j];
		  copy[j][i] = distmatrix[j][i];
		  between += copy[j][i];				/* initialize sum of between cluster dist	*/
		}
	  }
  if (j != nelements) { 
	for (i = 0; i < j; i++)
		 free(copy[i]);	
    free(clusterid);
    free(number);
    free(copy);
	free(map);
	return NULL;
  }

  for (n = nelements; n > 1; n--) { 
	/* Print Silhouette coefficient */
	if (w) {
		double betw = between/(n*(n-1)/2), with=within/w;
		sil[nelements-n-1] = (betw-with)/max(betw, with);
	}
    int is = 1, js = 0;
    double notused = find_closest_pair(n, distmatrix, &is, &js);


    /* Fix the distances and update between cluster distances */
    sum = number[is] + number[js];
	between -= distmatrix[is][js];
    for (j = 0; j < js; j++) {
	  between -= distmatrix[is][j]+distmatrix[js][j];
      distmatrix[js][j] = (distmatrix[is][j]*number[is] + distmatrix[js][j]*number[js]) / sum;
	  between += distmatrix[js][j];
    }
    for (j = js+1; j < is; j++) {
	  between -= distmatrix[is][j]+distmatrix[j][js];
      distmatrix[j][js] = (distmatrix[is][j]*number[is] + distmatrix[j][js]*number[js]) / sum;
	  between += distmatrix[j][js];
   }
    for (j = is+1; j < n; j++) {
	  between -= distmatrix[j][is]+distmatrix[j][js];
      distmatrix[j][js] = (distmatrix[j][is]*number[is] + distmatrix[j][js]*number[js]) / sum;
	  between += distmatrix[j][js];
   }

    for (j = 0; j < is; j++) distmatrix[is][j] = distmatrix[n-1][j];
    for (j = is+1; j < n-1; j++) distmatrix[j][is] = distmatrix[n-1][j];

    /* Update number of elements in the clusters */
    number[js] = sum;
    number[is] = number[n-1];

    /* Update clusterids */
    clusterid[js] = n-nelements-1;
    clusterid[is] = clusterid[n-1];
	
	/* Update within cluster distances */
	for (j=0; j<nelements; ++j) {
	  if (map[j]==js) 
		for (i=j+1; i<nelements; ++i) 
		  if (map[i] == is) {
			within += copy[i][j];
			++w;
		  }
	  if (map[j]==is) {
		for (i=j+1; i<nelements; ++i)
		  if (map[i] == js) {
			within += copy[i][j];
			++w;
		  }
		map[j] = js;
	  }
	  if (map[j]==n-1) 
		map[j] = is;
	}	
  }
  for (j = 0; j < nelements; j++)
	 free (copy[j]);
  free(clusterid);
  free(number);
  free(copy);
  free(map);

  
  SEXP retval;
  PROTECT(retval = NEW_NUMERIC(nelements-2));
  double *rsil = NUMERIC_POINTER(retval);
  for(int i=0; i<nelements-2; i++){rsil[i] = sil[i];}
//  DOUBLE_DATA(retval)[0] = 55;
  UNPROTECT(1);

  return retval;
}


R_CallMethodDef callMethods[] =
{
    {"palcluster", (DL_FUNC)&palcluster, 1},
    {NULL,NULL, 0}
};

void R_init_cluster(DllInfo *dll)
{
    R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
}





