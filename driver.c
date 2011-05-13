#include <stdio.h>

#define DIM 4

extern void LUdecomp2(float *, float *, int *, int);

float M[DIM][DIM];

float * pM[DIM];

float W[DIM*(4+DIM)];

int perm[DIM];

int main()
{
	int i;

	M[0][0] = 2.0;
	M[0][1] = 0.0;
	M[0][2] = 1.0;
	M[1][0] = 5.0;
	M[1][1] = 5.0;
	M[1][2] = 2.0;
	M[2][0] = 6.0;
	M[2][1] = 3.0;
	M[2][2] = 9.0;

	pM[0] = &M[0][0];
	pM[1] = &M[1][0];
	pM[2] = &M[2][0];

	puts("Old version:");

	LUdecomp(pM, perm, 3);

	for (i=0; i<3; i++) printf("%d\t%f\t%f\t%f\n", perm[i], M[i][0], M[i][1], M[i][2]);

	puts("New version:");

	LUdecomp(M, W, perm, 3);

//	for (i=0; i<3; i++) printf("%d\t%f\t%f\t%f\n", perm[i], M[i][0], M[i][1], M[i][2]);

	return(0);
}

int LUdecomp(float **A, int * perm, int N)
{
	int j,i,k;
	float bjjtest, bjjmax;
	int bjjrow;
	float * tempr;
	int tempi;

	/* init perm[] */
	for (j=0; j<N; j++)
		perm[j] = j;

	for (j=0; j<N; j++)
	{
		if (j != N-1)
		{
			bjjmax = 0.0;
			bjjrow = -1;

			for (i=j; i<N; i++) /* search max bjj */
			{
				bjjtest = A[i][j];
				for (k=0; k<j; k++)
					bjjtest -= A[i][k]*A[k][j];
				if (bjjtest < 0.0)
					bjjtest = -bjjtest;
				if (bjjtest > bjjmax)
				{
					bjjmax = bjjtest;
					bjjrow = i;
				}
			}

			/* exit if all bjj's are 0.0 */
			if (bjjrow == -1)
				return(1);

			/* swap pointers */
			tempr = A[bjjrow];
			A[bjjrow] = A[j];
			A[j] = tempr;

			/* swap perm[] entries */
			tempi = perm[bjjrow];
			perm[bjjrow] = perm[j];
			perm[j] = tempi;
		}
		/* actual decomposition */
		for (i=0; i<=j; i++)
			for (k=0; k<i; k++)
				A[i][j] -= A[i][k]*A[k][j];
		for ( ; i<N; i++)
		{
			for (k=0; k<j; k++)
				A[i][j] -= A[i][k]*A[k][j];
			A[i][j] /= A[j][j];

		}
	}
	return(0);
}
