#include <stdio.h>

#define DIM 4

extern float ** LUdecomp2(float *, float *, int *, int);

float M[DIM][DIM];
float T[DIM*(DIM+4)];
int perm[DIM];

int main()
{
	int i;
	float **pM;

	M[0][0] = 2.0;
	M[0][1] = 0.0;
	M[0][2] = 1.0;
	M[1][0] = 5.0;
	M[1][1] = 5.0;
	M[1][2] = 2.0;
	M[2][0] = 6.0;
	M[2][1] = 3.0;
	M[2][2] = 9.0;

	pM = LUdecomp2(&M[0][0], &T[0], &perm[0], 3);

	for (i=0; i<3; i++) printf("%f\t%f\t%f\n", pM[i][0], pM[i][1], pM[i][2]);

	return(0);
}
