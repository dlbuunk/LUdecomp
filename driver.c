#include <stdio.h>

extern void LUdecomp2(float **, int *, int);

float M[4][4]; /* alignment! */
float * pM[3];
int perm[3];

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

	LUdecomp2(pM, perm, 3);

	for (i=0; i<3; i++) printf("%f\t%f\t%f\n", pM[i][0], pM[i][1], pM[i][2]);

	return(0);
}
