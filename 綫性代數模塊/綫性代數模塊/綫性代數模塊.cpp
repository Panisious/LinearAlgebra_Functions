#include "LinearAlgebra.h"
#include "omp.h"
//	Notice: This part is used for debug	************************************

void randomvector(vector, int);
void randommatrix(matrix, int);
time_t randtime;
int c1, c2, c3;
//	************************************************************************


//Main function is used for test, the really important part is LinearAlgebra.h
int main()
{
	for (int X = 0; X < 10; X++)
	{
		c1 = clock();

		matrix m1 = Matrix_Create(4, "Matrix1");		

		randommatrix(m1, 100);

		double det=Matrix_Get_Determinant(m1);

		printf("%f\n", det);

		Matrix_Free(m1);

		c2 = clock() - c1;
		c3 += c2;
	}
	printf("\n\t\t\t\t---------------------------------------------------\n\n\t\t\t\t\t\tTotal Time: %.3fs\n\n\t\t\t\t---------------------------------------------------", c3 / 1000.0);
}


//	Notice: This part is used for debug	************************************

void randomvector(vector v, int threshold)
{
	time(&randtime);
	for (int i = randtime % 1000; i < 1000; i++)
		randtime = rand();
	for (int i = v->dim - 1; i > -1; i--)
		v->e[i] = rand() % threshold;
}
void randommatrix(matrix m,int threshold)
{
	time(&randtime);
	for (int i = randtime % 10000; i < 10000; i++)
		randtime = rand();
	int dim = m->dim;
	for(int row=0;row<dim;row++)
		for (int col = 0; col < dim; col++)
			m->e[row][col] = rand() % threshold;
}
//	************************************************************************