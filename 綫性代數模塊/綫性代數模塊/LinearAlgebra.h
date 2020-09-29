#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
//To judge if a pointer is NULL
#define IFNULL(a)	 if((a)==NULL)							\
					{										\
						printf("\n\nOut of memory!!!\n\n");	\
						exit(-1);							\
					}
//Storage the dimention of matrix
int dim;
//To calculate the determinant 
int counter_of_permutation = 1;
/*
						********************	Notice: This is the part of vector definition and operations	************************
*/
typedef struct tagVector {
	int dim;	//dimension
	char* name;
	double* e;	//elements
}*vector;
vector Vector_Create(int dimension, const char* name)
{
	vector v = (vector)malloc(sizeof(tagVector));
	IFNULL(v)
		v->dim = dimension;
	v->name = (char*)malloc(strlen(name) + 1);
	IFNULL(v->name)
		strcpy(v->name, name);
	v->e = (double*)malloc(sizeof(double) * dimension);
	IFNULL(v->e)
		printf("Vector %s created\n", v->name);
	for (int i = 0; i < dimension; i++)
		v->e[i] = 0;
	return v;
}
void Vector_Input(vector v, double* d)
{
	for (int i = 0; i < v->dim; i++)
		v->e[i] = d[i];
}
void Vector_Print(vector v)
{
	printf("Vector %s:[", v->name);
	for (int i = 0; i < v->dim - 1; i++)
		printf("%.3f, ", v->e[i]);
	printf("%.3f]T\n", v->e[v->dim - 1]);
}
vector Vector_Add_New(vector a, vector b)
{
	vector c = Vector_Create(a->dim, "Sum");
	for (int i = 0; i < c->dim; i++)
		c->e[i] += a->e[i] + b->e[i];
	return c;
}
vector Vector_Sub_New(vector a, vector b)
{
	vector c = Vector_Create(a->dim, "Diffrence");
	for (int i = 0; i < c->dim; i++)
		c->e[i] = a->e[i] - b->e[i];
	return c;
}
vector Vector_Multiply_by_Number_New(vector a, double d)
{
	vector c = Vector_Create(a->dim, "Product");
	for (int i = 0; i < c->dim; i++)
		c->e[i] += a->e[i] * d;
	return c;
}
inline void Vector_Add(vector a, vector b)
{
	for (int i = a->dim - 1; i >= 0; i--)
		a->e[i] += b->e[i];
}
inline void Vector_Sub(vector a, vector b)
{
	for (int i = a->dim - 1; i >= 0; i--)
		a->e[i] -= b->e[i];
}
inline void Vector_Multiply_by_Number(vector a, double k)
{
	for (int i = a->dim - 1; i >= 0; i--)
		a->e[i] *= k;
}
double Vector_Inner_Product(vector a, vector b)
{
	if (a->dim != b->dim)exit(-1);
	double sum = 0;
	for (int i = 0; i < a->dim; i++)
		sum += a->e[i] * b->e[i];
	return sum;
}
/*
						********************	Notice: This is the part of matrix basic definition and operations	************************
*/
typedef struct tagMatrix {
	int dim;	//dimension
	char* name;
	double** e;	//elements
	double det;
}*matrix;
matrix Matrix_Create(int dimension, const char* name)
{
	matrix m = (matrix)malloc(sizeof(tagMatrix));
	IFNULL(m)
		m->dim = dimension;
	if (name)
	{
		m->name = (char*)malloc(strlen(name) + 1);
		IFNULL(m->name)
			strcpy(m->name, name);
	}
	m->e = (double**)malloc(sizeof(double*) * dimension);
	IFNULL(m->e)
		for (int i = 0; i < dimension; i++)
		{
			m->e[i] = (double*)malloc(sizeof(double) * dimension);
			IFNULL(m->e[i])
				for (int j = 0; j < dimension; j++)
					m->e[i][j] = 0;
		}
	m->det = 0;
	//printf("Matrix %s created\n", m->name);
	return m;
}
void Matrix_Free(matrix m)
{
	//Safety
	if (!m)return;

	free(m->name);
	for (int i = 0; i < dim; i++)
		free(m->e[i]);
	free(m->e);
	free(m);
}
//Type 1: Print as Integer
//Type 2: Print as %.3f
//Type 3: Print as %.6f
void Matrix_Print(matrix m, int type)
{
	//Safety
	if (!m)
	{
		printf("\n\nMatrix not exist\n\n");
		return;
	}

	dim = m->dim;
	printf("\nMatrix %s:\n", m->name);
	if (type == 1)
		for (int row = 0; row < dim; row++)
		{
			printf("[ ");
			for (int col = 0; col < dim; col++)
			{
				if (m->e[row][col]<0 && m->e[row][col]>-1e-6)
					printf("%.0f ", -m->e[row][col]);
				else
					printf("%.0f ", m->e[row][col]);
			}
			printf("]\n");
		}
	else if (type == 2)
		for (int row = 0; row < dim; row++)
		{
			printf("[ ");
			for (int col = 0; col < dim; col++)
			{
				if (m->e[row][col]<0 && m->e[row][col]>-1e-6)
					printf("%.3f ", -m->e[row][col]);
				else
					printf("%.3f ", m->e[row][col]);
			}
			printf("]\n");
		}
	else if (type == 3)
		for (int row = 0; row < dim; row++)
		{
			printf("[ ");
			for (int col = 0; col < dim; col++)
			{
				if (m->e[row][col]<0 && m->e[row][col]>-1e-6)
					printf("%.6f ", -m->e[row][col]);
				else
					printf("%.6f ", m->e[row][col]);
			}
			printf("]\n");
		}
	printf("\n");
}
matrix Matrix_Add_New(matrix m, matrix b)
{
	//Safety
	if (!m||!b)
	{
		printf("\n\nMatrix not exist\n\n");
		return 0;
	}

	dim = m->dim;
	matrix c = Matrix_Create(m->dim, "Sum");
	for (int row = 0; row < dim; row++)
		for (int col = 0; col < dim; col++)
		{
			c->e[row][col] = m->e[row][col] + b->e[row][col];
		}
	return c;
}
matrix Matrix_Sub_New(matrix m, matrix b)
{
	//Safety
	if (!m || !b)
	{
		printf("\n\nMatrix not exist\n\n");
		return 0;
	}

	dim = m->dim;
	matrix c = Matrix_Create(m->dim, "Diffrence");
	for (int row = 0; row < dim; row++)
		for (int col = 0; col < dim; col++)
		{
			c->e[row][col] = m->e[row][col] - b->e[row][col];
		}
	return c;
}
//When k = 1, this function can be seen as Matrix_Copy
matrix Matrix_Multiply_by_Number_New(matrix m, double k)
{
	//Safety
	if (!m)
	{
		printf("\n\nMatrix not exist\n\n");
		return 0;
	}

	dim = m->dim;
	matrix c = Matrix_Create(dim, "Number_Product");
	for (int row = 0; row < dim; row++)
		for (int col = 0; col < dim; col++)
		{
			c->e[row][col] = m->e[row][col] * k;
		}
	return c;
}
matrix Matrix_Multiply_New(matrix m, matrix b)
{
	//Safety
	if (!m || !b)
	{
		printf("\n\nMatrix not exist\n\n");
		return 0;
	}

	dim = m->dim;
	double temp_sum;
	matrix c = Matrix_Create(dim, "Matrix_Product");
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			temp_sum = 0;
			for (int k = 0; k < dim; k++)
			{
				temp_sum += m->e[i][k] * b->e[k][j];
			}
			c->e[i][j] = temp_sum;
		}
	}
	return c;
}
void Matrix_Add(matrix m, matrix b)
{
	//Safety
	if (!m || !b)
	{
		printf("\n\nMatrix not exist\n\n");
		return;
	}

	dim = m->dim;
	for (int row = 0; row < dim; row++)
		for (int col = 0; col < dim; col++)
		{
			m->e[row][col] += b->e[row][col];
		}
}
void Matrix_Sub(matrix m, matrix b)
{
	//Safety
	if (!m || !b)
	{
		printf("\n\nMatrix not exist\n\n");
		return;
	}

	dim = m->dim;
	for (int row = 0; row < dim; row++)
		for (int col = 0; col < dim; col++)
		{
			m->e[row][col] -= b->e[row][col];
		}
}
void Matrix_Multiply_by_Number(matrix m, double k)
{
	//Safety
	if (!m)
	{
		printf("\n\nMatrix not exist\n\n");
		return;
	}

	dim = m->dim;
	for (int row = 0; row < dim; row++)
		for (int col = 0; col < dim; col++)
		{
			m->e[row][col] *= k;
		}
}
//Type 1: Row Permutation
//Type 2: Column Permutation
void Matrix_Permutation(matrix m, int position_1, int position_2, int type)
{
	//Safety
	if (!m)
	{
		printf("\n\nMatrix not exist\n\n");
		return;
	}

	if (position_1 == position_2)return;
	dim = m->dim;
	double temp;
	if (type == 1)
	{
		for (int i = 0; i < dim; i++)
		{
			temp = m->e[position_1][i];
			m->e[position_1][i] = m->e[position_2][i];
			m->e[position_2][i] = temp;
		}
	}
	else
	{
		for (int i = 0; i < dim; i++)
		{
			temp = m->e[i][position_1];
			m->e[i][position_1] = m->e[i][position_2];
			m->e[i][position_2] = temp;
		}
	}
	counter_of_permutation *= -1;
}
//Type 1: Row Permutation
//Type 2: Column Permutation
inline void Matrix_Elimination(matrix m, int source, int target, int start, int type, double rate)
{
	//Safety
	if (!m)
	{
		printf("\n\nMatrix not exist\n\n");
		return;
	}

	if (source == target)return;
	dim = m->dim;
	if (type == 1)
	{
		for (int i = start; i < dim; i++)
			m->e[target][i] -= m->e[source][i] * rate;
	}
	else
	{
		for (int i = start; i < dim; i++)
			m->e[i][target] -= m->e[i][source] * rate;
	}
	//Matrix_Print(m, 2);
}
matrix Matrix_Transform_to_Upper(matrix m)
{
	//Safety
	if (!m)
	{
		printf("\n\nMatrix not exist\n\n");
		return 0;
	}

	double rate;
	dim = m->dim;
	counter_of_permutation = 1;
	matrix I = Matrix_Create(dim, "Inverse");
	for (int i = 0; i < dim; i++)
		I->e[i][i] = 1;
	for (int col = 0; col < dim; col++)
	{
		int elem = 0, position = 0;
		for (int row = 0; row < dim; row++)
		{
			if (fabs(m->e[row][col]) > 1e-6)
			{
				elem++;
				position = row;
			}
		}
		if (elem == 0)
		{
			printf("\n\nThis is SINGULAR matrix\n\n");
			return 0;
		}
		else
		{
			if (fabs(m->e[col][col]) < 1e-6 && position != col)
			{
				Matrix_Permutation(m, col, position, 1);
				Matrix_Permutation(I, col, position, 1);
				counter_of_permutation *= -1;
			}
			for (int i = col; i < dim; i++)
			{
				rate = m->e[i][col] / m->e[col][col];
				Matrix_Elimination(m, col, i, col, 1, rate);
				Matrix_Elimination(I, col, i, 0, 1, rate);
			}
		}
	}
	double det = 1;
	for (int i = 0; i < dim; i++)
		det *= m->e[i][i];
	if (fabs(det) < 1e-6)
	{
		printf("\n\nThis is SINGULAR matrix\n\n");
		Matrix_Free(I);
		return 0;
	}
	m->det = det;
	return I;
}
void Matrix_Transform_Upper_to_Diagonal(matrix m, matrix I)
{
	//Safety
	if (!m || !I)
	{
		printf("\n\nMatrix not exist\n\n");
		return;
	}

	double rate;
	dim = m->dim;
	for (int col = dim - 1; col > 0; col--)
	{
		for (int row = col - 1; row >= 0; row--)
		{
			if (fabs(m->e[row][col]) > 1e-6)
			{
				rate = m->e[row][col] / m->e[col][col];
				Matrix_Elimination(m, col, row, col, 1, rate);
				Matrix_Elimination(I, col, row, 0, 1, rate);
			}
		}
	}
	for (int row = 0; row < dim; row++)
	{
		for (int col = 0; col < dim; col++)
		{
			I->e[row][col] /= m->e[row][row];
		}
		m->e[row][row] /= m->e[row][row];
	}
}
int Matrix_Check_if_I(matrix m)
{
	//Safety
	if (!m)
	{
		printf("\n\nMatrix not exist\n\n");
		return -1;
	}

	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			if (i == j)
			{
				if (fabs(m->e[i][i] - 1) > 1e-6)
					return 0;
			}
			else
				if (fabs(m->e[i][j]) > 1e-6)
					return 0;
	return 1;
}
/*
						********************	Notice: This is the part of Useful vector and matrix operations	************************
*/
matrix Matrix_Get_Inverse(matrix m)
{
	//Safety
	if (!m)
	{
		printf("\n\nMatrix not exist\n\n");
		return 0;
	}

	matrix b = Matrix_Multiply_by_Number_New(m, 1);
	matrix Inverse = Matrix_Transform_to_Upper(b);
	if (!Inverse)
		return 0;
	Matrix_Transform_Upper_to_Diagonal(b, Inverse);
	Matrix_Free(b);
	return Inverse;
}
double Matrix_Get_Determinant(matrix m)
{
	//Safety
	if (!m)
	{
		printf("\n\nMatrix not exist\n\n");
		return 0;
	}

	if (!m->det)
	{
		matrix temp = Matrix_Multiply_by_Number_New(m, 1);
		Matrix_Transform_to_Upper(temp);
		double det = temp->det;
		m->det = det;
		Matrix_Free(temp);
		return det;
	}
	return m->det;
}
matrix Matrix_Get_Power(matrix m, int power)
{
	//Safety
	if (!m)
	{
		printf("\n\nMatrix not exist\n\n");
		return 0;
	}

	matrix mPower = Matrix_Multiply_by_Number_New(m, 1);
	matrix temp, p_mPower;
	for (int i = 1; i < power; i++)
	{
		temp = Matrix_Multiply_New(mPower, m);
		p_mPower = mPower;
		mPower = temp;
		Matrix_Free(p_mPower);
	}
	return mPower;
}
matrix Linear_Combination_New(matrix m, vector v)
{
	//Safety
	{	
		if (!m)
		{
			printf("\n\nMatrix not exist\n\n");
			return 0;
		}
		if (!v)
		{
			printf("\n\nVector not exist\n\n");
			return 0;
		}
		int dim;
		if ((dim = m->dim) != v->dim)
		{
			printf("\n\nDimension not match:\nMatrix:%d\nVector:%d", m->dim, v->dim);
			return 0;
		}
	}

	matrix M = Matrix_Create(dim, "Combination");
	for (int row = 0; row < dim; row++)
	{
		for (int col = 0; col < dim; col++)
		{
			M->e[row][col] = m->e[row][col] * v->e[col];
		}
	}
	return M;
}
void Linear_Combination(matrix m, vector v)
{
	//Safety
	{	
		if (!m)
		{
			printf("\n\nMatrix not exist\n\n");
			return;
		}
		if (!v)
		{
			printf("\n\nVector not exist\n\n");
			return;
		}
		int dim;
		if ((dim = m->dim) != v->dim)
		{
			printf("\n\nDimension not match:\nMatrix:%d\nVector:%d", m->dim, v->dim);
			return;
		}
	}

	for (int row = 0; row < dim; row++)
	{
		for (int col = 0; col < dim; col++)
		{
			m->e[row][col] *= v->e[col];
		}
	}
}