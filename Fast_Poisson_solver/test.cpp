#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int Print_Matrix(double **A, int N);
int Print_Vector(double *x, int N);
int DST(double *x, int N);
int Print_Tensor(double ***A, int N);

int main()
{
	int i, j, k, L;
	double ***F, *f, *temp;
	k = 0;
	L = 3;
	temp = (double *) malloc(L*sizeof(double));
	F = (double ***) malloc(L*sizeof(double**));
	f = (double *) malloc(L*L*L*sizeof(double));
	for(i=0;i<L;i++)
	{
		F[i] = (double **) malloc(L*sizeof(double*));
		for(j=0;j<L;j++,k=k+L)
		{
			F[i][j] = f + k;
		}
	}

	for(i=0;i<L;i++)
	{
		for(j=0;j<L;j++)
		{
			for(k=0;k<L;k++)
			{
				F[i][j][k] = i + j + k;
			}
		}
	}
	
	Print_Tensor(F, L);
	for(k=0;k<L;k++) 
	{
		for(j=0;j<L;j++)
		{
			for(i=0;i<L;i++)
			{
				temp[i] = F[i][j][k];
			}
			DST(temp, L);
			for(i=0;i<L;i++)
			{
				F[i][j][k] = temp[i];
			}
		}
	}

	Print_Tensor(F, L);

	
	return 0;
	
}

int Print_Tensor(double ***A, int N)
{
	int i, j, k;
	for(k=0;k<N;k++)
	{
		printf("k = %d \n", k);
		for(j=0;j<N;j++)
		{
			for(i=0;i<N;i++)
			{
				printf("F[%d][%d][%d] = %f ", i, j, k, A[i][j][k]);
			}
			printf("\n");
		}
	}
}

int Print_Matrix(double **A, int N)
{
	int i, j;
	for(i=0;i<N;++i)
	{
		for(j=0;j<N;++j)
		{
			printf("%f ",A[i][j]);
		}
		printf("\n");
	}
}
int Print_Vector(double *x, int N)
{
	int i;
	for(i=0;i<N;i++) 
	{
		printf(" x[%d] = %f \n", i, x[i]);
	}
}
// Fast Fourier Transform in place for N = 2^p 
int DST(double *x, int N)
{
	int i, j, k, n, M, K;
	double t_r, t_i, *x_r, *x_i, *y_r, *y_i;
	
	K = 2*N + 2;	
	x_r = (double *) malloc(N*sizeof(double));
	x_i = (double *) malloc(N*sizeof(double));
	y_r = (double *) malloc(K*sizeof(double));
	y_i = (double *) malloc(K*sizeof(double));
	
	for(i=0;i<N;i++)
	{
		x_r[i] = x[i];
		x_i[i] = 0.0;
	}

	// expand y[n] to 2N+2-points from x[n]
	y_r[0] = y_i[0] = 0.0;
	y_r[N+1] = y_i[N+1] = 0.0;
	for(i=0;i<N;i++)
	{
		y_r[i+1] = x_r[i];
		y_i[i+1] = x_i[i];
		y_r[N+i+2] = -1.0*x_r[N-1-i];
		y_i[N+i+2] = -1.0*x_i[N-1-i];
	}
	
	
	i = j = 0;
	while(i < K)
	{
		if(i < j)
		{
			// swap y[i], y[j]
			t_r = y_r[i];
			t_i = y_i[i];
			y_r[i] = y_r[j];
			y_i[i] = y_i[j];
			y_r[j] = t_r;
			y_i[j] = t_i;
		}
		M = K/2;
		while(j >= M & M > 0)
		{
			j = j - M;
			M = M / 2;
		}
		j = j + M;		
		i = i + 1;
	}
	// Butterfly structure
	double theta, w_r, w_i;
	n = 2;
	while(n <= K)
	{
		for(k=0;k<n/2;k++)
		{
			theta = -2.0*k*M_PI/n;
			w_r = cos(theta);
			w_i = sin(theta);
			for(i=k;i<K;i+=n)
			{
				j = i + n/2;
				t_r = w_r * y_r[j] - w_i * y_i[j];
				t_i = w_r * y_i[j] + w_i * y_r[j];
				

				y_r[j] = y_r[i] - t_r;
				y_i[j] = y_i[i] - t_i;
				y_r[i] = y_r[i] + t_r;
				y_i[i] = y_i[i] + t_i;

			}
		}
		n = n * 2;
	}
	
	// After fft(y[k]), Y[k] = fft(y[k]), Sx[k] = i*Y[k+1]/2
	for(k=0;k<N;k++)
	{
		x[k] = -1.0*y_i[k+1]/2;
	}
	
}
