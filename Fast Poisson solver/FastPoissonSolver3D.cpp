#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int Print_Tensor(double ***A, int N);
int Exact_Solution(double ***U, int N);
int Exact_Source(double ***F, int N);
double Error(double ***X, double ***U, int N);
int DST(double *x, int N);
int Fast_Poisson_Solver3D(double ***F, double ***X, int N);

int main()
{
	//Solve the Poisson equation in 3D with Dirichlet boundary condition
	//Take the real solution U[x][y][z] = sin(pi*x)*sin(2*pi*y)*sin(pi*z)
	//Then laplacian U[x][y][z] = F[x][y][z] = -6*pi*pi*sin(pi*x)*sin(2*pi*y)*sin(pi*z)
	int i, j, k, N, L, p; 
	double *x, ***X, *u, ***U, *f, ***F;
	clock_t t1, t2;
	p = pow(2, 8);
	// Create memory for solving Ax = b
	for(N=4;N<=p;N*=2)
	{
		L = (N-1);
		// L*L*L is the total number of unknowns		
		// Prepare for three dimensional unknown F
		k=0;
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
		k=0;
		X = (double ***) malloc(L*sizeof(double**));
		x = (double *) malloc(L*L*L*sizeof(double));
		for(i=0;i<L;i++)
		{
			X[i] = (double **) malloc(L*sizeof(double*));
			for(j=0;j<L;j++,k=k+L)
			{
				X[i][j] = x + k;
			}
		}
				
		// Prepare for three dimensional unknowns U 
		k=0;
		U = (double ***) malloc(L*sizeof(double**));
		u = (double *) malloc(L*L*L*sizeof(double));
		for(i=0;i<L;i++)
		{
			U[i] = (double **) malloc(L*sizeof(double*));
			for(j=0;j<L;j++,k=k+L)
			{
				U[i][j] = u + k;
			}
		}
		
		Exact_Solution(U, L);
		Exact_Source(F, L);
		t1 = clock();
		Fast_Poisson_Solver3D(F, X, L);
		t2 = clock();
		//Print_Tensor(X,L);
		//Print_Tensor(U,L);
		printf("Fast Poisson Solver: %f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
		printf("For N = %d, error = %e \n", N, Error(X, U, L));

		free(X);
		free(U);
		free(F);
		//system("pause");
	}
	return 0;
}
int Exact_Solution(double ***U, int N)
{
	// Put the exact solution 
	int i,j,k;
	double x, y, z, h;
	h = 1.0/(N+1);
	for(i=0;i<N;++i)
	{
		for(j=0;j<N;++j)
		{
			for(k=0;k<N;k++)
			{
				x = (i+1)*h;
				y = (j+1)*h;
				z = (k+1)*h;
				U[i][j][k] = sin(M_PI*x)*sin(2*M_PI*y)*sin(M_PI*z);
			}

		}
	}
}
int Exact_Source(double ***F, int N)
{
	int i,j,k;
	double x, y, z, h;
	h = 1.0/(N+1);
	for(i=0;i<N;++i)
	{
		for(j=0;j<N;++j)
		{
			for(k=0;k<N;k++)
			{
				x = (i+1)*h;
				y = (j+1)*h;
				z = (k+1)*h;
				F[i][j][k] = -(1.0+4.0+1.0)*h*h*M_PI*M_PI*sin(M_PI*x)*sin(2*M_PI*y)*sin(M_PI*z);	
			}
		}
	}	
}

int Print_Tensor(double ***A, int N)
{
	int i, j, k;
	for(k=0;k<N;++k)
	{
		printf("\n k = %d \n \n", k);
		for(i=0;i<N;++i)
		{
			for(j=0;j<N;j++)
			{
				printf("%f ",A[i][j][k]);
			}
			printf("\n");
		}
	}
}

double Error(double ***X, double ***U, int N)
{
	// return max_i |X[i][j][k] - U[i][j][k]|
	int i, j, k;
	double v = 0.0, e;
	
	for(i=0;i<N;++i)
	{
		for(j=0;j<N;j++)
		{
			for(k=0;k<N;k++)
			{
				e = fabs(X[i][j][k] - U[i][j][k]);
				if(e>v) v = e;
			}
		}
	}
	return v;
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
	free(x_r);
	free(x_i);
	free(y_r);
	free(y_i);
}

int Fast_Poisson_Solver3D(double ***F, double ***X, int N)
{
	int i, j, k;
	double h, s, S, *lamda, *temp;
	
	//The iDST's scalar
	s = 2.0/(N+1);
	S = s*s*s;
	//Set tmep to save vector
	temp = (double *) malloc(N*sizeof(double));
	//Lamda for eigenvalues	
	lamda = (double *) malloc(N*sizeof(double));
	h = 1.0/(N+1);
	for(i=0;i<N;i++)
	{
		lamda[i] = 2 - 2*cos((i+1)*M_PI*h);
	}
	//Do the discrete sine transform about F
	for(k=0;k<N;k++) 
	{
		for(j=0;j<N;j++)
		{
			for(i=0;i<N;i++)
			{
				temp[i] = F[i][j][k];
			}
			DST(temp, N);
			for(i=0;i<N;i++)
			{
				F[i][j][k] = temp[i];
			}
		}
		
	}
	for(k=0;k<N;k++) 
	{
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				temp[j] = F[i][j][k];
			}
			DST(temp, N);
			for(j=0;j<N;j++)
			{
				F[i][j][k] = temp[j];
			}
		}
		
	}
	for(i=0;i<N;i++) 
	{
		for(j=0;j<N;j++)
		{
			for(k=0;k<N;k++)
			{
				temp[k] = F[i][j][k];
			}
			DST(temp, N);
			for(k=0;k<N;k++)
			{
				F[i][j][k] = temp[k];
			}
		}
		
	}

	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++) 
		{
			for(k=0;k<N;k++)
			{
				X[i][j][k] = -S*F[i][j][k]/(lamda[i] + lamda[j] + lamda[k]);				
			}

		}
	}

	for(k=0;k<N;k++) 
	{
		for(j=0;j<N;j++)
		{
			for(i=0;i<N;i++)
			{
				temp[i] = X[i][j][k];
			}
			DST(temp, N);
			for(i=0;i<N;i++)
			{
				X[i][j][k] = temp[i];
			}
		}
		
	}
	for(k=0;k<N;k++) 
	{
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				temp[j] = X[i][j][k];
			}
			DST(temp, N);
			for(j=0;j<N;j++)
			{
				X[i][j][k] = temp[j];
			}
		}
		
	}
	for(i=0;i<N;i++) 
	{
		for(j=0;j<N;j++)
		{
			for(k=0;k<N;k++)
			{
				temp[k] = X[i][j][k];
			}
			DST(temp, N);
			for(k=0;k<N;k++)
			{
				X[i][j][k] = temp[k];
			}
		}
		
	}
}




