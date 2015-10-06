#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int FFTr3(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int Initial(double *x, double *y, int N);
//int Print_Complex_Vector(double *y_r, double *y_i, int N);
int Generate_N(int p, int q, int r);

int main()
{
	int k, n, p, q, r, N;
	double *y_r, *y_i, *x_r, *x_i, w_r, w_i;
	clock_t t1, t2;
	
	printf("Please input p q r=");
	scanf("%d %d %d", &p, &q, &r);
	N = Generate_N(p, q, r);
	printf("N=2^%d 3^%d 5^%d = %d\n",p,q,r,N);
	
	x_r = (double *) malloc(N*sizeof(double));
	x_i = (double *) malloc(N*sizeof(double));
	y_r = (double *) malloc(N*sizeof(double));
	y_i = (double *) malloc(N*sizeof(double));
	
	Initial(x_r, x_i, N);
	t1 = clock();
	FFTr3(x_r, x_i, y_r, y_i, N);
	t2 = clock();
	printf("Fast FT3: %f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
//	Print_Complex_Vector(y_r, y_i, N);
	
	return 0;
} 

int Initial(double *x, double *y, int N)
{
	int n;
	for(n=0;n<N;++n)
	{
		x[n] = n;
		y[n] = 0.0;
	}
}

//int Print_Complex_Vector(double *y_r, double *y_i, int N)
//{
//	int n;
//	for(n=0;n<N;++n)
//	{
//		if (y_i[n] >=0) printf("%d : %f +%f i\n", n, y_r[n], y_i[n]);
//		else printf("%d : %f %f i\n", n, y_r[n], y_i[n]);
//	}
//	return 0;
//}

int Generate_N(int p, int q, int r)
{
	int N = 1;
	for(;p>0;p--) N*=2;
	for(;q>0;q--) N*=3;
	for(;r>0;r--) N*=5;
	return N;
}

int FFTr3(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	// input : x = x_r + i * x_i
	// output: y = y_r + i * y_i
	int k,n;
	for(n=0;n<N;++n)
	{
		y_r[n] = x_r[n];
		y_i[n] = x_i[n];
	}
	int i, j, M;
	double t_r, t_i;
	i = j = 0;
	while(i < N)
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
		M = N/3;
		while(j >= 2*M & M > 0)
		{
			j = j - 2 * M;
			M = M / 3;
		}
		j = j + M;		
		i = i + 1;
	}
	// Butterfly structure
	double theta, w_br, w_bi, w_cr, w_ci, a_r, a_i, b_r, b_i, c_r, c_i, s;
	int l;
	n = 3;
	s = sqrt(3)/2;
	while(n <= N)
	{
		for(k=0;k<n/3;k++)
		{
			theta = -2.0*k*M_PI/n;
			w_br = cos(theta);
			w_bi = sin(theta);
			w_cr = cos(2*theta);
			w_ci = sin(2*theta);
			for(i=k;i<N;i+=n)
			{
				j = i + n/3;
				l = j + n/3;
				a_r = y_r[i];
                a_i = y_i[i];
				b_r = w_br * y_r[j] - w_bi * y_i[j];
				b_i = w_br * y_i[j] + w_bi * y_r[j];
				c_r = w_cr * y_r[l] - w_ci * y_i[l];
				c_i = w_cr * y_i[l] + w_ci * y_r[l];
				
                y_r[i] = a_r + b_r + c_r; 
                y_i[i] = a_i + b_i + c_i;
                y_r[j] = a_r - (b_r + c_r)/2 + (b_i - c_i) * s;
                y_i[j] = a_i - (b_i + c_i)/2 - (b_r - c_r) * s;
                y_r[l] = a_r - (b_r + c_r)/2 - (b_i - c_i) * s;
                y_i[l] = a_i - (b_i + c_i)/2 + (b_r - c_r) * s;

                
			}
		}
		n = n * 3;
	}
	
}


