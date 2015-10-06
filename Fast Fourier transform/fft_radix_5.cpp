#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int FFTr5(double *x_r, double *x_i, double *y_r, double *y_i, int N);
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
	FFTr5(x_r, x_i, y_r, y_i, N);
	t2 = clock();
	printf("Fast FT5: %f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
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

int FFTr5(double *x_r, double *x_i, double *y_r, double *y_i, int N)
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
		M = N / 5;
		while(j >= 4 * M & M > 0)
		{
			j = j - 4 * M;
			M = M / 5;
		}
		j = j + M;		
		i = i + 1;
	}
	// Butterfly structure
	double a, b, c, d, theta, w_br, w_bi, w_cr, w_ci, w_dr, w_di, w_er, w_ei, a_r, a_i, b_r, b_i, c_r, c_i, d_r, d_i, e_r, e_i;
	int f, g, h;
	n = 5;
	a = cos(2.0*M_PI/5);
	b = sin(2.0*M_PI/5);
	c = -cos(4.0*M_PI/5);
	d = sin(4.0*M_PI/5);
	while(n <= N)
	{
		for(k=0;k<n/5;k++)
		{
			theta = -2.0*k*M_PI/n;
			w_br = cos(theta);
			w_bi = sin(theta);
			w_cr = cos(2.0*theta);
			w_ci = sin(2.0*theta);
			w_dr = cos(3.0*theta);
			w_di = sin(3.0*theta);
			w_er = cos(4.0*theta);
			w_ei = sin(4.0*theta);
			
			for(i=k;i<N;i+=n)
			{
				j = i + n/5;
				f = j + n/5;
				g = f + n/5;
				h = g + n/5;
				
				a_r = y_r[i];
				a_i = y_i[i];
				b_r = (y_r[j] * w_br - y_i[j] * w_bi);
				b_i = (y_r[j] * w_bi + y_i[j] * w_br);
				c_r = (y_r[f] * w_cr - y_i[f] * w_ci);
				c_i = (y_r[f] * w_ci + y_i[f] * w_cr);
				d_r = (y_r[g] * w_dr - y_i[g] * w_di);
				d_i = (y_r[g] * w_di + y_i[g] * w_dr);
				e_r = (y_r[h] * w_er - y_i[h] * w_ei);
				e_i = (y_r[h] * w_ei + y_i[h] * w_er);
				
				y_r[i] = a_r + b_r + c_r + d_r + e_r; 
				y_i[i] = a_i + b_i + c_i + d_i + e_i;
				y_r[j] = a_r + (a * (b_r + e_r) - c * (c_r + d_r)) + (b * (b_i - e_i) + d * (c_i - d_i));
				y_i[j] = a_i + (a * (b_i + e_i) - c * (c_i + d_i)) - (b * (b_r - e_r) + d * (c_r - d_r));
				y_r[f] = a_r - (c * (b_r + e_r) - a * (c_r + d_r)) + (d * (b_i - e_i) - b * (c_i - d_i));
				y_i[f] = a_i - (c * (b_i + e_i) - a * (c_i + d_i)) - (d * (b_r - e_r) - b * (c_r - d_r));
				y_r[g] = a_r - (c * (b_r + e_r) - a * (c_r + d_r)) - (d * (b_i - e_i) - b * (c_i - d_i));
				y_i[g] = a_i - (c * (b_i + e_i) - a * (c_i + d_i)) + (d * (b_r - e_r) - b * (c_r - d_r));
				y_r[h] = a_r + (a * (b_r + e_r) - c * (c_r + d_r)) - (b * (b_i - e_i) + d * (c_i - d_i));
				y_i[h] = a_i + (a * (b_i + e_i) - c * (c_i + d_i)) + (b * (b_r - e_r) + d * (c_r - d_r));

              
			}
		}
		n = n * 5;
	}
	
}

