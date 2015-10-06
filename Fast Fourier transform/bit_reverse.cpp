#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int bitreverse(int p, int q, int r, double *x_r, double *x_i, double *y_r, double *y_i, int N);
int Initial(double *x, double *y, int N);
int Print_Complex_Vector(double *y_r, double *y_i, int N);
int Generate_N(int p, int q, int r);
int Print_base(double *x, int N, int base);

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
	bitreverse(p,q,r,x_r, x_i, y_r, y_i, N);
	t2 = clock();
	printf("Fast bitreverse: %f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	
//	Print_Complex_Vector(y_r, y_i, N);
	
	return 0;
} 

int initial(int *x, int N)
{
	int i;
	for(i=0;i<N;i++) x[i] = i;
	
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

int Print_Complex_Vector(double *y_r, double *y_i, int N)
{
	int n;
	for(n=0;n<N;++n)
	{
		if (y_i[n] >=0) printf("%d : %f +%f i\n", n, y_r[n], y_i[n]);
		else printf("%d : %f %f i\n", n, y_r[n], y_i[n]);
	}
	return 0;
}

int Print_base(int *x, int N, int base)
{
	int i;
	for(i=0;i<N;i++) printf("base%d [%d]: %d \n", base, i, x[i]);
}

int Generate_N(int p, int q, int r)
{
	int N = 1;
	for(;p>0;p--) 
	N*=2;	
	for(;q>0;q--) 
	N*=3;
	for(;r>0;r--) 
	N*=5;
	return N;
}


int bitreverse(int p, int q, int r, double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	// input : x = x_r + i * x_i
	// output: y = y_r + i * y_i
	int k,n;
	for(n=0;n<N;++n)
	{
		y_r[n] = x_r[n];
		y_i[n] = x_i[n];
	}
	int i, j, base2N, base3N, base5N, M, *base2, *base3, *base5, temp;
	
	base2N = pow(2,p);
	base3N = pow(3,q);
	base5N = pow(5,r);
	
	base2 = (int *) malloc(base2N*sizeof(int));
	base3 = (int *) malloc(base3N*sizeof(int));
	base5 = (int *) malloc(base5N*sizeof(int));
	
	initial(base2, base2N);
	initial(base3, base3N);
	initial(base5, base5N);
	
	i = j = 0;
	if (p>0)
	{
		while(i < base2N)
		{
			if(i < j)
			{
				temp = base2[i];
				base2[i] = base2[j];
				base2[j] = temp;
			}
			M = base2N / 2;
			while(j >=M & M>0)
			{
				j = j - M;
				M = M/2;
			}
			j = j + M;
			i = i + 1;
		}
	}

	if(q>0)
	{
		i = j = 0;
		while(i < base3N)
		{
			if(i < j)
			{
				temp = base3[i];
				base3[i] = base3[j];
				base3[j] = temp;
			}
			M = base3N / 3;
			while(j >= 2*M & M>0)
			{
				j = j - 2*M;
				M = M/3;
			}
			j = j + M;
			i = i + 1;
		}		
	}

	if(r>0)
	{
		i = j = 0;
		while(i < base5N)
		{
			if(i < j)
			{
				temp = base5[i];
				base5[i] = base5[j];
				base5[j] = temp;
			}
			M = base5N / 5;
			while(j >= 4*M & M>0)
			{
				j = j - 4*M;
				M = M/5;
			}
			j = j + M;
			i = i + 1;
		}	
	}

	int *base35, base35N, *base235, base235N;
	base35N = base3N*base5N;
	base235N = base2N*base35N;
	base35 = (int *) malloc(base35N*sizeof(int));
	base235 = (int *) malloc(base235N*sizeof(int));
	for(i=0;i<base35N;i++) base35[i] = 0;
	for(i=0;i<base235N;i++) base235[i] = 0;
	
	for(i=0;i<base3N;i++) 
	{
		base3[i] = base3[i]*base5N;
		base35[i] = base3[i];
	}

	for(i=1;i<base5N;i++) 
	{
		k = i*base3N;
		for(j=0;j<base3N;j++)
		{
			temp = k + j;
			base35[temp] = base35[j] + base5[i];
		}
	}
//	for(i=0;i<base35N;i++) printf("%d: %d \n", i, base35[i]);	
	for (i=0;i<base2N;i++)
	{
		base2[i] = base2[i]*base35N;
		base235[i] = base2[i];
	}
	for(i=1;i<base35N;i++) 
	{
		k = i*base2N;
		for(j=0;j<base2N;j++)
		{
			temp = k + j;
			base235[temp] = base235[j] + base35[i];
		}
	}
	
//	Print_base(base235,base235N,235);
	
}
