#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void FFT(double *x_r, double *x_i, double *y_r, double *y_i, double *z_r, double *z_i, int N);
void Initial(double *x, double *y, int N);
int Print_Complex_Vector(double *y_r, double *y_i, int N);

int main()
{
	int i, m, N;
	double *y_r, *y_i, *x_r, *x_i, *z_r, *z_i, everage, *t;
	clock_t t1, t2;

	
	printf("Please input N (And such that 2N + 2 = 2^p*3^q*5^r) = ");
	scanf("%d",&N);
	printf("N = %d \n", N);
	
	m = 3;
	x_r = (double *) malloc(N*sizeof(double));
	x_i = (double *) malloc(N*sizeof(double));	
	z_r = (double *) malloc((2*N+2)*sizeof(double));
	z_i = (double *) malloc((2*N+2)*sizeof(double));	
	y_r = (double *) malloc((2*N+2)*sizeof(double));
	y_i = (double *) malloc((2*N+2)*sizeof(double));
	t = (double *) malloc(m*sizeof(double));
	
	Initial(x_r, x_i, N);
	for(i=0;i<m;i++)
	{
		t1 = clock();
		FFT(x_r, x_i, y_r, y_i, z_r, z_i, N);
		t2 = clock();
		t[i] = 1.0*(t2-t1)/CLOCKS_PER_SEC;
	}
	everage = 0;
	for(i=0;i<m;i++)
	{
		everage = everage + t[i];	
	}
	everage = everage / m;
	printf(" everage time: %f \n", everage);
	printf(" data[%d] = %f , data[%d] = %f \n", 0, z_r[0], N-1, z_r[N-1]);
/*	for(i=0;i<N;i++)
	{
		printf("%d: %f \n", i, z_r[i]);
	}
*/
//	Print_Complex_Vector(y_r, y_i, N);
	
	return 0;
} 

void Initial(double *x, double *y, int N)
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


// Consider Cx[k] = \sum_{n=0}^{N-1} 2x[n]cos(pi*k*(2n+1)/2N) for 0 <= k < N
void FFT(double *x_r, double *x_i, double *y_r, double *y_i,double *z_r, double *z_i, int N)
{
	int i, j, p, q, r, pp, qq, bse, L, M;
	L = 2*N + 2;
	// Prime factorization of L about 2,3,5
	p = q = r = 0;
	while(L > 1)
	{
		if(L % 2 == 0)
		{
			L = L/2;
			p = p + 1;
		}
		if(L % 3 == 0)
		{
			L = L/3;
			q = q + 1;
		}
		if(L % 5 == 0)
		{
			L = L/5;
			r = r + 1;
		}
	}
	
	L = 2*N + 2;
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
	
	for(i=0;i<L;i++)
	{
		z_r[i] = y_r[i];
		z_i[i] = y_i[i];
	}
	// Bit-reverse
	j=0;
	for(i=1;i<L;i++)
	{
		M = L;
		if(p==0)
		{
			if(q==0)
			{
				//(0,0,r)
				j = j + M/5;
				while(j >= M & M > 1)
				{
					j=j-M;
					M=M/5;
					j=j+M/5;
				}
			}
			else
			{
				if(r==0)
				{
					//(0,q,0)
					j=j+M/3;
					while(j>=M)
					{
						j=j-M;
						M=M/3;
						j=j+M/3;
					}
				}
				else
				{
					//(0,q,r);
					bse=3;
					qq=q;
					j=j+M/bse;
					while(j>=M)
					{
						j=j-M;
						M=M/bse;
						qq=qq-1;
						if(qq==0) bse=5;
						j=j+M/bse;
					}
				}
			}
		}
		else
		{
			if(q==0)
			{
				if(r==0)
				{
					//(p,0,0)
					j=j+M/2;
					while(j>=M)
					{
						j=j-M;
						M=M/2;
						j=j+M/2;
					}
				}
				else
				{
					//(p,0,r)
					bse=2;
					pp=p;
					j=j+M/bse;
					while(j>=M)
					{
						j=j-M;
						M=M/bse;
						pp=pp-1;
						if(pp==0) bse=5;
						j=j+M/bse;
					}
				}
			}
			else
			{
				if(r==0)
				{
					//(p,q,0)
					bse=2;
					pp=p;
					j=j+M/bse;
					while(j>=M)
					{
						j=j-M;
						M=M/bse;
						pp=pp-1;
						if(pp==0) bse=3;
						j=j+M/bse;
					}
				}
				else
				{
					//(p,q,r)
					pp=p;
					qq=q;
					j=j+M/2;
					while(j>=M)
					{
						j=j-M;
						M=M/2;
						pp=pp-1;
						if(pp>0) j=j+M/2;
						else
						{
							j=j+M/3;
							while(j>=M)
							{
								j=j-M;
								M=M/3;
								qq=qq-1;
								if(qq>0) j=j+M/3;
								else
								{
									j=j+M/5;
									while(j>=M)
									{
										j=j-M;
										M=M/3;
										j=j+M/5;
									}
								}
							}
						}
					}
				}
			}
		}
		y_r[i]=z_r[j];
		y_i[i]=z_i[j];
	}

	
	// Butterfly structure			 
	int P, Q, R;
	P = pow(2,p);
	Q = pow(3,q);		

	//p != 0
	if(p > 0)
	{
		double theta, w_br, w_bi, a_r, a_i;
		int k, n, pp;
		pp = 0;
		n = 2;
		while(pp < p & p > 0)
		{
			for(k=0;k<n/2;k++)
			{
				theta = -2.0*k*M_PI/n;
				w_br = cos(theta);
				w_bi = sin(theta);
				for(i=k;i<L;i+=n)
				{
					j = i + n/2;
					a_r = w_br * y_r[j] - w_bi * y_i[j];
					a_i = w_br * y_i[j] + w_bi * y_r[j];
					
	
					y_r[j] = y_r[i] - a_r;
					y_i[j] = y_i[i] - a_i;
					y_r[i] = y_r[i] + a_r;
					y_i[i] = y_i[i] + a_i;
	
				}
			}
			n = n * 2;
			pp = pp + 1;
		}
	}
	
	// q != 0
	if(q > 0)
	{				
		double theta, w_br, w_bi, w_cr, w_ci, a_r, a_i, b_r, b_i, c_r, c_i, s;
		int l, k, n, m, qq;
		qq = 0;
		n = 3*P;
		s = sqrt(3)/2;
		while(qq < q & q > 0)
		{
			for(k=0;k<n/3;k++)
			{
				theta = -2.0*k*M_PI/n;
				w_br = cos(theta);
				w_bi = sin(theta);
				w_cr = cos(2*theta);
				w_ci = sin(2*theta);
				for(i=k;i<L;i+=n)
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
			qq = qq + 1;
		}
	}
	
	if(r > 0)
	{
		double a, b, c, d, theta, w_br, w_bi, w_cr, w_ci, w_dr, w_di, w_er, w_ei, a_r, a_i, b_r, b_i, c_r, c_i, d_r, d_i, e_r, e_i;
		int f, g, h, n, k, rr;
		rr = 0;
		n = 5*Q*P;
		a = cos(2.0*M_PI/5);
		b = sin(2.0*M_PI/5);
		c = -cos(4.0*M_PI/5);
		d = sin(4.0*M_PI/5);
		while(rr < r & r > 0)
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
				
				for(i=k;i<L;i+=n)
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
			rr = rr + 1;
		}
	}
	
	// After fft(y[k]), Y[k] = fft(y[k]), Sx[k] = i*Y[k+1]/2
	int k;
	for(k=0;k<N;k++)
	{
		z_r[k] = -1.0*y_i[k+1]/2;
	}
	
}



	

