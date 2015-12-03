#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int FFT(double *x_r, double *x_i, double *y_r, double *y_i, int p, int q, int r);
int Initial(double *x, double *y, int N);
int Print_Complex_Vector(double *y_r, double *y_i, int N);
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
	FFT(x_r, x_i, y_r, y_i, p, q, r);
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

int FFT(double *x_r, double *x_i, double *y_r, double *y_i, int p, int q, int r)
{
	int i, j, pp, qq, bse, N, M;
	N = 1;
	N = Generate_N(p, q, r);
/*
	for(i=0;i<N;i++)
	{
		y_r[i]=x_r[i];
		y_i[i]=x_i[i];
	}
	j=0;
	for(i=1;i<N;i++)
	{
		M=N;
		if(p==0)
		{
			if(q==0)
			{
				//(0,0,r)
				j=j+M/5;
				while(j>=M & M>1)
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
		y_r[i]=x_r[j];
		y_i[i]=x_i[j];
	}
*/	
	if(p==0)
	{
		if(q==0)
		{
			if(r==0)
			{
				y_r[0]=x_r[0];
				y_i[0]=x_i[0];
				return 0;
			}
			else
			{
				int k,n,i,j;
				double *u_r,*u_i,*v_r,*v_i,w_r,w_i,t_r,t_i,s_r,s_i,r_r,r_i,q_r,q_i,tmp;
				u_r=(double *)malloc(N*sizeof(double));
				u_i=(double *)malloc(N*sizeof(double));
				v_r=(double *)malloc(N*sizeof(double));
				v_i=(double *)malloc(N*sizeof(double));
				for(i=0;i<5;i++)
				{
					for(j=0;j<N/5;j++)
					{
						u_r[j+i*N/5]=x_r[j*5+i];
						u_i[j+i*N/5]=x_i[j*5+i];
					}
				}
				
				for(i=0;i<5;i++) FFT(u_r+i*N/5,u_i+i*N/5,v_r+i*N/5,v_i+i*N/5,0,0,r-1);

				t_r=cos(-2.0*M_PI/5);
				t_i=sin(-2.0*M_PI/5);
				q_r=cos(-2.0*M_PI/N); 
				q_i=sin(-2.0*M_PI/N);
				s_r=1.0;
				s_i=0.0;
				y_r[0]=v_r[0];
				y_i[0]=v_i[0];
				for(int n=1;n<5;n++)
				{
					y_r[0]=y_r[0]+v_r[n*N/5];
					y_i[0]=y_i[0]+v_i[n*N/5];
				}
				r_r=1.0;
				r_i=0.0;
				for(int i=1;i<5;i++)
				{
					y_r[i*N/5]=v_r[0];
					y_i[i*N/5]=v_i[0];
					tmp=r_r;
					r_r=r_r*t_r-r_i*t_i;
					r_i=tmp*t_i+r_i*t_r;
					w_r=1.0;
					w_i=0.0;
					for(int n=1;n<5;n++)
					{
						tmp=w_r;
						w_r=w_r*r_r-w_i*r_i;
						w_i=tmp*r_i+w_i*r_r;
						y_r[i*N/5]=y_r[i*N/5]+w_r*v_r[n*N/5]-w_i*v_i[n*N/5];
						y_i[i*N/5]=y_i[i*N/5]+w_r*v_i[n*N/5]+w_i*v_r[n*N/5];
					}
				}
				
				for(int k=1;k<N/5;k++)
				{
					tmp=s_r;
					s_r=s_r*q_r-s_i*q_i;
					s_i=tmp*q_i+s_i*q_r;
					y_r[k]=v_r[k];
					y_i[k]=v_i[k];
					w_r=1;
					w_i=0;
					for(int n=1;n<5;n++)
					{
						tmp=w_r;
						w_r=w_r*s_r-w_i*s_i;
						w_i=tmp*s_i+w_i*s_r;
						y_r[k]=y_r[k]+w_r*v_r[k+n*N/5]-w_i*v_i[k+n*N/5];
						y_i[k]=y_i[k]+w_r*v_i[k+n*N/5]+w_i*v_r[k+n*N/5];
					}
					r_r=1;
					r_i=0;
					for(int i=1;i<5;i++)
					{
						y_r[k+i*N/5]=v_r[k];
						y_i[k+i*N/5]=v_i[k];
						tmp=r_r;
						r_r=r_r*t_r-r_i*t_i;
						r_i=tmp*t_i+r_i*t_r;
						w_r=1;
						w_i=0;
						for(int n=1;n<5;n++)
						{
							tmp=w_r;
							w_r=w_r*s_r-w_i*s_i;
							w_i=tmp*s_i+w_i*s_r;
							tmp=w_r;
							w_r=w_r*r_r-w_i*r_i;
							w_i=tmp*r_i+w_i*r_r;
							y_r[k+i*N/5]=y_r[k+i*N/5]+w_r*v_r[k+n*N/5]-w_i*v_i[k+n*N/5];
							y_i[k+i*N/5]=y_i[k+i*N/5]+w_r*v_i[k+n*N/5]+w_i*v_r[k+n*N/5];
						}
					}
				}
				free(u_r);
				free(u_i);
				free(v_r);
				free(v_i);
			}
		}
		else//N剩3的次方和5的次方 
		{
			int k,n,i,j;
			double *u_r,*u_i,*v_r,*v_i,w_r,w_i,t_r,t_i,s_r,s_i,r_r,r_i,q_r,q_i,tmp;
			u_r=(double *)malloc(N*sizeof(double));
			u_i=(double *)malloc(N*sizeof(double));
			v_r=(double *)malloc(N*sizeof(double));
			v_i=(double *)malloc(N*sizeof(double));
			for(i=0;i<3;i++)
			{
				for(j=0;j<N/3;j++)
				{
					u_r[j+i*N/3]=x_r[j*3+i];
					u_i[j+i*N/3]=x_i[j*3+i];
				}
			}
			
			for(i=0;i<3;i++) FFT(u_r+i*N/3,u_i+i*N/3,v_r+i*N/3,v_i+i*N/3,0,q-1,r);

			t_r=cos(-2.0*M_PI/3);
			t_i=sin(-2.0*M_PI/3);
			q_r=cos(-2.0*M_PI/N); 
			q_i=sin(-2.0*M_PI/N);
			s_r=1.0;
			s_i=0.0;
			y_r[0]=v_r[0];
			y_i[0]=v_i[0];
			for(int n=1;n<3;n++)
			{
				y_r[0]=y_r[0]+v_r[n*N/3];
				y_i[0]=y_i[0]+v_i[n*N/3];
			}
			r_r=1.0;
			r_i=0.0;
			for(int i=1;i<3;i++)
			{
				y_r[i*N/3]=v_r[0];
				y_i[i*N/3]=v_i[0];
				tmp=r_r;
				r_r=r_r*t_r-r_i*t_i;
				r_i=tmp*t_i+r_i*t_r;
				w_r=1.0;
				w_i=0.0;
				for(int n=1;n<3;n++)
				{
					tmp=w_r;
					w_r=w_r*r_r-w_i*r_i;
					w_i=tmp*r_i+w_i*r_r;
					y_r[i*N/3]=y_r[i*N/3]+w_r*v_r[n*N/3]-w_i*v_i[n*N/3];
					y_i[i*N/3]=y_i[i*N/3]+w_r*v_i[n*N/3]+w_i*v_r[n*N/3];
				}
			}
			
			for(int k=1;k<N/3;k++)
			{
				tmp=s_r;
				s_r=s_r*q_r-s_i*q_i;
				s_i=tmp*q_i+s_i*q_r;
				y_r[k]=v_r[k];
				y_i[k]=v_i[k];
				w_r=1;
				w_i=0;
				for(int n=1;n<3;n++)
				{
					tmp=w_r;
					w_r=w_r*s_r-w_i*s_i;
					w_i=tmp*s_i+w_i*s_r;
					y_r[k]=y_r[k]+w_r*v_r[k+n*N/3]-w_i*v_i[k+n*N/3];
					y_i[k]=y_i[k]+w_r*v_i[k+n*N/3]+w_i*v_r[k+n*N/3];
				}
				r_r=1;
				r_i=0;
				for(int i=1;i<3;i++)
				{
					y_r[k+i*N/3]=v_r[k];
					y_i[k+i*N/3]=v_i[k];
					tmp=r_r;
					r_r=r_r*t_r-r_i*t_i;
					r_i=tmp*t_i+r_i*t_r;
					w_r=1;
					w_i=0;
					for(int n=1;n<3;n++)
					{
						tmp=w_r;
						w_r=w_r*s_r-w_i*s_i;
						w_i=tmp*s_i+w_i*s_r;
						tmp=w_r;
						w_r=w_r*r_r-w_i*r_i;
						w_i=tmp*r_i+w_i*r_r;
						y_r[k+i*N/3]=y_r[k+i*N/3]+w_r*v_r[k+n*N/3]-w_i*v_i[k+n*N/3];
						y_i[k+i*N/3]=y_i[k+i*N/3]+w_r*v_i[k+n*N/3]+w_i*v_r[k+n*N/3];
					}
				}
			}
			free(u_r);
			free(u_i);
			free(v_r);
			free(v_i);
		}
	}
	else
	{
		int k,n,i,j;
		double *u_r,*u_i,*v_r,*v_i,w_r,w_i,q_r,q_i,tmp;
		u_r=(double *)malloc(N*sizeof(double));
		u_i=(double *)malloc(N*sizeof(double));
		v_r=(double *)malloc(N*sizeof(double));
		v_i=(double *)malloc(N*sizeof(double));
		for(i=0;i<2;i++)
		{
			for(j=0;j<N/2;j++)
			{
				u_r[j+i*N/2]=x_r[j*2+i];
				u_i[j+i*N/2]=x_i[j*2+i];
			}
		}
		
		for(i=0;i<2;i++) FFT(u_r+i*N/2,u_i+i*N/2,v_r+i*N/2,v_i+i*N/2,p-1,q,r);

		q_r=cos(-2.0*M_PI/N); 
		q_i=sin(-2.0*M_PI/N);
		w_r=1.0;
		w_i=0.0;
		y_r[0]=v_r[0]+v_r[N/2];
		y_i[0]=v_i[0]+v_i[N/2];
		y_r[N/2]=v_r[0]-v_r[N/2];
		y_i[N/2]=v_i[0]-v_i[N/2];
		for(int k=1;k<N/2;k++)
		{
			tmp=w_r;
			w_r=w_r*q_r-w_i*q_i;
			w_i=tmp*q_i+w_i*q_r;
			y_r[k]=v_r[k]+w_r*v_r[k+N/2]-w_i*v_i[k+N/2];
			y_i[k]=v_i[k]+w_r*v_i[k+N/2]+w_i*v_r[k+N/2];
			y_r[k+N/2]=v_r[k]-w_r*v_r[k+N/2]+w_i*v_i[k+N/2];
			y_i[k+N/2]=v_i[k]-w_r*v_i[k+N/2]-w_i*v_r[k+N/2];
		}
		free(u_r);
		free(u_i);
		free(v_r);
		free(v_i);
	}
	return 0;

}
