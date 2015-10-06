#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
int cti(char *a,int *x,int N);
int itc(char *a,int *x,int N);
int FFT(int *x,int *y,int N);
int iFFT(int *x,int *y,int N);

int main()
{
	int i,N=32,Na,Nb,Nc,x[32],X[32],y[32],Y[32],z[32],Z[32];
	char a[32],b[32],c[32];
	printf("a = ");
	scanf("%s",a);
	printf("b = ");
	scanf("%s",b);
	Na=strlen(a);
	Nb=strlen(b);
	cti(a,x,Na);
	cti(b,y,Nb);
	clock_t t1,t2;
	
	t1=clock();
	FFT(x,X,N);
	FFT(y,Y,N);
	for(i=0;i<32;i++)
	{
		Z[i]=X[i]*Y[i];
		Z[i]=Z[i]%1409;
	}
	iFFT(Z,z,N);
	for(i=0;i<32;i++)
	{
		z[i]=z[i]*1365;
		z[i]=z[i]%1409;
	}
	for(i=0;i<32;i++)
	{
		if(z[i]>9)
		{
			z[i+1]+=z[i]/10;
			z[i]=z[i]%10;
		}
	}
	for(i=0;i<32;i++) printf("%d",z[i]);
	printf("\n");
	if(z[Na+Nb-1]==0) Nc=Na+Nb-1;
	else Nc=Na+Nb;
	itc(c,z,Nc);
	printf("%s\n",c);
	system("pause");
	t2=clock();
//	for(i=0;i<N;i++) printf("%d\n",y[i]);

	system("pause");
	return 0;
}

int cti(char *a,int *x,int N)
{
	int i;
	char t;
	for(i=0;i<N/2;i++)
	{
		t=a[N-1-i];
		a[N-1-i]=a[i];
		a[i]=t;
	}
	for(i=0;i<N;i++)
		a[i]-=48;
	for(i=N;i<32;i++)
		a[i]=0;
	for(i=0;i<32;i++)
	{
		x[i]=a[i];
	}
	return 0;
}

int itc(char *a,int *x,int N)
{
	int i;
	char t;
	for(i=0;i<32;i++)
	{
		a[i]=x[i];
	}
	for(i=0;i<N;i++)
		a[i]+=48;
	for(i=0;i<N/2;i++)
	{
		t=a[N-1-i];
		a[N-1-i]=a[i];
		a[i]=t;
	}
	return 0;
}

int FFT(int *x,int *y,int N)
{
	if(N==1)
	{
		y[0]=x[0];
		return 0;
	}
		int k,n,i,j,l;
		int *u,*v,w,q;
		u=(int *)malloc(N*sizeof(int));
		v=(int *)malloc(N*sizeof(int));
		for(i=0;i<2;i++)
		{
			for(j=0;j<N/2;j++)
			{
				u[j+i*N/2]=x[j*2+i];
			}
		}
		for(i=0;i<2;i++) FFT(u+i*N/2,v+i*N/2,N/2);
		
		q=1;
		for(i=0;i<32/N;i++)
		{
			q=q*141;
			q=q%1409;
		}
		w=1;
		y[0]=v[0]+v[N/2];
		y[0]=y[0]%1409;
		y[N/2]=v[0]+1408*v[N/2];
		y[N/2]=y[N/2]%1409;
		for(int k=1;k<N/2;k++)
		{
			w=w*q;
			w=w%1409;
			y[k]=v[k]+w*v[k+N/2];
			y[k]=y[k]%1409;
			l=w*1408;
			l=l%1409;
			y[k+N/2]=v[k]+l*v[k+N/2];
			y[k+N/2]=y[k+N/2]%1409;
		}
		free(u);
		free(v);
	return 0;
}

int iFFT(int *x,int *y,int N)
{
	if(N==1)
	{
		y[0]=x[0];
		return 0;
	}
		int k,n,i,j,l;
		int *u,*v,w,q;
		u=(int *)malloc(N*sizeof(int));
		v=(int *)malloc(N*sizeof(int));
		for(i=0;i<2;i++)
		{
			for(j=0;j<N/2;j++)
			{
				u[j+i*N/2]=x[j*2+i];
			}
		}
		for(i=0;i<2;i++) iFFT(u+i*N/2,v+i*N/2,N/2);
		
		q=1;
		for(i=0;i<32/N;i++)
		{
			q=q*10;
			q=q%1409;
		}
		w=1;
		y[0]=v[0]+v[N/2];
		y[0]=y[0]%1409;
		y[N/2]=v[0]+1408*v[N/2];
		y[N/2]=y[N/2]%1409;
		for(k=1;k<N/2;k++)
		{
			w=w*q;
			w=w%1409;
			y[k]=v[k]+w*v[k+N/2];
			y[k]=y[k]%1409;
			l=w*1408;
			l=l%1409;
			y[k+N/2]=v[k]+l*v[k+N/2];
			y[k+N/2]=y[k+N/2]%1409;
		}
		free(u);
		free(v);
	return 0;
}
