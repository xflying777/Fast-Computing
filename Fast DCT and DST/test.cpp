#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main()
{
	int p, q, r, N, L;
	N = 5;
	L = 2*N+2;
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
	
	printf("p = %d , q = %d , r = %d \n", p, q , r);
	return 0;
}
