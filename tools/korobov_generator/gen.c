#include <stdio.h>
#include <math.h>
#include <malloc/malloc.h>

int opt_coeff1(int kp, int m) 
{
	double c, p;
	int iz1, iz, is, ik, i, d;
	double h, b, hm;
	int *a =  (int *)malloc(sizeof(int)*m);
	
	for (iz = 1; iz < (kp-1)/2; iz++) {
	
		h = 0;

		a[0] = 1;
		for (i = 1; i < m; i++) 
			a[i] = ((a[i-1] * iz) % kp);
		
		for (ik = 1; ik <= (kp-1)/2; ik++) {
			b = 1;
			for (is = 0; is < m; is++) {
				d = ((a[is] * ik) % kp);
				b = b * (kp - 2*d)/kp;
			}
			h = h + b*b;
		}
		
		if (iz == 1) {
			hm = h;
			iz1 = iz;
		}
//		if ((h - hm) < 1e-10) printf("!a = %d\n", iz);
		if (h < hm) {
			hm = h;
			iz1 = iz;
//			printf("!!!%d!!!, %.15f \n", iz, hm);			
		}
	}
	
	free(a);
	return iz1;	
}

int opt_coeff2(int kp1, int kp2, int m, int izz) 
{
	double c, p;
	int iz1, iz, is, ik, i;
	long long d, kp = (long long)kp1 * kp2;
	double h, b, hm;
	long long *a =  (long long *)malloc(sizeof(long long)*m);
	long long *zz = (long long *)malloc(sizeof(long long)*m);
	
	zz[0] = kp2;
	for (i = 1; i < m; i++) 
		zz[i] = ((zz[i-1] * izz) % kp);
		
	for (iz = 1; iz < kp2; iz++) {
	
		h = 0;

		a[0] = kp1;
		for (i = 1; i < m; i++) 
			a[i] = ((a[i-1] * iz) % kp);
		
		for (ik = 1; ik <= (kp-1)/2; ik++) {
			b = 1;
			for (is = 0; is < m; is++) {
				d = (((a[is] + zz[is]) * ik) % kp);
				b = b * (kp - 2*d)/kp;
			}
			h = h + b*b;
		}
		
		if (iz == 1) {
			hm = h;
			iz1 = iz;
		}
//		if ((h - hm) < 1e-10) printf("!b = %d\n", iz);
		if (h < hm) {
			hm = h;
			iz1 = iz;
//			printf("!!%d!!, %.15f \n", iz, hm);			
		}
	}
	
	free(a);
	free(zz);
	return iz1;	
}

int main(int argc, char* argv[])
{
	int s;
	int p1, p2, p, i, j;
	long long d, f;
	int a, b;

	if (argc == 3) {
		p = atoi(argv[2]);
		s = atoi(argv[1]);
		
		printf("p = %d, s = %d\n", p,  s);
	
		a = opt_coeff1(p, s);
		printf("a = %d\n", a);
		
		d = 1;
		for (i = 1; i < s; i++) {
			d = d * a;
			d = d % p;
			printf("%d ", d);
		}
		printf("\n");
	}
	else if (argc == 4) {
		p1 = atoi(argv[2]);
		p2 = atoi(argv[3]);
		s = atoi(argv[1]);
		
		printf("%d, %d, %d\n", p1, p2, s);
	
		a = opt_coeff1(p1, s);
		b = opt_coeff2(p1, p2, s, a);
		printf("a = %d, b = %d\n", a, b);
		
		i = 1;
		while ((((long long)i*p1*p2+1)%(p1+p2)) != 0) 
			i++;
			
		i = ((long long)i*p1*p2+1)/(p1+p2);
		d = ((long long)p2*i)%(p1*p2);
		f = ((long long)p1*i)%(p1*p2);
		printf("%d\n", i);
		for (i = 1; i < s; i++) {
			d = d * a;
			d = d % (p1*p2);
			f = f * b;
			f = f % (p1*p2);
			printf("%d ", (d+f)%(p1*p2));
		}
		printf("\n");
		
		
		
	}	
	else printf("arg1 - s\narg2 - p1\n[arg3] - p2\n");
/*	hm = (2*hm + 1) * pow(3, m) / kp;
	printf("H = %.15f, z = %d\n", hm, iz1);*/
}
