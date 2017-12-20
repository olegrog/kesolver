#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <malloc/malloc.h>

struct resdata_ {
	int a;
	double H;
};

typedef struct resdata_ resdata;

resdata opt_coeff1(int kp, int m, int n1, int n2) 
{
	double c, p;
	int iz1, iz, is, ik, i, d;
	double h, b, hm;
	int *a =  (int *)malloc(sizeof(int)*m);
	resdata res;
	
	for (iz = n1; iz < n2; iz++) {
	
		h = 0;

		a[0] = 1;
		for (i = 1; i < m; i++) 
			a[i] = (((long long)a[i-1] * iz) % kp);
		
		for (ik = 1; ik <= (kp-1)/2; ik++) {
			b = 1;
			for (is = 0; is < m; is++) {
				d = (((long long)a[is] * ik) % kp);
				b = b * (kp - 2*d)/kp;
			}
			h = h + b*b;
		}
		
		if (iz == n1) {
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
	res.a = iz1;
	res.H = hm;
	return res;	
}

int main(int argc, char* argv[])
{
	int s;
	int p, i, j, n, m, n1, n2;
	long long d, f;
	int a, b;
	resdata res, res1;
	MPI_Status stat;

	p = atoi(argv[2]);
	s = atoi(argv[1]);
		
	printf("%d, %d\n", p,  s);

	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &m);
	MPI_Comm_rank(MPI_COMM_WORLD, &n);

	if (n == 0) {
		printf("m = %d\n", m); 
		res.H = 1E10;
		for (i = 1; i < m; i++) {
			MPI_Recv(&res1, sizeof(resdata), MPI_CHAR, i, 0, MPI_COMM_WORLD, &stat); 
			if (res.H > res1.H)
				res = res1;
		}

		a = res.a;
		printf("s = %d, p = %d, a = %d\n", s, p, a);
		d = 1;
		for (i = 1; i < s; i++) {
			d = d * a;
			d = d % p;
			printf("%lld ", d);
		}
		printf("\n");
	}
	else {
		
		printf("Child %d\n", n);

	/*	if      (n == 1) { n1 = 1;       n2 = 125002; }
		else if (n == 2) { n1 = 125002;  n2 = 250004; }

		else if (n == 3) { n1 = 250005;  n2 = 500006; }
		else if (n == 4) { n1 = 500006;  n2 = 750007; }
		else if (n == 5) { n1 = 750008;  n2 = 1000008; }
		else if (n == 6) { n1 = 1000008; n2 = 1250008; }

		else if (n == 7) { n1 = 1250009; n2 = 1325009; }
		else if (n == 8) { n1 = 1325009; n2 = 1500009; }*/
		n1 = ((p+1)/2 * (n-1))/(m-1);
		n2 = ((p+1)/2 * n)/(m-1);

		res = opt_coeff1(p, s, n1, n2);
		printf("%d, %f\n", res.a, res.H);


		MPI_Send(&res, sizeof(resdata), MPI_CHAR, 0, 0, MPI_COMM_WORLD);

	}
	MPI_Finalize();
}
