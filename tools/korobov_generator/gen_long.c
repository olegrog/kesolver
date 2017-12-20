#include <stdio.h>
#include <math.h>
#include <malloc.h>

#define MAXINT 4294967296
#define MLONG 12


struct _verylong {
	unsigned long long *a;
	char n;
};

typedef struct _verylong verylong;

int m = 8;
int kp= 24041;
verylong *h, *b, *hm, *res, *w;
int *a;

void multx(verylong *a, unsigned long long x, verylong *res)
{
	int i, j = 0;
	unsigned long long w;
	for	(i = 0; i < a->n; i++) {
		w = a->a[i] * x + j;
		res->a[i] = w % MAXINT;
		j =  w / MAXINT;
	}
	if (j!=0) { res->n = a->n + 1; res->a[res->n-1] = j; }
	else res->n = a->n;
}

void multa(verylong *a, verylong *b, verylong *res)
{
	int i, j;
	unsigned long long w;
	for (i = 0; i < MLONG; i++) 
		res->a[i] = 0;

	for	(i = 0; i < a->n; i++) {
		for	(j = 0; j < b->n; j++) {
			w = b->a[j] * a->a[i];
			res->a[i+j] += w % MAXINT;
			res->a[i+j+1] +=  w / MAXINT;
		}
	}
	if (res->a[a->n+b->n-1] != 0) res->n = a->n + b->n;
	else res->n = a->n + b->n - 1;
}

void add(verylong *a, verylong *b, verylong *res)
{
	int i, j;
	unsigned long long w;
	verylong *min, *max;
	
	if (a->n < b->n) { min = a; max = b; }
	else { min = b; max = a; }

	j = 0;	
	for (i = 0; i < min->n; i++) {
		w = a->a[i] + b->a[i] + j;
		res->a[i] = w % MAXINT;
		j = w / MAXINT;
	}
	for (i = min->n; i < max->n; i++) {
		w = j + max->a[i];
		res->a[i] = w % MAXINT;
		j = w / MAXINT;
	}	
	
	if (j!=0) { res->n = max->n + 1; res->a[res->n] = j; }
	else res->n = max->n;
}

void print(verylong *a) {
	int i;
	for (i = 0; i < a->n; i++) 
		printf("%d ", a->a[i]);
	printf("\n");
}

int comp(verylong *a, verylong *b) 
{
	int i;
	if (a->n < b->n) return 0;
	else if (a->n > b->n) return 1;
	else {
		for (i = a->n-1; i >= 0; i--) 	
			if (a->a[i] < b->a[i]) 
				return 0; // a < b 
			else if (a->a[i] > b->a[i])
				return 1; // b > a
		return 1;   // a = b
	}
}

void H(int iz)
{
	int is, ik, i, d;
	
	h->a[1] = 0;	
	h->n = 0;
		
	a[0] = 1;
	for (i = 1; i < m; i++) 
		a[i] = ((a[i-1] * iz) % kp);
		
	for (ik = 1; ik <= (kp-1)/2; ik++) {
	
		b->a[0] = 1;				
		b->n = 1;
						
		for (is = 0; is < m; is++) {
			d = ((a[is] * ik) % kp);
				
			multx(b, abs(kp - 2*d), res);
			w = b;
			b = res;
			res = w;
		}
			
		multa(b, b, res);
		add(h, res, b);
		w = b;
		b = h;
		h = w;
//		h = h + b*b;
	}
	printf("%d ", iz);
	print(h);
}

int main(int argc, char* argv[])
{

	double c, p;
	int iz1, izz, i;

	h = (verylong *)malloc(sizeof(verylong));
	h->a = (long long *)malloc(sizeof(long long) * MLONG);
	hm = (verylong *)malloc(sizeof(verylong));
	hm->a = (long long *)malloc(sizeof(long long) * MLONG);
	b = (verylong *)malloc(sizeof(verylong));
	b->a = (long long *)malloc(sizeof(long long) * MLONG);
	res = (verylong *)malloc(sizeof(verylong));
	res->a = (long long *)malloc(sizeof(long long) * MLONG);

	
	if (argc > 1) kp = atoi(argv[1]);
	if (argc == 3) m = atoi(argv[2]);
	
	int *kq = (int *)malloc(sizeof(int)*m);
	a = (int *)malloc(sizeof(int)*m);
	
	printf("%d, %d\n", kp,  m);
	

	//H(17441);
	//H(1629);
	H(2189);
		
		
	/*	if (iz == 1) {
			for (i = 0; i < MLONG; i++) 
				hm->a[i] = h->a[i];
			hm->n = h->n;
		
			iz1 = iz;
		}
		if (comp(h, hm) == 0) {
			for (i = 0; i < MLONG; i++) 
				hm->a[i] = h->a[i];
			hm->n = h->n;	
			iz1 = iz;
			printf("!!!%d!!! \n", iz);			
		}
	}
	kq[0] = 1;
	for (is = 1; is < m; is++) 
		kq[is] = kq[is-1]*iz1;
	for  (is = 2; is < m; is++) 
		kq[is] = kq[is] - kq[is]/kp*kp;
	for (is = 0; is < m; is++)
		printf("%d ", kq[is]);*/
	printf("\n");
}
