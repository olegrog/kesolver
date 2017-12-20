#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

struct resdata_ {			// структура выходных данных фукции вычисления оптимальных коэффициентов
	int a;					// оптимальный коэффициент 
	double H;				// H(a) характеризует качество оптимального коэффициента
};

typedef struct resdata_ resdata;

//	Функция вычисления оптимальных коэффициентов
//	вычисляет значение функции H(a) для a из [n1, n2) и
//	находит а при котором на данном промежутке H(a) достигает минимума.
//  Эти данные выводятся в структуре resdata.
//  kp - простое число задающее количество узлов
//  m - размерность

resdata opt_coeff1(int kp, int m, int n1, int n2) 
{
	double p;
//  iz счетчик перебора по а
//  iz1 текущий оптимальный коэффициент 
//  is счетчик для вычисления произведения в сумме
//  ik счетчик для вычисления суммы
//  d вспомогательная переменная при счете {ik*iz^is % kp}
	int iz1, iz, is, ik, d;
//	h значение суммы
//	b вспомогательная переменная при счете произведения в сумме
//	hm текущий минимум
	double h, b, hm;
	int *a =  (int *)malloc(sizeof(int)*m);		// вспомогательный массив для {iz^is % kp}
//  структура результата
	resdata res;
	
	for (iz = n1; iz < n2; iz++) {
	
		h = 0; // пока сумма равна 0
		// рассчет вспомогательного массива
		a[0] = 1;
		for (is = 1; is < m; is++) 
			a[is] = (((long long)a[is-1] * iz) % kp);
		// рассчет суммы		
		for (ik = 1; ik <= (kp-1)/2; ik++) {
			// рассчет произведения в сумме
			b = 1;
			for (is = 0; is < m; is++) {
				d = (((long long)a[is] * ik) % kp);
				b = b * (kp - 2*d)/kp;
			}
			h = h + b*b;
		}
		
		// если текущий минимум еще не задан или он действительно минимум
		if ((iz == n1)||(h < hm)) {
			hm = h;
			iz1 = iz;
		}
	}
	
	free(a);	
	res.a = iz1;	// заполнение структуры с результатом
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

	p = atoi(argv[2]);	// считывание нач. данных		
	s = atoi(argv[1]);
		
	MPI_Init(&argc, &argv);				// распараллеливание при помощи MPI
	MPI_Comm_size(MPI_COMM_WORLD, &m);
	MPI_Comm_rank(MPI_COMM_WORLD, &n);

	if (n == 0) {
		// контролирующий поток

		res.H = 0;
		// в цикле ждем результата от каждого потока
		for (i = 1; i < m; i++) {
			MPI_Recv(&res1, sizeof(resdata), MPI_CHAR, i, 0, MPI_COMM_WORLD, &stat); 
			if ((res.H > res1.H) || (i == 1))	// определяем глобальный минимум
				res = res1;
		}

		a = res.a;
		printf("s = %d, p = %d, a = %d\n", s, p, a);
		// счет всех коэффициентов из a
		d = 1;
		for (i = 1; i < s; i++) {
			d = d * a;
			d = d % p;
			printf("%lld ", d);
		}
		printf("\n");
	}
	else {
		// поток счета
		
		n1 = ((p+1)/2 * (n-1))/(m-1);		// разделение поровну всего участка на куски для каждого потока 
		n2 = ((p+1)/2 * n)/(m-1);

		res = opt_coeff1(p, s, n1, n2);		// вычисление минимума на данном участке

		MPI_Send(&res, sizeof(resdata), MPI_CHAR, 0, 0, MPI_COMM_WORLD); // посылаем данные контролирующему потоку
	}
	MPI_Finalize();
}
