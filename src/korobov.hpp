#ifndef _KOROBOV_H_
#define _KOROBOV_H_

#include <cstdlib>
#include <climits>

namespace korobov {

	const int k8[][9] = {
		{    10007, 1,    1905,    6491,    6710,    3611,    4146,    2607,    2863 },
		{    20011, 1,    7191,    2057,    3758,    8928,    5960,   14809,   12988 },
		{    30011, 1,   13167,   26353,    2769,   26069,   14716,   14556,    8606 },
		{    50021, 1,   11281,    7537,   39218,   32534,   11977,    5816,   32765 }, 
		{   100003, 1,   20285,   68883,   49739,   25348,   68757,   93907,   46351 },
		{   200003, 1,   47369,  188507,   54145,  156036,  158419,   37051,   42494 },
		{   300017, 1,   81575,  103565,  136172,  101475,   54078,  262899,  170731 },
		{   500009, 1,   42535,  193663,  307439,  182488,  487373,   37415,  418387 },
		{  1000003, 1,  417564,  171019,  163483,  410620,  615303,  611111,  188073 },
		{  1500007, 1,  413996,  388189,  943278, 1496508,  434758,  733031,  985685 },
		{  2000003, 1,  832685, 1269182, 1228431,  532894,  174792,  458201,  527381 },
		{  3000017, 1,  368334,  166765, 2871452,  407635,  979274, 1865572, 2703215 },
		{  4000037, 1,   72362,  210611,   92212,  583028,  681897, 2974319, 1680656 },
		{  6000011, 1, 1323844, 1723313, 1392620,  251332, 5750225,  908859, 5328166 }, 
		{  8000009, 1,   93973, 6914802, 3957321,  907968, 4380879, 3879127, 4791477 }, 
		{ 10000019, 1, 1833663, 3609180, 7252140, 5522715, 8914182, 4652083, 6262402 }
	};

	const int k12[][13] = {
        {  200003, 1,  24955, 142686,  75721, 189214, 164546, 183840,  58386,    775, 139837, 179994,  82896 },
        {  300017, 1,  78601, 167137, 290958, 193899,  91716, 160840,  68494, 191846, 133009, 248027,  65567 },
        {  500009, 1,  70666,  93673, 377076, 472997, 204370, 250473, 106427, 135013, 156929, 345112, 245626 },  
        { 1000003, 1, 361746, 775939, 987421, 525484, 164791, 306250, 380148, 605860, 774062, 992216,  92349 }
    };

    const int k10[][11] = {
        {  200003, 1, 43298,  88685, 25533, 111253,  160142,  124312, 180243, 44354, 10686 }
    };

    template <int n> const int* kor(int line);
    template <> const int* kor<8>(int line) { return k8[line]; }
    template <> const int* kor<12>(int line) { return k12[line]; }
    template <> const int* kor<10>(int line) { return k10[line]; }

	inline double frac(double x) {
		if (x > 0) return (x - (int)x);
		else if (x == 0) return 0;
		else return (x - (int)x - 1);
	} 

	class Point {
		public:
			Point(const int* line, const double* shift, int s) :
					line(line), shift(shift), s(s) {}
			double operator[](size_t i) const {
				return frac(shift[i] + 
						static_cast<double>(line[i+1]) / 
						line[0] * (s + 1));
			}

		private:
			friend class Iterator;
			const int* line;
			const double* shift;
			int s;
	};

	class Iterator {
		public:
			Iterator(const int* line, const double* shift, int s = 0) : 
					point(line, shift, s) {}

			void operator++() { ++point.s; }
			bool operator!=(const Iterator& other) {
				return point.s != other.point.s;
			}
		
			Point& operator*() { return point; }

		private:
			Point point;
	};

    template <int s>
	class Grid {
		public:
			Grid(int size = 0) { resize(size); update(); }

			void resize(int size) {
//      		std::cout << "korobov_size = " << size << ' ';
	        	for (line = 0; kor<s>(line)[0] < size; ++line) {}
                sz = kor<s>(line)[0];
//	        	std::cout << "sz = " << sz << ' ' << std::endl;
                update();
	        }

			int size() const { return sz; }
			typedef Iterator iterator;
			iterator begin() {
				return iterator(kor<s>(line), random_shift);
			}
			iterator end() {
				return iterator(kor<s>(line), random_shift, kor<s>(line)[0]);
			}
			void update() {
				for (int i = 0; i < s; ++i) 
					random_shift[i] = static_cast<double>(std::rand())/RAND_MAX;
			}
		private:
			int sz, line;
			double random_shift[s];
	};


}

#endif
