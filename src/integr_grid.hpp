#ifndef _INTEGR_GRID_H_
#define _INTEGR_GRID_H_

#include <cstdlib>
#include <climits>
#include <limits>

namespace korobov {

const int dimension = 8;
const int coefficients[][1 + dimension] = {
    {    10007, 1,    1905,    6491,    6710,    3611,    4146,    2607,    2863 },
    {    20011, 1,    7191,    2057,    3758,    8928,    5960,   14809,   12988 },
    {    30011, 1,   13167,   26353,    2769,   26069,   14716,   14556,    8606 },
    {    50021, 1,   11281,    7537,   39218,   32534,   11977,    5816,   32765 }, 
    {   100003, 1,   20285,   68883,   49739,   25348,   68757,   93907,   46351 },
    {   200003, 1,   47369,  188507,   54145,  156036,  158419,   37051,   42494 },
    {   300017, 1,   81575,  103565,  136172,  101475,   54078,  262899,  170731 },
    {   400009, 1,  193141,  206577,  390670,  296791,   20804,   14959,  331221 },
    {   500009, 1,   42535,  193663,  307439,  182488,  487373,   37415,  418387 },
    {   750019, 1,   10525,  522832,  667416,  625465,  102362,  332766,  523439 },
    {  1000003, 1,  417564,  171019,  163483,  410620,  615303,  611111,  188073 },
    {  1500007, 1,  413996,  388189,  943278, 1496508,  434758,  733031,  985685 },
    {  2000003, 1,  832685, 1269182, 1228431,  532894,  174792,  458201,  527381 },
    {  3000017, 1,  368334,  166765, 2871452,  407635,  979274, 1865572, 2703215 },
    {  4000037, 1,   72362,  210611,   92212,  583028,  681897, 2974319, 1680656 },
    {  6000011, 1, 1323844, 1723313, 1392620,  251332, 5750225,  908859, 5328166 }, 
    {  8000009, 1,   93973, 6914802, 3957321,  907968, 4380879, 3879127, 4791477 }, 
    { 10000019, 1, 1833663, 3609180, 7252140, 5522715, 8914182, 4652083, 6262402 },
    { 1000000000, 1, 1, 1, 1, 1, 1, 1, 1}
};

inline double frac(double x) {
    if (x > 0) return (x - (int)x);
    else if (x == 0) return 0;
    else return (x - (int)x + 1);
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

class Grid {
    public:
        Grid(int size = 0);

        int size() const { return size_; }

        typedef Iterator iterator;
        typedef Point point;

        iterator begin() {
            return iterator(coefficients[line], random_shift);
        }
        iterator end() {
            return iterator(coefficients[line], random_shift, coefficients[line][0]);
        }

    private:
        int size_, line;
        double random_shift[dimension];
};

inline Grid::Grid(int size) {
    for (line = 0; coefficients[line][0] < size; ++line) {}
    size_ = coefficients[line][0];
    for (int i = 0; i < dimension; ++i) 
        random_shift[i] = static_cast<double>(std::rand())/RAND_MAX;
}

}

namespace crand {

const int dimension = 10;

class Point {
    public:
        Point() {
            for (int i = 0; i < dimension; ++i) 
                random[i] = static_cast<double>(std::rand())/RAND_MAX;
        }
        double operator[](size_t i) const {
            return random[i];
        }

    private:
        double random[dimension];
};

class Iterator {
    public:
        Iterator(int s = 0) : s(s) {}

        void   operator++ () { point = Point(); ++s; }
        bool   operator!= (const Iterator& other) { return s != other.s; }
        Point& operator*  () { return point; }

    private:
        Point point;
        int s;
};

class Grid {
    public:
        Grid(int size = 0) : size_(size) {}

        int size() const { return size_; }

        typedef Iterator iterator;
        typedef Point point;

        iterator begin() {
            return iterator(0);
        }
        iterator end() {
            return iterator(size_);
        }

    private:
        int size_;
};

}

namespace sobol {

unsigned int v[][32] = {
    {0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x08000000, 0x04000000, 0x02000000, 0x01000000,
     0x00800000, 0x00400000, 0x00200000, 0x00100000, 0x00080000, 0x00040000, 0x00020000, 0x00010000,
     0x00008000, 0x00004000, 0x00002000, 0x00001000, 0x00000800, 0x00000400, 0x00000200, 0x00000100,
     0x00000080, 0x00000040, 0x00000020, 0x00000010, 0x00000008, 0x00000004, 0x00000002, 0x00000001},
    {0x80000000, 0xC0000000, 0xA0000000, 0xF0000000, 0x88000000, 0xCC000000, 0xAA000000, 0xFF000000,
     0x80800000, 0xC0C00000, 0xA0A00000, 0xF0F00000, 0x88880000, 0xCCCC0000, 0xAAAA0000, 0xFFFF0000,
     0x80008000, 0xC000C000, 0xA000A000, 0xF000F000, 0x88008800, 0xCC00CC00, 0xAA00AA00, 0xFF00FF00,
     0x80808080, 0xC0C0C0C0, 0xA0A0A0A0, 0xF0F0F0F0, 0x88888888, 0xCCCCCCCC, 0xAAAAAAAA, 0xFFFFFFFF},
    {0x80000000, 0x40000000, 0xE0000000, 0xB0000000, 0x68000000, 0xF4000000, 0x86000000, 0x4F000000,
     0xE8800000, 0xB4400000, 0x66E00000, 0xFFB00000, 0x80E80000, 0x40B40000, 0xE0660000, 0xB0FF0000,
     0x68808000, 0xF4404000, 0x86E0E000, 0x4FB0B000, 0xE8E86800, 0xB4B4F400, 0x66668600, 0xFFFF4F00,
     0x80006880, 0x4000F440, 0xE00086E0, 0xB0004FB0, 0x6800E8E8, 0xF400B4B4, 0x86006666, 0x4F00FFFF},
    {0x80000000, 0xC0000000, 0xE0000000, 0x50000000, 0x38000000, 0xAC000000, 0x62000000, 0x93000000,
     0xDB800000, 0xFD400000, 0x5AE00000, 0x3DB00000, 0xB8080000, 0x6C0C0000, 0x820E0000, 0xC3050000,
     0xE3838000, 0x514AC000, 0x38E62000, 0xAEB93000, 0x6385B800, 0x9143D400, 0xD8EBAE00, 0xFEB6DB00,
     0x5B800080, 0x3D4000C0, 0xBAE000E0, 0x6DB00050, 0x80080038, 0xC00C00AC, 0xE00E0062, 0x50050093},
    {0x80000000, 0x40000000, 0xA0000000, 0x30000000, 0x78000000, 0xCC000000, 0xFA000000, 0x8D000000,
     0x58800000, 0xBDC00000, 0x21600000, 0x72F00000, 0xD8880000, 0xFDC40000, 0x816A0000, 0x42F30000,
     0xA08F8000, 0x31C8C000, 0x7B65A000, 0xCFFBD000, 0xF80A0800, 0x8C031C00, 0x5A07B600, 0xBD0CFF00,
     0x208F8080, 0x71C8C040, 0xDB65A0A0, 0xFFFBD030, 0x800A0878, 0x40031CCC, 0xA007B6FA, 0x300CFF8D},
    {0x80000000, 0xC0000000, 0x20000000, 0x10000000, 0x48000000, 0xEC000000, 0x32000000, 0x59000000,
     0xA0800000, 0xD0C00000, 0x68200000, 0xFC100000, 0x7A480000, 0xB5EC0000, 0x92B20000, 0x89990000,
     0xC8008000, 0x2C00C000, 0x12002000, 0x49001000, 0xE8804800, 0x3CC0EC00, 0x5A203200, 0xA5105900,
     0xDAC8A080, 0x652CD0C0, 0xFA926820, 0x7589FC10, 0xB248FA48, 0x99EC75EC, 0x80B2B2B2, 0xC0999999},
    {0x80000000, 0x40000000, 0x60000000, 0x70000000, 0xF8000000, 0xBC000000, 0xDA000000, 0xAD000000,
     0x5A800000, 0xED400000, 0x3AE00000, 0x9D300000, 0xC2180000, 0x218C0000, 0x18C20000, 0x8C210000,
     0x42188000, 0x618C4000, 0x78C26000, 0xFC217000, 0xBA187800, 0xDD8CFC00, 0xA2C2BA00, 0x5121DD00,
     0xE0982280, 0x30CC1140, 0x982280E0, 0xCC114030, 0x2280E098, 0x114030CC, 0x80E09822, 0x4030CC11},
    {0x80000000, 0xC0000000, 0x60000000, 0x90000000, 0x48000000, 0xE4000000, 0x56000000, 0x2B000000,
     0x70800000, 0x1C400000, 0xC8200000, 0x24300000, 0x36180000, 0xBB240000, 0x38920000, 0xF8790000,
     0x9E358000, 0x0F3AC000, 0x46842000, 0xA7631000, 0xF0800800, 0xDC400C00, 0xA8200600, 0xB4300900,
     0x7E180480, 0x5F240E40, 0x6E920560, 0xD37902B0, 0xEEB58708, 0x137AC1C4, 0x8EA42C82, 0x83531243},
    {0x80000000, 0xC0000000, 0xE0000000, 0x70000000, 0xA8000000, 0xF4000000, 0x6E000000, 0x13000000,
     0x1D800000, 0xBE400000, 0xEE200000, 0xD3300000, 0xFDB80000, 0xCE5C0000, 0x460A0000, 0x270D0000,
     0x93A38000, 0xDD58C000, 0x5B8D6000, 0x99629000, 0x7DB80800, 0x0E5C0C00, 0xA60A0E00, 0x570D0700,
     0x3BA38A80, 0x2958CF40, 0x358D66E0, 0x8A629130, 0x603809D8, 0xB01C07E4, 0x482A00E2, 0x843D0A33},
    {0x80000000, 0x40000000, 0xA0000000, 0xB0000000, 0xD8000000, 0xD4000000, 0x8A000000, 0x19000000,
     0x33800000, 0x99C00000, 0x72200000, 0x3D100000, 0xC1A80000, 0xE4EC0000, 0x13960000, 0x69E50000,
     0x0A2A8000, 0x593A4000, 0x9392E000, 0x29FF7000, 0xAA280800, 0xE92C0400, 0x4BB60A00, 0xFDF50B00,
     0x20028D80, 0xF0164D40, 0x7824E8A0, 0x640A7190, 0x522A8B38, 0xCD3A4D9C, 0xB992ED22, 0x80FF78D1},
    {0x80000000, 0x40000000, 0xE0000000, 0x30000000, 0xE8000000, 0xCC000000, 0x5E000000, 0x61000000,
     0x74800000, 0x09C00000, 0xFCA00000, 0xB5D00000, 0xAA980000, 0x28DC0000, 0x68220000, 0x8C2F0000,
     0xBE158000, 0x51274000, 0x9CB0A000, 0xC5E93000, 0xA2958800, 0xD4E74400, 0xDE10AE00, 0x21393300,
     0x948D8680, 0x39FB48C0, 0x1492ABE0, 0x79C63510, 0xF48001C8, 0x49C0085C, 0x1CA0042A, 0x85D00E4D},
    {0x80000000, 0xC0000000, 0xE0000000, 0xD0000000, 0x18000000, 0x8C000000, 0xB2000000, 0x09000000,
     0x75800000, 0xE8400000, 0xA7A00000, 0xF1700000, 0x2A180000, 0x45040000, 0x27BE0000, 0x31570000,
     0xCA0A8000, 0x95114000, 0x3F89E000, 0xBD4C5000, 0x78128800, 0x9C154C00, 0x4A37EE00, 0x551B5D00,
     0xDF980980, 0x6D4404C0, 0x601E0520, 0x10270D90, 0xF81286D8, 0x5C154644, 0xAA37E15A, 0x851B5F87},
    {0x80000000, 0xC0000000, 0xA0000000, 0x10000000, 0x78000000, 0x4C000000, 0xE2000000, 0x73000000,
     0xCD800000, 0x27400000, 0xD7A00000, 0xD8700000, 0x58080000, 0x9C340000, 0x3A360000, 0x2F170000,
     0x57A68000, 0x184FC000, 0xF82BE000, 0x8C251000, 0x422E8800, 0x633BCC00, 0xB5BDEA00, 0x6B421100,
     0x35800F80, 0xAB4008C0, 0x95A00420, 0xBB700630, 0xED880B58, 0xF77406B4, 0x0F96035A, 0x84670AB7}
};

template <typename T, int n>
class Array {
    public:
        T&      operator[](int i)       { return a[i]; }
        const T operator[](int i) const { return a[i]; }
    private:
        T a[n];
};

template <int n>
Array<double, n> sobol(unsigned int i)
{
    unsigned int z[n];
    for (int j = 0; j < n; ++j) 
        z[j] = 0;
    int j = 0;
    while (i > 0) {
        int x = i % 2;
        if (x > 0) 
            for (int k = 0; k < n; ++k) 
                z[k] ^= v[k][j];
        i /= 2;
        ++j;
    }
    Array<double, n> y;
    for (int k = 0; k < n; ++k) 
        y[k] = static_cast<double>(z[k]) / std::numeric_limits<unsigned int>::max();
    return y;
}

double frac(double x) {
    return x - floor(x);
}


template <int n>
class Iterator {
    public:
        Iterator(unsigned int i, const double* s) : i(i) {
            update();
            for (int j = 0; j < n; ++j) shift[j] = s[j];
        }
        Iterator& operator++() { 
            ++i;
            update();
            return *this;
        }
        bool operator!=(const Iterator& other) { return i != other.i; }
        const Array<double, n>& operator*() const { return point; }
    private:
        unsigned int i;
        double shift[n];
        Array<double, n> point;

        void update() { 
            Array<double, n> p = sobol<n>(i);
            for (int j = 0; j < n; ++j) 
                point[j] = frac(p[j] + shift[j]);
            
        }
};

template <int n>
class Grid {
    public:
        Grid(int size = 0) : size_(size) {
            for (int i = 0; i < n; ++i) 
                shift[i] = static_cast<double>(std::rand())/RAND_MAX;
            start = std::rand();
        }

        int size() const { return size_; }

        typedef Iterator<n> iterator;
        typedef Array<double, n>    point;

        iterator begin() { return iterator(start, shift); }
        iterator end()   { return iterator(size_+start, shift); }

    private:
        double shift[n];
        size_t size_;
        unsigned int start;
};

}

#endif
