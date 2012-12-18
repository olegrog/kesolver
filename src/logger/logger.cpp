#include <iomanip>
#include <ctime>
#include <sys/time.h>

#include "logger.hpp"

namespace logger {

const char* severity_names[] = {"[I]",
                                "[W]",
                                "[E]"};

std::ostringstream& Message::stream(Severity l)
{
    level = l;

    timeval tv;
    gettimeofday(&tv, NULL);

    tm* tm_p;
    tm_p = localtime(&tv.tv_sec);

    using namespace std;

    os << severity_names[static_cast<int>(level)] << ' ' 
       << setfill('0') 
       << setw(2) << tm_p->tm_hour << ':' 
       << setw(2) << tm_p->tm_min  << ':' 
       << setw(2) << tm_p->tm_sec  << '.' << tv.tv_usec 
       << " : ";

    return os;
}

Message::~Message()
{
    std::cerr << os.str() << std::endl;
}

} // namespace logger

