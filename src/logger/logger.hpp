#ifndef _LOGGER_HPP_
#define _LOGGER_HPP_

#include <iostream>
#include <sstream>

#define LOG(severity) \
    logger::Message().stream(logger::severity)

#ifdef DEBUG 
#define DLOG(severity) LOG(severity)
#else
#define DLOG(severity) \
    true ? (void)0 : logger::StreamVoidifier() & LOG(severity)
#endif

namespace logger {

enum Severity {
    INFO    = 0,
    WARNING = 1,
    ERROR   = 2
};

class Message {
    public:
        Message() {}

        ~Message();

        std::ostringstream& stream(Severity l);

    private:
        Severity level;
        std::ostringstream os;
        
        Message(const Message&);
        Message& operator=(const Message&);

};

class StreamVoidifier {
    public:
        StreamVoidifier() {}
        void operator& (std::ostream&) {}
};

} // namespace logger

#endif
