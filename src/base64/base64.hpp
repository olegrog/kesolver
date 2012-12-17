#ifndef _BASE64_HPP_
#define _BASE64_HPP_

#include <vector>
#include <string>

namespace base64 {

std::vector<unsigned char> decode(const std::string& input);

}

#endif // _BASE64_HPP_
