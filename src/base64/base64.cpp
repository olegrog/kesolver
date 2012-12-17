#include "base64.hpp"

namespace base64 {

char encoding_table[65] = 
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

unsigned char decoding_table[256];

struct decoding_table_maker {
    decoding_table_maker() {
        for (unsigned int i = 0; i < 256; ++i)
            decoding_table[i] = 0;
        for (unsigned char i = 0; i < 64; ++i)
            decoding_table[encoding_table[i]] = i;
    }
};

decoding_table_maker decoding_table_make;

std::vector<unsigned char> decode(const std::string& input)
{
    if (input.size() % 4 != 0) {
        // TODO: raise exception
    }

    std::vector<unsigned char> output;
    int output_length = input.size() / 4 * 3;
    if (input[input.size() - 1] == '=') --output_length;
    if (input[input.size() - 2] == '=') --output_length;

    output.resize(output_length);

    for (size_t i = 0, j = 0; i < input.size();) {

        unsigned int sextet_a = input[i] == '=' ? 0 & i++ : decoding_table[input[i++]];
        unsigned int sextet_b = input[i] == '=' ? 0 & i++ : decoding_table[input[i++]];
        unsigned int sextet_c = input[i] == '=' ? 0 & i++ : decoding_table[input[i++]];
        unsigned int sextet_d = input[i] == '=' ? 0 & i++ : decoding_table[input[i++]];

        unsigned int triple = (sextet_a << 3 * 6)
                            + (sextet_b << 2 * 6)
                            + (sextet_c << 1 * 6)
                            + (sextet_d << 0 * 6);

        if (j < output.size()) output[j++] = (triple >> 2 * 8) & 0xFF;
        if (j < output.size()) output[j++] = (triple >> 1 * 8) & 0xFF;
        if (j < output.size()) output[j++] = (triple >> 0 * 8) & 0xFF;
    }

    return output;
}

} // namespace base64
