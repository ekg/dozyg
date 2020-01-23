#include "utility.hpp"

namespace gyeet {

void allocate_file(const std::string& fname, size_t size) {
    FILE* f = fopen(fname.c_str(), "w");
    int r = fseek(f, size-1, SEEK_SET);
    char x = '\0';
    size_t w = fwrite(&x, 1, 1, f);
    r = fclose(f);
}

std::ifstream::pos_type filesize(const char* filename) {
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

}
