#pragma once

#include <cstdio>
#include <string>
#include <fstream>

namespace dozyg {

void allocate_file(const std::string& fname, size_t size);
std::ifstream::pos_type filesize(const char* filename);

}
