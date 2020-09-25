#ifndef DSGVG_HASH_HPP_INCLUDED
#define DSGVG_HASH_HPP_INCLUDED

#include <cstdint>

namespace dozyg {

uint32_t djb2_hash32(const char *str, uint64_t len);
uint64_t djb2_hash64(const char *str, uint64_t len);
uint64_t MurmurHash64A(const void * key, int len, uint64_t seed);

}

#endif
