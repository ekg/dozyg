#include "hash.hpp"

namespace gyeet {

uint32_t djb2_hash32(const char *str, uint64_t len) {
    uint32_t hash = 5381;
    while (len-- > 0) {
        hash = ((hash << 5) + hash) + *str++; /* hash * 33 + c */
    }
    return hash;
}

uint64_t djb2_hash64(const char *str, uint64_t len) {
    uint64_t hash = 5381;
    while (len-- > 0) {
        hash = ((hash << 5) + hash) + *str++; /* hash * 33 + c */
    }
    return hash;
}

}
