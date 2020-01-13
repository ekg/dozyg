#pragma once

#include <cstdint> 
#include "BooPHF.h"

namespace gyeet {

typedef boomphf::SingleHashFunctor<uint64_t>  hasher_t;
typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

}
