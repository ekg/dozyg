#pragma once

#include <cstdint> 
#include "BooPHF.h"

namespace dozyg {

template <typename T> class HashHasher {
public:
	// the class should have operator () with this signature :
	// BBhash will use the 'seed' paramater to generate two different hash values form this key.
     //then it will generate internally a sequence of hash values using xorshifts, using these two first hash values as starting point.
	uint64_t operator () (T key, T seed=0) const {
		return key ^= seed;
	}
};

typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
//typedef HashHasher<uint64_t> hasher_t;
typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

}
