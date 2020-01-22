#include <iostream>
#include <fstream>
#include <chrono>
#include <sdsl/bit_vectors.hpp>
#include "mmap_allocator.h"
#include "mmap_exception.h"
#include "mmappable_vector.h"
#include "ips4o.hpp"
#include "pmhf.hpp"
#include "algorithms/kmer.hpp"
#include "algorithms/hash.hpp"
#include "threads.hpp"
#include "utility.hpp"
//#include "mmap.hpp"

namespace gyeet {

using namespace mmap_allocator_namespace;

struct kmer_pos_t {
    uint64_t hash;
    int64_t begin;
    int64_t end;
};

struct kmer_start_end_t {
    int64_t begin;
    int64_t end;
};

class gyeet_index_t {
public:
    boophf_t* bphf;
    mmappable_vector<uint64_t, mmap_allocator<uint64_t>> kmer_pos_ptr;
    mmappable_vector<kmer_start_end_t, mmap_allocator<kmer_start_end_t>> kmer_pos_vec;
    gyeet_index_t(void) { }
    ~gyeet_index_t(void) {
        delete bphf;
    }
    void build(const HandleGraph& graph,
               const uint64_t& kmer_length,
               const uint64_t& max_furcations,
               const uint64_t& max_degree,
               const std::string& work_prefix);
    void save(const std::string& outfile);
    void load(const std::string& outfile);
    void for_values_of(const std::string& seq, const std::function<void(const kmer_start_end_t& v)>& lambda);
    //(kmer_start_end_t*&
    //void for_each_match(
};

}
