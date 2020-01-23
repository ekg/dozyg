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
#include "dna.hpp"
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

struct node_ref_t {
    uint64_t idx = 0;  // index among edges
    uint64_t count_prev = 0;
};

class gyeet_index_t {
public:
    // total sequence length of the graph
    uint64_t seq_length = 0;
    // forward sequence of the graph, stored here for fast access during alignment
    std::string seq_fwd;
    // reverse complemented sequence of the graph, for fast access during alignment
    std::string seq_rev;
    // lets us map between our seq vector and handles (when our input graph is compacted!)
    sdsl::bit_vector::rank_1_type seq_bv_rank;
    // edge count
    uint64_t n_edges = 0;
    // what's the most-efficient graph topology we can store?
    mmappable_vector<handle_t> edges;
    // node count
    uint64_t n_nodes = 0;
    // refer to ranges in edges
    mmappable_vector<node_ref_t> node_ref;
    // number of kmers in the index
    uint64_t n_kmers = 0;
    // our kmer hash table
    boophf_t* bphf;
    // our kmer reference table (maps from bphf to index in kmer_pos_vec)
    mmappable_vector<uint64_t> kmer_pos_ref;
    // our kmer positions
    mmappable_vector<kmer_start_end_t> kmer_pos_table;
    gyeet_index_t(void) { }
    ~gyeet_index_t(void) {
        delete bphf;
    }
    void build(const HandleGraph& graph,
               const uint64_t& kmer_length,
               const uint64_t& max_furcations,
               const uint64_t& max_degree,
               const std::string& out_prefix);
    //void index_graph(const HandleGraph& graph);
    void load(const std::string& in_prefix);
    void for_values_of(const std::string& seq, const std::function<void(const kmer_start_end_t& v)>& lambda);
    // todo -- approximate graph alignment / GAF format
    //std::vector<handle_t> graph_alignment(...)
    // todo -- graph aware alignment scoring
    //uint64_t score_alignment(...)
};

}
