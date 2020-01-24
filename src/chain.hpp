#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include "index.hpp"

namespace gyeet {

struct anchor_t {
    seq_pos_t query_begin = 0;
    seq_pos_t query_end = 0;
    seq_pos_t target_begin = 0;
    seq_pos_t target_end = 0;
    double max_chain_score = 0;
};

struct chain_t {
    std::vector<anchor_t*> anchors;
    double score = 0;
    bool is_secondary = false;
};

std::vector<chain_t>
chains(std::vector<anchor_t>& anchors,
       const uint64_t& seed_length,
       const uint64_t& max_gap,
       const uint64_t bandwidth = 50);

double score_anchors(const anchor_t& a,
                     const anchor_t& b,
                     const uint64_t& seed_length,
                     const uint64_t& max_gap);

}
