#pragma once

#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <IITree.h> // cgranges
#include "index.hpp"

namespace gyeet {

struct anchor_t {
    seq_pos_t query_begin = 0;
    seq_pos_t query_end = 0;
    seq_pos_t target_begin = 0;
    seq_pos_t target_end = 0;
    double max_chain_score = 0;
    anchor_t* best_predecessor = nullptr;
    anchor_t(const seq_pos_t& qb,
             const seq_pos_t& qe,
             const seq_pos_t& tb,
             const seq_pos_t& te)
        : query_begin(qb)
        , query_end(qe)
        , target_begin(tb)
        , target_end(te) { }
};

std::vector<anchor_t> anchors_for_query(const gyeet_index_t& index,
                                        const char* seq,
                                        const size_t& len);

struct chain_t {
    std::vector<anchor_t*> anchors;
    double score = 0;
    double mapping_quality = std::numeric_limits<double>::min();
    bool is_secondary = false;
    bool processed(void) {
        return mapping_quality != std::numeric_limits<double>::min();
    }
};

struct chain_node_t {
    chain_t* chain = nullptr;
    chain_node_t* best_predecessor = nullptr;
    double max_superchain_score = 0;
    chain_node_t(chain_t* c) : chain(c) { }
};

std::vector<chain_t>
chains(std::vector<anchor_t>& anchors,
       const uint64_t& seed_length,
       const uint64_t& max_gap,
       const uint64_t bandwidth = 50,
       const double secondary_chain_threshold = 0.5,
       const double max_mapq = 60);

double score_anchors(const anchor_t& a,
                     const anchor_t& b,
                     const uint64_t& seed_length,
                     const uint64_t& max_gap);

struct superchain_t {
    std::vector<chain_t*> chains;
    double score = 0;
    /*
    double mapping_quality = std::numeric_limits<double>::min();
    //bool is_secondary = false;
    bool processed(void) {
        return mapping_quality != std::numeric_limits<double>::min();
    }
    */
};

std::vector<superchain_t>
superchains(std::vector<chain_t>& chains,
            const uint64_t bandwidth = 50);

double score_chain_nodes(const chain_node_t& a,
                         const chain_node_t& b);

}
