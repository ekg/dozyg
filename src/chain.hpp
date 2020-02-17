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
    // inner target boundaries
    seq_pos_t target_begin = 0;
    seq_pos_t target_end = 0;
    // query boundaries are fixed
    seq_pos_t query_begin(void) const { return anchors.front()->query_begin; }
    seq_pos_t query_end(void) const { return anchors.back()->query_end; }
    void compute_boundaries(const uint64_t& seed_length, const double& mismatch_rate) {
        // find the longest contiguous range covered by the chain
        // where the size is no more than some function of our seed_length * number of anchors * 1+mismatch_rate
        
        // find the inner boundaries
        auto& first_target_end = anchors.front()->target_end;
        auto& last_target_begin = anchors.back()->target_begin;
        // can we linearly extend to the start, end?
        auto& first_target_begin = anchors.front()->target_begin;
        auto& last_target_end = anchors.back()->target_end;
        //if (first_target_begin < last_target_end) {
        if (seq_pos::is_rev(first_target_begin) == seq_pos::is_rev(last_target_end) &&
            first_target_begin < last_target_end && score * (1+mismatch_rate) > last_target_end - first_target_begin) {
            target_begin = first_target_begin;
            target_end = last_target_end;
        } else if (seq_pos::is_rev(first_target_end) == seq_pos::is_rev(last_target_begin)
                   && target_begin < target_end) {
            target_begin = first_target_end;
            target_end = last_target_begin;
        } else {
            // kill chains with a single anchor that jumps nonlinearly
            /*
            std::cerr << "chain failure " << std::endl;
            for (auto& anchor : anchors) {
                std::cerr << &anchor << " " << "score=" << anchor->max_chain_score << ","
                          << "(" << seq_pos::to_string(anchor->query_begin) << ".."
                          << seq_pos::to_string(anchor->query_end) << "),"
                          << "(" << seq_pos::to_string(anchor->target_begin) << ".."
                          << seq_pos::to_string(anchor->target_end) << ")" << std::endl;
                    //<< "\"" << anchor.best_predecessor << "\" -> \"" << &anchor << "\";" << std::endl;
            }
            */
            score = -std::numeric_limits<double>::max();
        }
    }
};

struct chain_node_t {
    chain_t* chain = nullptr;
    chain_node_t* best_predecessor = nullptr;
    double max_superchain_score = 0;
    bool used = false;
    chain_node_t(chain_t* c) : chain(c) { }
};

std::vector<chain_t>
chains(std::vector<anchor_t>& anchors,
       const uint64_t& seed_length,
       const uint64_t& max_gap,
       const double& mismatch_rate,
       const uint64_t bandwidth = 50,
       const double secondary_chain_threshold = 0.5,
       const double max_mapq = 60);

double score_anchors(const anchor_t& a,
                     const anchor_t& b,
                     const uint64_t& seed_length,
                     const uint64_t& max_gap);

uint64_t chain_query_length(const chain_t& chain);

struct superchain_t {
    std::vector<chain_t*> chains;
    double score = 0;
    bool is_secondary = true;
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
            const uint64_t& kmer_length,
            const double& mismatch_rate,
            const uint64_t bandwidth = 1000);

double score_chain_nodes(const chain_node_t& a,
                         const chain_node_t& b,
                         const uint64_t& kmer_length);

}
