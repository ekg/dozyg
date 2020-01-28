#pragma once

#include <iostream>
#include "edlib.h"
#include "index.hpp"
#include "chain.hpp"

namespace gyeet {

struct alignment_t {
    std::string query_name;
    uint64_t query_length = 0;
    uint64_t query_begin = 0;
    uint64_t query_end = 0;
    uint64_t target_begin = 0;
    uint64_t target_end = 0;
    uint64_t anchor_count = 0;
    double chain_score = 0;
    bool is_secondary = false;
    double mapping_quality = 0;
    uint64_t edit_distance = 0;
    uint64_t score = 0;
    std::vector<handle_t> path;
    std::string cigar;
};


void for_handle_at_anchor_begin_in_chain(
    const chain_t& chain,
    const gyeet_index_t& index,
    const std::function<void(const handle_t&)>& func);

void write_chain_gaf(
    std::ostream& out,
    const chain_t& chain,
    const gyeet_index_t& index,
    const std::string& query_name,
    const uint64_t& query_length);

void write_alignment_gaf(
    std::ostream& out,
    const alignment_t& aln,
    const gyeet_index_t& index);

alignment_t align(
    const std::string& query_name,
    const uint64_t& query_total_length,
    const char* query,
    const chain_t& chain,
    const gyeet_index_t& index);

std::vector<handle_t> alignment_path(
    alignment_t& aln,
    const gyeet_index_t& index,
    const unsigned char* const alignment,
    const int alignmentLength);

std::string alignment_cigar(
    const unsigned char* const alignment,
    const int alignmentLength,
    const bool& extended_cigar);

}
