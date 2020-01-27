#pragma once

#include <iostream>
#include "edlib.h"
#include "index.hpp"
#include "chain.hpp"

namespace gyeet {

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

struct alignment_t {
    uint64_t query_begin = 0;
    uint64_t query_end = 0;
    uint64_t target_begin = 0;
    uint64_t target_end = 0;
    uint64_t score = 0;
    double mapping_quality = 0;
    std::vector<handle_t> path;
    std::string cigar;
};

alignment_t align(
    const char* query,
    const chain_t& chain,
    const gyeet_index_t& index);

void fill_alignment_path(
    alignment_t& aln,
    const gyeet_index_t& index,
    const unsigned char* const alignment,
    const int alignmentLength);

void fill_alignment_cigar(
    alignment_t& aln,
    const unsigned char* const alignment,
    const int alignmentLength,
    const bool& extended_cigar);


}
