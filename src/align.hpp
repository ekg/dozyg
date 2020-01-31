#pragma once

#include <iostream>
#include "edlib.h"
#include "index.hpp"
#include "chain.hpp"

namespace gyeet {

typedef std::vector<std::pair<uint32_t, char>> cigar_t;
typedef std::vector<handle_t> path_t;

struct alignment_t {
    std::string query_name;
    // total query length
    uint64_t query_length = 0;
    // start and end on query
    uint64_t query_begin = 0;
    uint64_t query_end = 0;
    // relative to the path through the graph
    uint64_t target_begin = 0;
    uint64_t target_end = 0;
    // offsets relative to the handle we start and end in
    uint64_t first_handle_from_begin = 0;
    uint64_t last_handle_to_end = 0;
    uint64_t anchor_count = 0;
    double chain_score = 0;
    bool is_secondary = false;
    double mapping_quality = 0;
    uint64_t edit_distance = 0;
    int64_t score = 0;
    path_t path;
    cigar_t cigar;
    //std::string cigar;
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
    const gyeet_index_t& index,
    const uint64_t& extra_bp,
    const uint64_t& max_edit_distance);

bool has_matches(const cigar_t& cigar);
uint64_t insertion_length(const cigar_t& cigar);
int64_t score_cigar(
    const cigar_t& cigar,
    int64_t match = 1,
    int64_t mismatch = 2,
    int64_t gap_open = 2,
    int64_t gap_extend = 1);
std::ostream& operator<<(std::ostream& out, const cigar_t& cigar);
void extend_cigar_string(std::string& cigar_str, const cigar_t& cigar);
void extend_cigar(cigar_t& cigar, const cigar_t& extension);
void extend_cigar(cigar_t& cigar, const uint32_t& len, const char& type);

void graph_relativize(
    alignment_t& aln,
    seq_pos_t query_pos,
    seq_pos_t target_pos,
    const gyeet_index_t& index,
    const unsigned char* const alignment,
    const int alignmentLength,
    const bool extended_cigar);

uint64_t edit_distance_estimate(const chain_t& chain, const double& max_mismatch_rate);

alignment_t superalign(
    const std::string& query_name,
    const uint64_t& query_total_length,
    const char* query,
    const superchain_t& superchain,
    const gyeet_index_t& index,
    const uint64_t& extra_bp,
    const double& max_mismatch_rate);

std::string alignment_cigar(
    const unsigned char* const alignment,
    const int alignmentLength,
    const bool extended_cigar);

}
