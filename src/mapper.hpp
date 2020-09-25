#pragma once

#include <vector>
#include <string>
#include <thread>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "atomic_queue.h"
#include "index.hpp"
#include "chain.hpp"
#include "align.hpp"

namespace dozyg {

struct seq_record_t {
    std::string name;
    std::string seq;
    seq_record_t(const std::string& n, const std::string& s)
        : name(n), seq(s) { }
};

// load into this
typedef atomic_queue::AtomicQueue<seq_record_t*, 2 << 16> seq_atomic_queue_t;
// results into this, write out
typedef atomic_queue::AtomicQueue<std::string*, 2 << 16> gaf_atomic_queue_t;

void for_each_seq_in_file(
    const std::string& filename,
    const std::function<void(const std::string&, const std::string&)>& func);

void reader_thread(
    const std::vector<std::string>& files,
    seq_atomic_queue_t& seq_queue,
    std::atomic<bool>& reader_done);

void writer_thread(
    std::ostream& out,
    gaf_atomic_queue_t& gaf_queue,
    const std::vector<std::atomic<bool>>& working);

bool still_working(
    const std::vector<std::atomic<bool>>& working);

std::string map_seq(
    const std::string& name,
    const std::string& query,
    const dozyg_index_t& index,
    const uint64_t& max_chain_gap,
    const double& mismatch_rate,
    const uint64_t& chain_min_n_anchors,
    const double& chain_overlap_max,
    const uint64_t& align_best_n,
    const bool& write_alignments,
    const bool& write_chains,
    const bool& write_superchains);

void worker_thread(
    uint64_t tid,
    std::atomic<bool>& is_working,
    seq_atomic_queue_t& seq_queue,
    gaf_atomic_queue_t& gaf_queue,
    const std::atomic<bool>& reader_done,
    const dozyg_index_t& index,
    uint64_t max_chain_gap,
    double mismatch_rate,
    uint64_t chain_min_n_anchors,
    double chain_overlap_max,
    uint64_t align_best_n,
    bool write_alignment,
    bool write_chains,
    bool write_superchains);

void map_reads(
    const std::vector<std::string>& input_files,
    const dozyg_index_t& index,
    const uint64_t& max_chain_gap,
    const double& mismatch_rate,
    const uint64_t& chain_min_n_anchors,
    const double& chain_overlap_max,
    const uint64_t& align_best_n,
    const uint64_t& nthreads,
    const bool& write_alignment,
    const bool& write_chains,
    const bool& write_superchains);


}
