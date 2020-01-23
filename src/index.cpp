#include "index.hpp"

namespace gyeet {

void gyeet_index_t::build(const HandleGraph& graph,
                          const uint64_t& kmer_length,
                          const uint64_t& max_furcations,
                          const uint64_t& max_degree,
                          const std::string& out_prefix) {

    // stash our node count (helps when loading our index files
    n_nodes = graph.get_node_count();
    // collect our sequences
    graph.for_each_handle(
        [&](const handle_t& h) {
            seq_length += graph.get_length(h);
        });
    seq_fwd.reserve(seq_length);
    graph.for_each_handle(
        [&](const handle_t& h) {
            seq_fwd.append(graph.get_sequence(h));
        });
    seq_rev = dna::reverse_complement(seq_fwd);
    // mark our node starts in fwd
    sdsl::bit_vector seq_bv(seq_length);
    uint64_t seq_idx = 0;
    graph.for_each_handle(
        [&](const handle_t& h) {
            seq_bv[seq_idx] = 1;
            seq_idx += graph.get_length(h);
        });
    sdsl::util::assign(seq_bv_rank, sdsl::bit_vector::rank_1_type(&seq_bv));
    // save our sequences
    std::string seq_fwd_idx = out_prefix + ".sqf";
    std::ofstream seq_fwd_f(seq_fwd_idx.c_str());
    seq_fwd_f << seq_fwd;
    seq_fwd_f.close();
    std::string seq_rev_idx = out_prefix + ".sqr";
    std::ofstream seq_rev_f(seq_rev_idx.c_str());
    seq_rev_f << seq_rev;
    seq_rev_f.close();
    // and our rank structure
    std::string seq_bv_rank_idx = out_prefix + ".sbv";
    std::ofstream seq_bv_rank_f(seq_bv_rank_idx.c_str());
    seq_bv_rank.serialize(seq_bv_rank_f);
    seq_bv_rank_f.close();

    // follow forward and reverse edges
    std::string edge_idx = out_prefix + ".gye";
    std::ofstream edge_f(edge_idx.c_str());
    std::string node_ref_idx = out_prefix + ".gyn";
    std::ofstream node_ref_f(node_ref_idx.c_str());
    uint64_t ref_idx = 0;
    graph.for_each_handle(
        [&](const handle_t& h) {
            node_ref_t ref = { ref_idx++, 0 };
            graph.follow_edges(
                h, true,
                [&](const handle_t& p) {
                    edge_f.write((char*)&p, sizeof(p));
                    ++ref.count_prev;
                });
            n_edges += ref.count_prev;
            node_ref_f.write((char*)&ref, sizeof(ref));
            graph.follow_edges(
                h, false,
                [&](const handle_t& n) {
                    edge_f.write((char*)&n, sizeof(n));
                    ++n_edges;
                });
        });
    {
        // write a marker reference, to simplify counting of edges
        node_ref_t ref = { ref_idx, 0 };
        node_ref_f.write((char*)&ref, sizeof(ref));
    }

    edge_f.close();
    node_ref_f.close();

    edges.mmap_file(edge_idx.c_str(), READ_WRITE_SHARED, 0, n_edges);
    node_ref.mmap_file(node_ref_idx.c_str(), READ_WRITE_SHARED, 0, n_nodes+1);

    // map from kmer hash idx to kmer start/end in graph sequence vector
    std::string kmer_pos_idx = out_prefix + ".kpos";
    std::ofstream kmer_pos_f(kmer_pos_idx.c_str());

    std::string kmer_set_idx = out_prefix + ".kset";
    std::ofstream kmer_set_f(kmer_set_idx.c_str());
        
    uint64_t n_kmers = 0;
    algorithms::for_each_kmer(
        graph, kmer_length, max_furcations, max_degree,
        [&](const kmer_t& kmer) {
            uint64_t hash = djb2_hash64(kmer.seq.c_str());
            auto& b = kmer.begin;
            int64_t start_pos = seq_bv_rank(id(b)) + (is_rev(b) ? graph.get_length(graph.get_handle(id(b))) - offset(b) : offset(b));
            if (is_rev(b)) start_pos = -start_pos;
            auto& e = kmer.end;
            int64_t end_pos = seq_bv_rank(id(e)) + (is_rev(e) ? graph.get_length(graph.get_handle(id(e))) - offset(e) : offset(e));
            if (is_rev(e)) end_pos = -end_pos;
            auto v = std::make_pair(start_pos, end_pos);
            kmer_pos_t p = { hash, start_pos, end_pos };
            //kmer_map.add(bbidx, std::make_pair(start_pos, end_pos));
#pragma omp critical (kmer_pos_write)
            {
                ++n_kmers;
                kmer_pos_f.write((char*)&p, sizeof(kmer_pos_t));
            }
#pragma omp critical (kmer_set_write)
            {
                kmer_set_f.write((char*)&hash, sizeof(uint64_t));
            }
        });
    kmer_pos_f.close();
    kmer_set_f.close();
        
    //ska::flat_hash_map<uint32_t, uint32_t> kmer_table;
    //std::vector<uint64_t> kmers;

    //std::cerr << std::endl;
    mmappable_vector<uint64_t> kmer_set;
    kmer_set.mmap_file(kmer_set_idx.c_str(), READ_WRITE_SHARED, 0, n_kmers);
    //mmap_buffer_t kmer_set_buf = open_mmap_buffer(kmer_set_idx.c_str());
    ips4o::parallel::sort(kmer_set.begin(), kmer_set.end());
    //std::cerr << "total kmers " << kmer_set.size() << std::endl;
    kmer_set.erase(std::unique(kmer_set.begin(), kmer_set.end()), kmer_set.end());
    //std::cerr << "unique kmers " << kmer_set.size() << std::endl;

    double gammaFactor = 8.0;
    bphf = new boomphf::mphf<uint64_t,hasher_t>(kmer_set.size(),
                                                kmer_set,
                                                get_thread_count(),
                                                gammaFactor,
                                                false,false);

    // remove the kmer set working file
    kmer_set.munmap_file();
    std::remove(kmer_set_idx.c_str());

    // rename our kmers
    // build our sequence space index
    mmappable_vector<kmer_pos_t> kmer_pos;
    kmer_pos.mmap_file(kmer_pos_idx.c_str(), READ_WRITE_SHARED, 0, n_kmers);
    // set our kmer hashes
#pragma omp parallel for
    for (uint64_t i = 0; i < n_kmers; ++i) {
        auto& kp = kmer_pos[i];
        kp.hash = bphf->lookup(kp.hash);
    }

    ips4o::parallel::sort(
        kmer_pos.begin(), kmer_pos.end(),
        [](const kmer_pos_t& a, const kmer_pos_t& b) {
            return a.hash < b.hash
                            || a.hash == b.hash
                            && (a.begin < b.begin
                                || a.begin == b.begin && a.end < b.end);
        });

    // now we can iterate through our keys/values in order, storing them somewhere
    std::string kmer_pos_table_idx = out_prefix + ".kpv";
    std::ofstream kmer_pos_table_f(kmer_pos_table_idx.c_str());
    std::string kmer_pos_ref_idx = out_prefix + ".kpp";
    std::ofstream kmer_pos_ref_f(kmer_pos_ref_idx.c_str());

    uint64_t last_hash = std::numeric_limits<uint64_t>::max();
    uint64_t marker_idx = 0;
    for (uint64_t i = 0; i < n_kmers; ++i) {
        auto& kp = kmer_pos[i];
        if (kp.hash != last_hash) {
            kmer_pos_ref_f.write((char*)&marker_idx, sizeof(uint64_t));
        }
        kmer_start_end_t p = { kp.begin, kp.end };
        kmer_pos_table_f.write((char*)&p, sizeof(kmer_start_end_t));
        ++marker_idx;
    }
    kmer_pos_table_f.close();
    // write the total in the last entry
    kmer_pos_ref_f.write((char*)&marker_idx, sizeof(uint64_t));
    kmer_pos_ref_f.close();

    // remove the temporary file of kmer positions
    kmer_pos.munmap_file();
    std::remove(kmer_pos_idx.c_str());

    //std::cerr << "querying kmers" << std::endl;
    /*
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for (auto& hash : kmer_set) {
        uint64_t k = bphf->lookup(hash);
        uint64_t c = kmer_pos_ref[k+1] - kmer_pos_ref[k];
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto used_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    std::cerr << "done with " << kmer_set.size() << " @ " << (double)used_time/(double)kmer_set.size() << "ns/kmer" << std::endl;
    */

    // write the bbhash index
    std::string bbhash_out_idx = out_prefix + ".bbx";
    std::ofstream f(bbhash_out_idx.c_str());
    bphf->save(f);

    /*
      algorithms::for_each_kmer(graph, args::get(kmer_length), [&](const kmer_t& kmer) {
      uint64_t hash = djb2_hash64(kmer.seq.c_str());
      bphf->lookup(hash);
      #pragma omp atomic
      ++seen_kmers;
      });
      std::cerr << "done with " << seen_kmers << " kmers" << std::endl;
    */

}

void gyeet_index_t::load(const std::string& in_prefix) {
    bphf = new boomphf::mphf<uint64_t,hasher_t>();
    std::string bbx_in_idx = in_prefix + ".bbx";
    std::ifstream f(bbx_in_idx.c_str());
    bphf->load(f);
    f.close();
    std::string kmer_pos_table_idx = in_prefix + ".kpv";
    n_kmers = filesize(kmer_pos_table_idx.c_str());
    kmer_pos_table.mmap_file(kmer_pos_table_idx.c_str(), READ_WRITE_SHARED, 0, n_kmers);
    std::string kmer_pos_ref_idx = in_prefix + ".kpp";
    kmer_pos_ref.mmap_file(kmer_pos_ref_idx.c_str(), READ_WRITE_SHARED, 0, n_kmers+1);
}

void gyeet_index_t::for_values_of(const std::string& seq, const std::function<void(const kmer_start_end_t& v)>& lambda) {
    uint64_t hash = djb2_hash64(seq.c_str());
    uint64_t k = bphf->lookup(hash);
    uint64_t b = kmer_pos_ref[k];
    uint64_t e = kmer_pos_ref[k+1];
    for (uint64_t i = b; i < e; ++i) {
        lambda(kmer_pos_table[i]);
    }
}

}
