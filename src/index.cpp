#include "index.hpp"

namespace gyeet {

void gyeet_index_t::build(const HandleGraph& graph,
                          const uint64_t& kmer_length,
                          const uint64_t& max_furcations,
                          const uint64_t& max_degree,
                          const std::string& work_prefix) {

    uint64_t total_seq_length = 0;
    graph.for_each_handle(
        [&](const handle_t& h) {
            total_seq_length += graph.get_length(h);
        });
    sdsl::bit_vector seq_bv(total_seq_length);
    uint64_t seq_idx = 0;
    graph.for_each_handle(
        [&](const handle_t& h) {
            seq_bv[seq_idx] = 1;
            seq_idx += graph.get_length(h);
        });
    sdsl::bit_vector::rank_1_type seq_bv_rank;
    sdsl::util::assign(seq_bv_rank, sdsl::bit_vector::rank_1_type(&seq_bv));

    // map from kmer hash idx to kmer start/end in graph sequence vector
    std::string kmer_pos_idx = work_prefix + ".kmer_pos";
    std::ofstream kmer_pos_f(kmer_pos_idx.c_str());

    std::string kmer_set_idx = work_prefix + ".kmer_set";
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
    mmappable_vector<uint64_t, mmap_allocator<uint64_t>> kmer_set;
    kmer_set.mmap_file(kmer_set_idx.c_str(), READ_WRITE_SHARED, 0, n_kmers);
    //mmap_buffer_t kmer_set_buf = open_mmap_buffer(kmer_set_idx.c_str());
    ips4o::parallel::sort(kmer_set.begin(), kmer_set.end());
    std::cerr << "total kmers " << kmer_set.size() << std::endl;
    kmer_set.erase(std::unique(kmer_set.begin(), kmer_set.end()), kmer_set.end());
    std::cerr << "unique kmers " << kmer_set.size() << std::endl;

    double gammaFactor = 8.0;
    bphf = new boomphf::mphf<uint64_t,hasher_t>(kmer_set.size(),
                                                kmer_set,
                                                get_thread_count(),
                                                gammaFactor,
                                                false,false);

    //kmers.clear();
    // rename our kmers

    // build our sequence space index
    mmappable_vector<kmer_pos_t, mmap_allocator<kmer_pos_t>> kmer_pos;
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
    std::string kmer_pos_vec_idx = work_prefix + ".kmer_pos_vec";
    std::ofstream kmer_pos_vec_f(kmer_pos_vec_idx.c_str());
    std::string kmer_pos_ptr_idx = work_prefix + ".kmer_pos_ptr";
    std::ofstream kmer_pos_ptr_f(kmer_pos_ptr_idx.c_str());

    uint64_t last_hash = std::numeric_limits<uint64_t>::max();
    uint64_t marker_idx = 0;
    for (uint64_t i = 0; i < n_kmers; ++i) {
        auto& kp = kmer_pos[i];
        if (kp.hash != last_hash) {
            kmer_pos_ptr_f.write((char*)&marker_idx, sizeof(uint64_t));
        }
        kmer_start_end_t p = { kp.begin, kp.end };
        kmer_pos_vec_f.write((char*)&p, sizeof(kmer_start_end_t));
        ++marker_idx;
    }
    kmer_pos_vec_f.close();
    // write the total in the last entry
    kmer_pos_ptr_f.write((char*)&marker_idx, sizeof(uint64_t));
    kmer_pos_ptr_f.close();

    kmer_pos_vec.mmap_file(kmer_pos_vec_idx.c_str(), READ_WRITE_SHARED, 0, n_kmers);
    kmer_pos_ptr.mmap_file(kmer_pos_ptr_idx.c_str(), READ_WRITE_SHARED, 0, kmer_set.size()+1);

    //std::cerr << "querying kmers" << std::endl;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for (auto& hash : kmer_set) {
        uint64_t k = bphf->lookup(hash);
        uint64_t c = kmer_pos_ptr[k+1] - kmer_pos_ptr[k];
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto used_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    std::cerr << "done with " << kmer_set.size() << " @ " << (double)used_time/(double)kmer_set.size() << "ns/kmer" << std::endl;
        
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

void gyeet_index_t::save(const std::string& outfile) {
    //std::string bphf_out = outfile + ".bbh";
    std::ofstream f(outfile.c_str());
    bphf->save(f);
}

void gyeet_index_t::load(const std::string& outfile) {
}

void gyeet_index_t::for_values_of(const std::string& seq, const std::function<void(const kmer_start_end_t& v)>& lambda) {
    
}

}
