#include "index.hpp"

namespace gyeet {

void gyeet_index_t::build(const HandleGraph& graph,
                          const uint64_t& kmer_length,
                          const uint64_t& max_furcations,
                          const uint64_t& max_degree,
                          const std::string& out_prefix) {

    // stash our node count
    n_nodes = graph.get_node_count();

    // collect our sequences
    graph.for_each_handle(
        [&](const handle_t& h) {
            seq_length += graph.get_length(h);
        });

    // mark our node starts in fwd
    seq_bv.resize(seq_length);

    // follow forward and reverse edges
    std::string seq_fwd_filename = out_prefix + ".sqf";
    std::ofstream seq_fwd_f(seq_fwd_filename.c_str());
    std::string seq_rev_filename = out_prefix + ".sqr";
    std::ofstream seq_rev_f(seq_rev_filename.c_str());
    std::string edge_filename = out_prefix + ".gye";
    std::ofstream edge_f(edge_filename.c_str());
    std::string node_ref_filename = out_prefix + ".gyn";
    std::ofstream node_ref_f(node_ref_filename.c_str());
    uint64_t seq_idx = 0;
    uint64_t ref_idx = 0;
    graph.for_each_handle(
        [&](const handle_t& h) {
            seq_bv[seq_idx] = 1;
            const std::string seq = graph.get_sequence(h);
            seq_fwd_f << seq;
            seq_rev_f << seq;
            node_ref_t ref = { seq_idx, ref_idx++, 0 };
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
            seq_idx += seq.length();
        });
    {
        // write a marker reference, to simplify counting of edges
        node_ref_t ref = { seq_idx, ref_idx, 0 };
        node_ref_f.write((char*)&ref, sizeof(ref));
    }

    // save our rank structure and bitvector
    sdsl::util::assign(seq_bv_rank, sdsl::bit_vector::rank_1_type(&seq_bv));
    std::string seq_bv_filename = out_prefix + ".sbv";
    std::ofstream seq_bv_f(seq_bv_filename.c_str());
    seq_bv.serialize(seq_bv_f);
    seq_bv_rank.serialize(seq_bv_f);
    seq_bv_f.close();

    // close our graph files
    seq_fwd_f.close();
    seq_rev_f.close();
    edge_f.close();
    node_ref_f.close();

    // load our sequence vectors
    seq_fwd.mmap_file(seq_fwd_filename.c_str(), READ_WRITE_SHARED, 0, seq_length);
    seq_rev.mmap_file(seq_rev_filename.c_str(), READ_WRITE_SHARED, 0, seq_length);

    // and reverse complement the revcomp one
    dna::reverse_complement_in_place((char*)&*seq_rev.begin(), seq_length);

    // open the graph vectors
    edges.mmap_file(edge_filename.c_str(), READ_WRITE_SHARED, 0, n_edges);
    node_ref.mmap_file(node_ref_filename.c_str(), READ_WRITE_SHARED, 0, n_nodes+1);

    // map from kmer hash idx to kmer start/end in graph sequence vector
    std::string kmer_pos_filename = out_prefix + ".kpos";
    std::ofstream kmer_pos_f(kmer_pos_filename.c_str());
    std::string kmer_set_filename = out_prefix + ".kset";
    std::ofstream kmer_set_f(kmer_set_filename.c_str());
        
    algorithms::for_each_kmer(
        graph, kmer_length, max_furcations, max_degree,
        [&](const kmer_t& kmer) {
            uint64_t hash = djb2_hash64(kmer.seq.c_str());
            seq_pos_t begin_pos = get_seq_pos(kmer.begin.handle) + kmer.begin.pos;
            seq_pos_t end_pos = get_seq_pos(kmer.end.handle) + kmer.end.pos;
            kmer_pos_t p = { hash, begin_pos, end_pos };
#pragma omp critical (kmer_pos_write)
            {
                ++n_kmer_positions;
                kmer_pos_f.write((char*)&p, sizeof(kmer_pos_t));
            }
#pragma omp critical (kmer_set_write)
            {
                kmer_set_f.write((char*)&hash, sizeof(uint64_t));
            }
        });
    kmer_pos_f.close();
    kmer_set_f.close();

    mmappable_vector<uint64_t> kmer_set;
    kmer_set.mmap_file(kmer_set_filename.c_str(), READ_WRITE_SHARED, 0, n_kmer_positions);
    ips4o::parallel::sort(kmer_set.begin(), kmer_set.end());
    //std::cerr << "total kmers " << kmer_set.size() << std::endl;
    kmer_set.erase(std::unique(kmer_set.begin(), kmer_set.end()), kmer_set.end());
    n_kmers = kmer_set.size();
    //std::cerr << "unique kmers " << kmer_set.size() << std::endl;

    double gammaFactor = 8.0;
    bphf = new boomphf::mphf<uint64_t,hasher_t>(n_kmers,
                                                kmer_set,
                                                get_thread_count(),
                                                gammaFactor,
                                                false,false);

    // remove the kmer set working file
    kmer_set.munmap_file();
    std::remove(kmer_set_filename.c_str());

    // rename our kmers
    // build our sequence space index
    mmappable_vector<kmer_pos_t> kmer_pos;
    kmer_pos.mmap_file(kmer_pos_filename.c_str(), READ_WRITE_SHARED, 0, n_kmer_positions);
    // set our kmer hashes
#pragma omp parallel for
    for (uint64_t i = 0; i < n_kmer_positions; ++i) {
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
    std::string kmer_pos_table_filename = out_prefix + ".kpv";
    std::ofstream kmer_pos_table_f(kmer_pos_table_filename.c_str());
    std::string kmer_pos_ref_filename = out_prefix + ".kpp";
    std::ofstream kmer_pos_ref_f(kmer_pos_ref_filename.c_str());

    uint64_t last_hash = std::numeric_limits<uint64_t>::max();
    uint64_t marker_idx = 0;
    for (uint64_t i = 0; i < n_kmer_positions; ++i) {
        auto& kp = kmer_pos[i];
        if (kp.hash != last_hash) {
            //std::cerr << "kmer_pos_ref_f " << marker_idx << " " << kp.hash << std::endl;
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
    std::remove(kmer_pos_filename.c_str());

    // write the bbhash index
    std::string bbhash_out_filename = out_prefix + ".bbx";
    std::ofstream f(bbhash_out_filename.c_str(), std::ios::binary);
    bphf->save(f);

}

void gyeet_index_t::load(const std::string& in_prefix) {

    // load our sequences
    std::string seq_fwd_filename = in_prefix + ".sqf";
    std::string seq_rev_filename = in_prefix + ".sqr";
    seq_length = filesize(seq_fwd_filename.c_str());
    seq_fwd.mmap_file(seq_fwd_filename.c_str(), READ_ONLY, 0, seq_length);
    seq_rev.mmap_file(seq_rev_filename.c_str(), READ_ONLY, 0, seq_length);

    // load our handle rank structure
    std::string seq_bv_filename = in_prefix + ".sbv";
    std::ifstream seq_bv_f(seq_bv_filename.c_str());
    seq_bv.load(seq_bv_f);
    seq_bv_rank.load(seq_bv_f, &seq_bv);
    
    // load our graph topology
    std::string edge_filename = in_prefix + ".gye";
    n_edges = filesize(edge_filename.c_str()) / sizeof(handle_t);
    edges.mmap_file(edge_filename.c_str(), READ_ONLY, 0, n_edges);
    std::string node_ref_filename = in_prefix + ".gyn";
    node_ref.mmap_file(node_ref_filename.c_str(), READ_ONLY, 0, n_nodes+1);

    // load our kmer table
    std::string kmer_pos_table_filename = in_prefix + ".kpv";
    n_kmer_positions = filesize(kmer_pos_table_filename.c_str()) / sizeof(kmer_start_end_t);
    kmer_pos_table.mmap_file(kmer_pos_table_filename.c_str(), READ_ONLY, 0, n_kmer_positions);
    std::string kmer_pos_ref_filename = in_prefix + ".kpp";
    n_kmers = filesize(kmer_pos_ref_filename.c_str())-1;
    kmer_pos_ref.mmap_file(kmer_pos_ref_filename.c_str(), READ_ONLY, 0, n_kmers+1);

    // load our bbhash
    bphf = new boomphf::mphf<uint64_t,hasher_t>();
    std::string bbx_in_filename = in_prefix + ".bbx";
    std::ifstream f(bbx_in_filename.c_str(), std::ios::binary);
    bphf->load(f);
}

void gyeet_index_t::for_values_of(const std::string& seq, const std::function<void(const kmer_start_end_t& v)>& lambda) {
    uint64_t hash = djb2_hash64(seq.c_str());
    uint64_t k = bphf->lookup(hash);
    uint64_t b = kmer_pos_ref[k];
    uint64_t e = kmer_pos_ref[k+1];
    std::cerr << "hash " << hash << " " << k << " " << b << " " << e << std::endl;
    for (uint64_t i = b; i < e; ++i) {
        lambda(kmer_pos_table[i]);
    }
}

size_t gyeet_index_t::get_length(const handle_t& h) {
    uint64_t i = handle_rank(h);
    return node_ref[i+1].seq_idx - node_ref[i].seq_idx;
}

bool gyeet_index_t::is_reverse(const handle_t& h) {
    return handle_is_rev(h);
}

seq_pos_t gyeet_index_t::get_seq_pos(const handle_t& h) {
    return seq_pos::encode(node_ref[handle_rank(h)].seq_idx, is_reverse(h));
}

}
