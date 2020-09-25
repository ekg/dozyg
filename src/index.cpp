#include "index.hpp"

namespace dozyg {

dozyg_index_t::dozyg_index_t(void) { }

dozyg_index_t::~dozyg_index_t(void) {
    if (loaded) {
        delete bphf;
        seq_fwd.munmap_file();
        seq_rev.munmap_file();
        edges.munmap_file();
        node_ref.munmap_file();
        kmer_pos_ref.munmap_file();
        kmer_pos_table.munmap_file();
    }
}

void dozyg_index_t::build(const HandleGraph& graph,
                          const uint64_t& _kmer_length,
                          const uint64_t& max_furcations,
                          const uint64_t& max_degree,
                          const double& sampling_rate,
                          const std::string& out_prefix) {

    // stash our node count
    n_nodes = graph.get_node_count();

    // and our kmer size
    kmer_length = _kmer_length;

    // calculate the sampling mod
    if (sampling_rate < 1.0) {
        sampling_mod = std::round(1.0 / sampling_rate);
    } else {
        assert(sampling_rate == 1.0);
        sampling_mod = 0;
    }

    // collect our sequences
    graph.for_each_handle(
        [&](const handle_t& h) {
            seq_length += graph.get_length(h);
        });

    // mark our node starts in fwd
    seq_bv = sdsl::bit_vector(seq_length+1); // allocates to 0

    // follow forward and reverse edges
    std::string seq_fwd_filename = out_prefix + ".sqf";
    std::ofstream seq_fwd_f(seq_fwd_filename.c_str(), std::ios::binary | std::ios::trunc);
    std::string seq_rev_filename = out_prefix + ".sqr";
    std::ofstream seq_rev_f(seq_rev_filename.c_str(), std::ios::binary | std::ios::trunc);
    std::string edge_filename = out_prefix + ".gye";
    std::ofstream edge_f(edge_filename.c_str(), std::ios::binary | std::ios::trunc);
    std::string node_ref_filename = out_prefix + ".gyn";
    std::ofstream node_ref_f(node_ref_filename.c_str(), std::ios::binary | std::ios::trunc);
    uint64_t seq_idx = 0;
    //uint64_t ref_idx = 0;
    graph.for_each_handle(
        [&](const handle_t& h) {
            seq_bv[seq_idx] = 1;
            const std::string seq = graph.get_sequence(h);
            seq_fwd_f << seq;
            seq_rev_f << seq;
            node_ref_t ref = { seq_idx, n_edges, 0 };
            graph.follow_edges(
                h, true,
                [&](const handle_t& p) {
                    edge_f.write(reinterpret_cast<char const*>(&p), sizeof(p));
                    ++ref.count_prev;
                });
            n_edges += ref.count_prev;
            node_ref_f.write(reinterpret_cast<char const*>(&ref), sizeof(ref));
            graph.follow_edges(
                h, false,
                [&](const handle_t& n) {
                    edge_f.write(reinterpret_cast<char const*>(&n), sizeof(n));
                    ++n_edges;
                });
            //ref_idx = n_edges;
            seq_idx += seq.length();
        });
    {
        // write a marker reference, to simplify counting of edges
        node_ref_t ref = { seq_idx, n_edges, 0 };
        node_ref_f.write(reinterpret_cast<char const*>(&ref), sizeof(ref));
    }
    assert(seq_idx == seq_length);
    seq_bv[seq_idx] = 1; // end marker
    // save our rank structure and bitvector
    sdsl::util::assign(seq_bv_rank, sdsl::bit_vector::rank_1_type(&seq_bv));
    std::string seq_bv_filename = out_prefix + ".sbv";
    std::ofstream seq_bv_f(seq_bv_filename.c_str(), std::ios::binary | std::ios::trunc);
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
    dna::reverse_complement_in_place(reinterpret_cast<char*>(&*seq_rev.begin()), seq_length);

    // open the graph vectors
    edges.mmap_file(edge_filename.c_str(), READ_WRITE_SHARED, 0, n_edges);
    node_ref.mmap_file(node_ref_filename.c_str(), READ_WRITE_SHARED, 0, n_nodes+1);

    // map from kmer hash idx to kmer start/end in graph sequence vector
    std::string kmer_pos_filename = out_prefix + ".kpos";
    std::ofstream kmer_pos_f(kmer_pos_filename.c_str(), std::ios::binary | std::ios::trunc);
    std::string kmer_set_filename = out_prefix + ".kset";
    std::ofstream kmer_set_f(kmer_set_filename.c_str(), std::ios::binary | std::ios::trunc);
        
    algorithms::for_each_kmer(
        graph, kmer_length, max_furcations, max_degree,
        [&](const kmer_t& kmer) {
            if (kmer.seq.find('N') == std::string::npos) {
                uint64_t hash = to_key(kmer.seq.c_str(), kmer.seq.length());
                //std::cerr << "kmer " << kmer.seq << " " << hash << std::endl;
                if (keep_key(hash)) {
                    //std::cerr << "recording!" << std::endl;
                    seq_pos_t begin_pos = get_seq_pos(kmer.begin.handle) + kmer.begin.pos;
                    seq_pos_t end_pos = get_seq_pos(kmer.end.handle) + kmer.end.pos;
                    kmer_pos_t p = { hash, begin_pos, end_pos };
#pragma omp critical (kmer_pos_write)
                    {
                        ++n_kmer_positions;
                        kmer_pos_f.write(reinterpret_cast<char*>(&p), sizeof(kmer_pos_t));
                    }
#pragma omp critical (kmer_set_write)
                    {
                        kmer_set_f.write(reinterpret_cast<char*>(&hash), sizeof(uint64_t));
                    }
                }
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
    // take uniques because our kmer generation seems to make a few non-unique kmer positions
    kmer_pos.erase(
        std::unique(
            kmer_pos.begin(), kmer_pos.end(),
            [](const kmer_pos_t& a,
               const kmer_pos_t& b) {
                return a.hash == b.hash
                    && a.begin == b.begin
                    && a.end == a.end;
            }), kmer_pos.end());
    n_kmer_positions = kmer_pos.size();

    // now we can iterate through our keys/values in order, storing them somewhere
    std::string kmer_pos_table_filename = out_prefix + ".kpv";
    std::ofstream kmer_pos_table_f(kmer_pos_table_filename.c_str(), std::ios::binary | std::ios::trunc);
    std::string kmer_pos_ref_filename = out_prefix + ".kpp";
    std::ofstream kmer_pos_ref_f(kmer_pos_ref_filename.c_str(), std::ios::binary | std::ios::trunc);

    uint64_t last_hash = std::numeric_limits<uint64_t>::max();
    uint64_t marker_idx = 0;
    for (uint64_t i = 0; i < n_kmer_positions; ++i) {
        auto& kp = kmer_pos[i];
        //std::cerr << "kmer_pos " << i << " hash is " << kp.hash << std::endl;
        if (kp.hash != last_hash) {
            //std::cerr << "kmer_pos_ref_f " << marker_idx << " " << kp.hash << std::endl;
            kmer_pos_ref_f.write(reinterpret_cast<char*>(&marker_idx), sizeof(uint64_t));
            last_hash = kp.hash;
        }
        kmer_start_end_t p = { kp.begin, kp.end };
        kmer_pos_table_f.write(reinterpret_cast<char*>(&p), sizeof(kmer_start_end_t));
        ++marker_idx;
    }
    kmer_pos_table_f.close();
    // write the total in the last entry
    kmer_pos_ref_f.write(reinterpret_cast<char*>(&marker_idx), sizeof(uint64_t));
    kmer_pos_ref_f.close();

    // remove the temporary file of kmer positions
    kmer_pos.munmap_file();
    std::remove(kmer_pos_filename.c_str());

    // write the bbhash index
    std::string bbhash_out_filename = out_prefix + ".bbx";
    std::ofstream bbhash_f(bbhash_out_filename.c_str(), std::ios::binary | std::ios::trunc);
    bphf->save(bbhash_f);
    bbhash_f.close();

    // write our metadata
    std::string metadata_filename = out_prefix + ".mtd";
    std::ofstream metadata_f(metadata_filename.c_str(), std::ios::binary | std::ios::trunc);
    metadata_f.write(reinterpret_cast<char*>(&seq_length), sizeof(seq_length));
    metadata_f.write(reinterpret_cast<char*>(&kmer_length), sizeof(kmer_length));
    metadata_f.write(reinterpret_cast<char*>(&sampling_mod), sizeof(sampling_mod));
    metadata_f.write(reinterpret_cast<char*>(&n_nodes), sizeof(n_nodes));
    metadata_f.write(reinterpret_cast<char*>(&n_edges), sizeof(n_edges));
    metadata_f.write(reinterpret_cast<char*>(&n_kmers), sizeof(n_kmers));
    metadata_f.write(reinterpret_cast<char*>(&n_kmer_positions), sizeof(n_kmer_positions));
    metadata_f.close();

    // todo
    // update the metadata
    // concatenate all the files together
    // remove the temporary files

}

void dozyg_index_t::load(const std::string& in_prefix) {

    // load our metadata
    std::string metadata_filename = in_prefix + ".mtd";
    std::ifstream metadata_f(metadata_filename.c_str(), std::ios::binary);
    metadata_f.read(reinterpret_cast<char*>(&seq_length), sizeof(seq_length));
    metadata_f.read(reinterpret_cast<char*>(&kmer_length), sizeof(kmer_length));
    metadata_f.read(reinterpret_cast<char*>(&sampling_mod), sizeof(sampling_mod));
    metadata_f.read(reinterpret_cast<char*>(&n_nodes), sizeof(n_nodes));
    metadata_f.read(reinterpret_cast<char*>(&n_edges), sizeof(n_edges));
    metadata_f.read(reinterpret_cast<char*>(&n_kmers), sizeof(n_kmers));
    metadata_f.read(reinterpret_cast<char*>(&n_kmer_positions), sizeof(n_kmer_positions));
    metadata_f.close();
    
    // load our sequences
    std::string seq_fwd_filename = in_prefix + ".sqf";
    std::string seq_rev_filename = in_prefix + ".sqr";
    seq_fwd.mmap_file(seq_fwd_filename.c_str(), READ_ONLY, 0, seq_length);
    seq_rev.mmap_file(seq_rev_filename.c_str(), READ_ONLY, 0, seq_length);

    // load our handle rank structurenn
    std::string seq_bv_filename = in_prefix + ".sbv";
    std::ifstream seq_bv_f(seq_bv_filename.c_str(), std::ios::binary);
    seq_bv.load(seq_bv_f);
    seq_bv_rank.load(seq_bv_f, &seq_bv);
    
    // load our graph topology
    std::string edge_filename = in_prefix + ".gye";
    edges.mmap_file(edge_filename.c_str(), READ_ONLY, 0, n_edges);
    std::string node_ref_filename = in_prefix + ".gyn";
    node_ref.mmap_file(node_ref_filename.c_str(), READ_ONLY, 0, n_nodes+1);

    // load our kmer table
    std::string kmer_pos_table_filename = in_prefix + ".kpv";
    kmer_pos_table.mmap_file(kmer_pos_table_filename.c_str(), READ_ONLY, 0, n_kmer_positions);
    std::string kmer_pos_ref_filename = in_prefix + ".kpp";
    kmer_pos_ref.mmap_file(kmer_pos_ref_filename.c_str(), READ_ONLY, 0, n_kmers+1);

    // load our bbhash
    bphf = new boomphf::mphf<uint64_t,hasher_t>();
    std::string bbx_in_filename = in_prefix + ".bbx";
    std::ifstream f(bbx_in_filename.c_str(), std::ios::binary);
    bphf->load(f);

    loaded = true;
}

// get a key representation of a sequence kmer
uint64_t dozyg_index_t::to_key(const char* seq, const size_t& len) const {
    return djb2_hash64(seq, len);
}

uint64_t dozyg_index_t::to_key(const std::string& seq) const {
    return to_key(seq.c_str(), seq.length());
}

// when sampling kmers, would this key pass our filter?
bool dozyg_index_t::keep_key(const uint64_t& key) const {
    return sampling_mod == 0 || key % sampling_mod == 0;
}

// iterate over values of a given key, if we would expect it in the index
void dozyg_index_t::for_values_of(const char* seq, const size_t& len, const std::function<void(const kmer_start_end_t& v)>& lambda) const {
    assert(len == kmer_length);
    uint64_t key = to_key(seq, len);
    if (keep_key(key)) {
        uint64_t idx = bphf->lookup(key);
        if (idx != std::numeric_limits<uint64_t>::max()) {
            uint64_t b = kmer_pos_ref[idx];
            uint64_t e = kmer_pos_ref[idx+1];
            for (uint64_t i = b; i < e; ++i) {
                lambda(kmer_pos_table[i]);
            }
        }
    }
}

// specialization for string input
void dozyg_index_t::for_values_of(const std::string& seq, const std::function<void(const kmer_start_end_t& v)>& lambda) const {
    for_values_of(seq.c_str(), seq.length(), lambda);
}

// node length in the graph
size_t dozyg_index_t::get_length(const handle_t& h) const {
    uint64_t i = handle_rank(h);
    return node_ref[i+1].seq_idx - node_ref[i].seq_idx;
}

// if a given handle is reverse (just use the set bit)
bool dozyg_index_t::is_reverse(const handle_t& h) const {
    return handle_is_rev(h);
}

// returns our sequence position type for the given strand
seq_pos_t dozyg_index_t::get_seq_pos(const handle_t& h) const {
    if (is_reverse(h)) {
        // flip our handle start onto the reverse strand
        return seq_pos::encode(seq_length
                               - node_ref[handle_rank(h)].seq_idx
                               - get_length(h),
                               is_reverse(h));
    } else {
        return seq_pos::encode(node_ref[handle_rank(h)].seq_idx,
                               is_reverse(h));
    }
}

handle_t dozyg_index_t::get_handle_at(const seq_pos_t& pos) const {
    bool is_rev = seq_pos::is_rev(pos);
    uint64_t offset = seq_pos::offset(pos);
    // the forward/reverse conversion is somewhat tricky because rank(N) yields the number of set bits in [0..N)
    if (is_rev) {
        return make_handle(seq_bv_rank(seq_length - offset) - 1, is_rev);
    } else {
        return make_handle(seq_bv_rank(offset + 1) - 1, is_rev);
    }
}

const char* dozyg_index_t::get_target(const seq_pos_t& pos) const {
    return seq_pos::is_rev(pos) ?
        &seq_rev[seq_pos::offset(pos)]
        : &seq_fwd[seq_pos::offset(pos)];
}

nid_t dozyg_index_t::get_id(const handle_t& h) const {
    return handle_rank(h) + 1;
}

void gyeet_index_t::follow_edges(
    const handle_t& h,
    bool go_left,
    const std::function<void(const handle_t&)>& func) const {
    uint64_t i = handle_rank(h);
    auto& node = node_ref[i];
    auto& next_edge_idx = node_ref[i+1].edge_idx;
    if (!handle_is_rev(h)) {
        if (go_left) {
            for (uint64_t j = node.edge_idx; j < node.edge_idx + node.count_prev; ++j) {
                func(edges[j]);
            }
        } else {
            for (uint64_t j = node.edge_idx + node.count_prev; j < next_edge_idx; ++j) {
                func(edges[j]);
            }
        }
    } else {
        if (!go_left) {
            for (uint64_t j = node.edge_idx; j < node.edge_idx + node.count_prev; ++j) {
                func(flip_handle(edges[j]));
            }
        } else {
            for (uint64_t j = node.edge_idx + node.count_prev; j < next_edge_idx; ++j) {
                func(flip_handle(edges[j]));
            }
        }
    }
}

}
