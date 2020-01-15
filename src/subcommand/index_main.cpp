#include "subcommand.hpp"
#include <chrono>
#include "odgi.hpp"
#include "xg.hpp"
#include "algorithms/kmer.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/hash.hpp"
#include "pmhf.hpp"
#include "mmmultiset.hpp"
#include "mmmultimap.hpp"
#include "mmap_allocator.h"
#include "mmap_exception.h"
#include "mmappable_vector.h"
#include "utility.hpp"
#include "ips4o.hpp"

namespace gyeet {

using namespace gyeet::subcommand;

int main_index(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "gyeet index";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("build a bbhash based kmer index of the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> graph_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> idx_out_file(parser, "FILE", "save our kmer index to this file", {'o', "out"});
    args::ValueFlag<uint64_t> kmer_length(parser, "K", "the length of the kmers to generate", {'k', "kmer-length"});
    args::ValueFlag<uint64_t> max_furcations(parser, "N", "break at edges that would be induce this many furcations in a kmer", {'e', "max-furcations"});
    args::ValueFlag<uint64_t> max_degree(parser, "N", "remove nodes that have degree greater that this level", {'D', "max-degree"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});
    args::ValueFlag<std::string> work_prefix(parser, "FILE", "write temporary files with this prefix", {'w', "work-prefix"});
    args::Flag kmers_stdout(parser, "", "write the kmers to stdout", {'c', "stdout"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    odgi::graph_t graph;
    //xg::XG graph;
    assert(argc > 0);
    assert(args::get(kmer_length));
    std::string infile = args::get(graph_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }
    }
    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    }
    std::vector<std::vector<kmer_t>> buffers(get_thread_count());
    if (args::get(kmers_stdout)) {
        algorithms::for_each_kmer(
            graph, args::get(kmer_length), args::get(max_furcations), args::get(max_degree),
            [&](const kmer_t& kmer) {
                int tid = omp_get_thread_num();
                auto& buffer = buffers.at(tid);
                buffer.push_back(kmer);
                if (buffer.size() > 1e5) {
#pragma omp critical (cout)
                    {
                        for (auto& kmer : buffer) {
                            std::cout << kmer << "\n";
                        }
                        buffer.clear();
                    }
                }
            });
        for (auto& buffer : buffers) {
            for (auto& kmer : buffer) {
                std::cout << kmer << "\n";
            }
            buffer.clear();
        }
        std::cout.flush();
    } else {
        std::string kmer_set_idx = args::get(work_prefix) + ".kmer_set";
        mmmulti::set<uint64_t> kmer_set(kmer_set_idx);
        //ska::flat_hash_map<uint32_t, uint32_t> kmer_table;
        //std::vector<uint64_t> kmers;
        uint64_t seen_kmers = 0;
        algorithms::for_each_kmer(
            graph, args::get(kmer_length), args::get(max_furcations), args::get(max_degree),
            [&](const kmer_t& kmer) {
                uint64_t hash = djb2_hash64(kmer.seq.c_str());
                kmer_set.append(hash);
            });
        //std::cerr << std::endl;
        //ips4o::sort(kmers.begin(), kmers.end());

        std::string kmer_vec_idx = args::get(work_prefix) + ".kmer_vec";
        uint64_t n_kmers = 0;
        kmer_set.index();
        kmer_set.for_each_unique_value(
            [&](const uint64_t& k) {
                ++n_kmers;
            });
        //nodes_.mmap_file(filename.c_str(), READ_WRITE_SHARED, 0, record_count())
        mmap_allocator_namespace::mmappable_vector<uint64_t> kmers_vec; //, mmap_allocator<uint64_t>>
        allocate_file(kmer_vec_idx, n_kmers * sizeof(uint64_t));
        kmers_vec.mmap_file(kmer_vec_idx.c_str(), mmap_allocator_namespace::READ_WRITE_SHARED, 0, n_kmers);
        //kmers_vec.
        uint64_t i = 0;
        kmer_set.for_each_unique_value(
            [&](const uint64_t& k) {
                kmers_vec[i++] = k;
            });
        //std::remove(kmer_set_idx.c_str());
        //std::sort(kmers.begin(), kmers.end());
        //kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
        double gammaFactor = 4.0;
        boophf_t* bphf = new boomphf::mphf<uint64_t,hasher_t>(kmers_vec.size(),
                                                              kmers_vec,
                                                              get_thread_count(),
                                                              gammaFactor,
                                                              false,false);

        //kmers.clear();

        // build our sequence space index
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

        // the way we store kmers in our index
        struct kmer_pos_t {
            uint64_t kmer_hash;
            int64_t begin;
            int64_t end;
        };

        // map from kmer hash idx to kmer start/end in graph sequence vector
        //mmmulti::map<uint64_t, std::pair<int64_t, int64_t>> kmer_map(kmer_map_idx);
        std::string kmer_pos_idx = args::get(work_prefix) + ".kmer_pos";
        ofstream kmer_pos_f(kmer_pos_idx.c_str());
        uint64_t n_positions = 0;
        algorithms::for_each_kmer(
            graph, args::get(kmer_length), args::get(max_furcations), args::get(max_degree),
            [&](const kmer_t& kmer) {
                uint64_t hash = djb2_hash64(kmer.seq.c_str());
                //kmer_set.append(hash);
                uint64_t bbidx = bphf->lookup(hash);
                auto& b = kmer.begin;
                int64_t start_pos = seq_bv_rank(id(b)) + (is_rev(b) ? graph.get_length(graph.get_handle(id(b))) - offset(b) : offset(b));
                if (is_rev(b)) start_pos = -start_pos;
                auto& e = kmer.end;
                int64_t end_pos = seq_bv_rank(id(e)) + (is_rev(e) ? graph.get_length(graph.get_handle(id(e))) - offset(e) : offset(e));
                if (is_rev(e)) end_pos = -end_pos;
                auto v = std::make_pair(start_pos, end_pos);
                kmer_pos_t p = { bbidx, start_pos, end_pos };
                //kmer_map.add(bbidx, std::make_pair(start_pos, end_pos));
#pragma omp critical (kmer_pos_write)
                {
                    ++n_positions;
                    kmer_pos_f.write((char*)&p, sizeof(kmer_pos_t));
                }
            });

        kmer_pos_f.close();
        mmap_allocator_namespace::mmappable_vector<kmer_pos_t> kmers_pos; //, mmap_allocator<uint64_t>>;
        kmers_pos.mmap_file(kmer_pos_idx.c_str(), mmap_allocator_namespace::READ_WRITE_SHARED, 0, n_positions);
        ips4o::parallel::sort(
            kmers_pos.begin(), kmers_pos.end(),
            [](const kmer_pos_t& a, const kmer_pos_t& b) {
                return a.kmer_hash < b.kmer_hash
                                     || a.kmer_hash == b.kmer_hash
                                     && (a.begin < b.begin
                                         || a.begin == b.begin && a.end < b.end);
            });

        // now we can iterate through our keys/values in order, storing them somewhere
        /*
        for (auto& p : kmers_pos) {
            std::cerr << p.kmer_hash << " " << p.begin << " " << p.end << std::endl;
        }
        */
        /*
        //std::string kmer_pos_idx = args::get(work_prefix) + ".kmer_pos";
        //allocate_file(kmer_pos_idx, n_kmers * sizeof(uint64_t));
        */
        // TODO here we gotta build our postion index

        //std::cerr << "querying kmers" << std::endl;
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
        for (auto& hash : kmers_vec) {
            bphf->lookup(hash);
        }
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        auto used_time = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        std::cerr << "done with " << kmers_vec.size() << " @ " << (double)used_time/(double)kmers_vec.size() << "ns/kmer" << std::endl;
        
        /*
        algorithms::for_each_kmer(graph, args::get(kmer_length), [&](const kmer_t& kmer) {
                uint64_t hash = djb2_hash64(kmer.seq.c_str());
                bphf->lookup(hash);
#pragma omp atomic
                ++seen_kmers;
            });
        std::cerr << "done with " << seen_kmers << " kmers" << std::endl;
        */
        if (!args::get(idx_out_file).empty()) {
            ofstream f(args::get(idx_out_file).c_str());
            bphf->save(f);
        }
        delete bphf;
    }
    return 0;
}

static Subcommand gyeet_index("index", "process and dump the kmers of the graph",
                              PIPELINE, 3, main_index);


}
