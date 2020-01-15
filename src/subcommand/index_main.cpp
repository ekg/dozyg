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
#include "mmap_allocator.h"
#include "mmap_exception.h"
#include "mmappable_vector.h"
#include "utility.hpp"

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
        algorithms::for_each_kmer(graph, args::get(kmer_length), args::get(max_furcations), args::get(max_degree), [&](const kmer_t& kmer) {
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
        algorithms::for_each_kmer(graph, args::get(kmer_length), args::get(max_furcations), args::get(max_degree), [&](const kmer_t& kmer) {
                uint64_t hash = djb2_hash64(kmer.seq.c_str());
                //if (hash % 31 == 0) {
//#pragma omp critical (kmers)
                kmer_set.append(hash);
                //}
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
