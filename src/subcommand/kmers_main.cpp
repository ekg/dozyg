#include "subcommand.hpp"
#include "odgi.hpp"
#include "xg.hpp"
#include "algorithms/kmer.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/hash.hpp"
#include "pmhf.hpp"
#include <chrono>

namespace gyeet {

using namespace gyeet::subcommand;

int main_kmers(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "gyeet kmers";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("show and characterize the kmer space of the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> graph_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> idx_out_file(parser, "FILE", "save our kmer index to this file", {'o', "out"});
    args::ValueFlag<uint64_t> kmer_length(parser, "K", "the length of the kmers to generate", {'k', "kmer-length"});
    args::ValueFlag<uint64_t> max_furcations(parser, "N", "break at edges that would be induce this many furcations in a kmer", {'e', "max-furcations"});
    args::ValueFlag<uint64_t> max_degree(parser, "N", "remove nodes that have degree greater that this level", {'D', "max-degree"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});
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
        //ska::flat_hash_map<uint32_t, uint32_t> kmer_table;
        std::vector<uint64_t> kmers;
        uint64_t seen_kmers = 0;
        algorithms::for_each_kmer(graph, args::get(kmer_length), args::get(max_furcations), args::get(max_degree), [&](const kmer_t& kmer) {
                //int tid = omp_get_thread_num();
//#pragma omp atomic
                uint64_t hash = djb2_hash64(kmer.seq.c_str());
                //if (hash % 31 == 0) {
#pragma omp critical (kmers)
                kmers.push_back(hash);
                //}
                /*
#pragma omp critical (cerr)
                if (seen_kmers % 100000 == 0) {
                    std::cerr << seen_kmers << " " << kmer_table.size() << "\r";
                } else {
                    ++seen_kmers;
                }
                */
            });
        //std::cerr << std::endl;
        //ips4o::sort(kmers.begin(), kmers.end());
        std::sort(kmers.begin(), kmers.end());
        kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
        double gammaFactor = 4.0;
        boophf_t* bphf = new boomphf::mphf<uint64_t,hasher_t>(kmers.size(),
                                                              kmers,
                                                              get_thread_count(),
                                                              gammaFactor,
                                                              false,false);
        //kmers.clear();
        //std::cerr << "querying kmers" << std::endl;
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
        for (auto& hash : kmers) {
            bphf->lookup(hash);
        }
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        auto used_time = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        std::cerr << "done with " << kmers.size() << " @ " << (double)used_time/(double)kmers.size() << "ns/kmer" << std::endl;
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

static Subcommand gyeet_kmers("kmers", "process and dump the kmers of the graph",
                              PIPELINE, 3, main_kmers);


}
