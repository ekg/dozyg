#include "subcommand.hpp"
#include <chrono>
#include "odgi.hpp"
#include "xg.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "index.hpp"

namespace dozyg {

using namespace dozyg::subcommand;

int main_query(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "dozyg query";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("query a pmhf-based kmer positional index of the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> graph_in_file(parser, "FILE", "load the graph from this file", {'i', "in"});
    args::ValueFlag<std::string> idx_in_file(parser, "FILE", "load the index from this prefix", {'x', "index"});
    //args::ValueFlag<uint64_t> kmer_length(parser, "K", "the length of the kmers to use for queries", {'k', "kmer-length"});
    args::ValueFlag<uint64_t> max_furcations(parser, "N", "break at edges that would be induce this many furcations in a kmer", {'e', "max-furcations"});
    args::ValueFlag<uint64_t> max_degree(parser, "N", "remove nodes that have degree greater that this level", {'D', "max-degree"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});
    //args::ValueFlag<std::string> work_prefix(parser, "FILE", "write temporary files with this prefix", {'w', "work-prefix"});
    //args::Flag kmers_stdout(parser, "", "write the kmers to stdout", {'c', "stdout"});

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
    assert(argc > 0);
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

    std::string idx_prefix;
    if (args::get(idx_in_file).empty()) {
        if (infile == "-") {
            std::cerr << "[dozyg query] Error: reading graph from stdin but no input file specified with -o" << std::endl;
            return 4;
        } else {
            idx_prefix = infile;
        }
    } else {
        idx_prefix = args::get(idx_in_file);
    }

    dozyg_index_t index;
    index.load(idx_prefix);

    algorithms::for_each_kmer(
        graph,
        index.kmer_length,
        args::get(max_furcations),
        args::get(max_degree),
        [&](const kmer_t& kmer) {
            if (kmer.seq.find('N') == std::string::npos) {
                //std::cerr << "kmer " << kmer.seq << std::endl;
                seq_pos_t begin_pos = index.get_seq_pos(kmer.begin.handle) + kmer.begin.pos;
                seq_pos_t end_pos = index.get_seq_pos(kmer.end.handle) + kmer.end.pos;
                //std::cerr << "begin/end " << seq_pos::to_string(begin_pos) << " " << seq_pos::to_string(end_pos) << std::endl;
                // check if we can find it in the index
                if (index.keep_key(index.to_key(kmer.seq))) {
                    bool found = false;
                    index.for_values_of(
                        kmer.seq,
                        [&](const kmer_start_end_t& p) {
                            //std::cerr << "--> begin/end " << seq_pos::to_string(p.begin) << " " << seq_pos::to_string(p.end) << std::endl;
                            if (p.begin == begin_pos && p.end == end_pos) {
                                //std::cerr << "found!" << std::endl;
                                found = true;
                            }
                        });
                    assert(found);
                }
            }
        });

    return 0;
}

static Subcommand dozyg_query("query", "test a dozyg index",
                              PIPELINE, 3, main_query);


}
