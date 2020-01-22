#include "subcommand.hpp"
#include <chrono>
#include "odgi.hpp"
#include "xg.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "index.hpp"

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
    args::ValueFlag<std::string> idx_out_file(parser, "FILE", "save our index with this prefix (defaults to input file name and path)", {'o', "out"});
    args::ValueFlag<uint64_t> kmer_length(parser, "K", "the length of the kmers to generate", {'k', "kmer-length"});
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

    //odgi::graph_t graph;
    xg::XG graph;
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

    std::string idx_prefix;
    if (args::get(idx_out_file).empty()) {
        if (infile == "-") {
            std::cerr << "[gyeet index] Error: reading graph from stdin but no output file specified with -o" << std::endl;
            return 4;
        } else {
            idx_prefix = infile;
        }
    } else {
        idx_prefix = args::get(idx_out_file);
    }

    gyeet_index_t index;
    index.build(graph,
                args::get(kmer_length),
                args::get(max_furcations),
                args::get(max_degree),
                idx_prefix);

    return 0;
}

static Subcommand gyeet_index("index", "process and dump the kmers of the graph",
                              PIPELINE, 3, main_index);


}
