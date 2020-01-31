#include "subcommand.hpp"
#include <chrono>
#include "odgi.hpp"
#include "xg.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "index.hpp"
#include "chain.hpp"
#include "align.hpp"

namespace gyeet {

using namespace gyeet::subcommand;

int main_map(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "gyeet map";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("map sequences to a graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> idx_in_file(parser, "FILE", "load the index from this prefix", {'i', "index"});
    args::ValueFlag<std::string> query_seq(parser, "SEQ", "query one sequence", {'s', "one-sequence"});
    args::ValueFlag<uint64_t> max_gap_length(parser, "INT", "maximum gap length in chaining", {'G', "max-gap-length"});
    args::ValueFlag<double> max_mismatch_rate(parser, "FLOAT", "maximum allowed mismatch rate (default 0.1)", {'r', "max-mismatch-rate"});
    args::ValueFlag<uint64_t> threads(parser, "INT", "number of threads to use", {'t', "threads"});
    
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

    assert(argc > 0);

    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    }

    std::string idx_prefix;
    if (args::get(idx_in_file).empty()) {
        std::cerr << "[gyeet map] Error: an index basename is required (-i)" << std::endl;
    } else {
        idx_prefix = args::get(idx_in_file);
    }

    gyeet_index_t index;
    index.load(idx_prefix);

    const uint64_t& kmer_length = index.kmer_length;
    uint64_t max_gap = args::get(max_gap_length)
        ? args::get(max_gap_length) : 1000;

    double mismatch_rate = args::get(max_mismatch_rate)
        ? args::get(max_mismatch_rate)
        : 0.1;

    if (!args::get(query_seq).empty()) {
        const std::string& query = args::get(query_seq);
        auto anchors = anchors_for_query(index,
                                         query.c_str(),
                                         query.length());
        auto query_chains = chains(anchors,
                                   kmer_length,
                                   max_gap);
        auto query_superchains = superchains(query_chains);
        std::string query_name = "unknown";
        for (auto& superchain : query_superchains) {
            alignment_t aln = superalign(
                query_name,
                query.length(),
                query.c_str(),
                superchain,
                index,
                index.kmer_length,
                mismatch_rate);
            write_alignment_gaf(std::cout, aln, index);
        }
    }

    return 0;
}

static Subcommand gyeet_map("map", "map sequences to an index",
                            PIPELINE, 3, main_map);


}
