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
    args::ValueFlag<uint64_t> max_gap_length(parser, "N", "maximum gap length in chaining", {'G', "max-gap-length"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});
    
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

    if (!args::get(query_seq).empty()) {
        const std::string& query = args::get(query_seq);
        auto anchors = anchors_for_query(index,
                                         query.c_str(),
                                         query.length());
        auto query_chains = chains(anchors,
                                   kmer_length,
                                   max_gap);
        for (auto& chain : query_chains) {
            std::cout << chain.score << " "
                      << chain.mapping_quality << " "
                      << chain.anchors.size() << " "
                      << chain.is_secondary
                      << std::endl;
            nid_t last = 0;
            for (auto& a : chain.anchors) {
                std::cout << "," << a->query_begin << "/" << seq_pos::to_string(a->target_begin);
            }
            std::cout << std::endl;
            for_handle_in_chain(
                chain, index,
                [&last](const handle_t& h) {
                    nid_t curr = to_id(h);
                    if (curr != last) {
                        std::cout << (handle_is_rev(h) ? "<" : ">") << curr;
                        last = curr;
                    }
                });
        }
    }

    return 0;
}

static Subcommand gyeet_map("map", "map sequences to an index",
                            PIPELINE, 3, main_map);


}
