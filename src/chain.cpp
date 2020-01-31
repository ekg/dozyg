#include "chain.hpp"

namespace gyeet {


std::vector<anchor_t>
anchors_for_query(
    const gyeet_index_t& index,
    const char* seq,
    const size_t& len) {
    std::vector<anchor_t> anchors;
    // for each kmer
    const uint64_t& kmer_length = index.kmer_length;
    for (uint64_t i = 0; i <= len-kmer_length; ++i) {
        //std::cerr << std::string(seq+i, kmer_length) << " " << i << " " << std::endl;
        index.for_values_of(
            seq + i, kmer_length,
            [&](const kmer_start_end_t& v) {
                /*
                std::cerr << std::string(seq+i, kmer_length) << " " << i << " "
                          << seq_pos::offset(v.begin) << " " << seq_pos::offset(v.end)
                          << std::endl;
                */
                anchors.emplace_back(seq_pos::encode(i, false),
                                     seq_pos::encode(i+kmer_length, false),
                                     v.begin, v.end);
            });
    }
    return anchors;
}

std::vector<chain_t>
chains(std::vector<anchor_t>& anchors,
       const uint64_t& seed_length,
       const uint64_t& max_gap,
       const uint64_t bandwidth,
       const double secondary_chain_threshold,
       const double max_mapq) {
    // sort the anchors by their ending position in the target
    std::sort(anchors.begin(), anchors.end(),
              [](const anchor_t& a,
                 const anchor_t& b) {
                  return a.target_end < b.target_end;
              });
    // calculate max chaining scores under our bandwidth
    // we walk from the last to first anchor
    for (int64_t i = 0; i < anchors.size(); ++i) {
        anchor_t& anchor_i = anchors[i];
        //std::cerr << "anchor_i " << anchor_i.query_begin << ".." << anchor_i.query_end << std::endl;
        anchor_i.max_chain_score = seed_length;
        int64_t min_j = bandwidth > i ? 0 : i - bandwidth;
        for (int64_t j = i-1; j >= min_j; --j) {
            anchor_t& anchor_j = anchors[j];
            //std::cerr << "anchor_j " << anchor_j.query_begin << ".." << anchor_j.query_end << std::endl;
            double proposed_score = score_anchors(anchor_j,
                                                  anchor_i,
                                                  seed_length,
                                                  max_gap);
            //std::cerr << "proposed " << proposed_score << " vs " << anchor_i.max_chain_score << std::endl;
            //std::cerr << "diff " << proposed_score - anchor_i.max_chain_score << std::endl;
            if (proposed_score > anchor_i.max_chain_score) {
                //std::cerr << "taken!" << std::endl;
                anchor_i.max_chain_score = proposed_score;
                anchor_i.best_predecessor = &anchor_j;
            }
        }
    }

    std::ofstream out("chains.dot");
    out << "digraph G {" << std::endl;
    for (auto& anchor : anchors) {
        out << "\"" << &anchor << "\" "
            << "[shape=box label=\"" << "score=" << anchor.max_chain_score << ","
            << "(" << seq_pos::to_string(anchor.query_begin) << ".."
            << seq_pos::to_string(anchor.query_end) << "),"
            << "(" << seq_pos::to_string(anchor.target_begin) << ".."
            << seq_pos::to_string(anchor.target_end) << ")\"];" << std::endl
            << "\"" << anchor.best_predecessor << "\" -> \"" << &anchor << "\";" << std::endl;
    }
    out << "}" << std::endl;
    out.close();

    // collect chains
    std::vector<chain_t> chains;
    int64_t i = anchors.size()-1;
    while (i >= 0) {
        anchor_t* a = &anchors[i];
        //std::cerr << "best predecessor " << a->best_predecessor << " "<< a->max_chain_score <<std::endl;
        if (a->best_predecessor != nullptr
            && a->max_chain_score > seed_length) { //!curr_chain) {
            chains.emplace_back();
            auto& curr_chain = chains.back();
            curr_chain.anchors.push_back(a);
            curr_chain.score = a->max_chain_score;
            do {
                anchor_t* b = a->best_predecessor;
                curr_chain.anchors.push_back(b);
                a->best_predecessor = nullptr; // mark done
                a = b; // swap
            } while (a->best_predecessor != nullptr);
            // investigate if a deque would be faster
            std::reverse(curr_chain.anchors.begin(),
                         curr_chain.anchors.end());
        }
        --i;
    }
    //std::cerr << "chain count " <<chains.size() <<std::endl;
    // sort the chains by score, descending
    std::sort(chains.begin(), chains.end(),
              [](const chain_t& a,
                 const chain_t& b) {
                  return a.score > b.score;
              });
    /*
    uint64_t j = 0;
    for (auto& chain : chains) {
        std::cerr << "chain " << ++j << " " << &chain << " " << chain.anchors.size() << " " << chain.score << " " << chain.mapping_quality << " " << chain.is_secondary << " " << chain.processed() << std::endl;
    }
    */
    // find the primary chains by examining their overlaps in the query
    IITree<seq_pos_t, chain_t*> tree;
    for (auto& chain : chains) {
        const seq_pos_t& query_begin = chain.anchors.front()->query_begin;
        const seq_pos_t& query_end = chain.anchors.back()->query_end;
        tree.add(query_begin, query_end, &chain);
    }
    tree.index();
    for (auto& chain : chains) {
        if (!chain.processed()) {
            const seq_pos_t& chain_begin = chain.anchors.front()->query_begin;
            const seq_pos_t& chain_end = chain.anchors.back()->query_end;
            uint64_t chain_length = chain_end - chain_begin;
            chain_t* best_secondary = nullptr;
            std::vector<size_t> ovlp;
            tree.overlap(chain_begin, chain_end, ovlp);
            for (auto& idx : ovlp) {
                chain_t* other_chain = tree.data(idx);
                if (other_chain != &chain && other_chain->score <= chain.score) {
                    const seq_pos_t& other_begin = tree.start(idx);
                    const seq_pos_t& other_end = tree.end(idx);
                    uint64_t other_length = other_end - other_begin;
                    seq_pos_t ovlp_begin = std::max(chain_begin, other_begin);
                    seq_pos_t ovlp_end = std::min(chain_end, other_end);
                    uint64_t ovlp_length = ovlp_end - ovlp_begin;
                    if ((double)ovlp_length > (double)other_length * secondary_chain_threshold) {
                        // this chain is secondary
                        other_chain->mapping_quality = 0;
                        other_chain->is_secondary = true;
                    }
                    if (best_secondary == nullptr
                        || best_secondary->score < other_chain->score) {
                        best_secondary = other_chain;
                    }
                }
            }
            if (best_secondary == nullptr) {
                chain.mapping_quality = max_mapq;
            } else {
                chain.mapping_quality =
                    40 * (1 - best_secondary->score / chain.score)
                    * std::min(1.0, (double)chain.anchors.size()/10.0)
                    * log(chain.score);
            }
        }
    }
    return chains;
}

double score_anchors(const anchor_t& a,
                     const anchor_t& b,
                     const uint64_t& seed_length,
                     const uint64_t& max_gap) {
    if (a.query_end >= b.query_end) {
        return -std::numeric_limits<double>::max();
    } else {
        uint64_t query_length = std::min(b.query_begin - a.query_begin,
                                         b.query_end - a.query_end);
        uint64_t query_overlap = (a.query_end > b.query_begin ? a.query_end - b.query_begin : 0);
        uint64_t target_length = std::min(b.target_begin - a.target_begin,
                                          b.target_end - a.target_end);
        //std::cerr << "query_length " << query_length << " target_length " << target_length << std::endl;
        if (std::max(query_length, target_length) > max_gap) {
            return -std::numeric_limits<double>::max();
            //return std::numeric_limits<double>::min();
        } else {
            //std::cerr << "query_length " << query_length << " target_length " << target_length << std::endl;
            uint64_t gap_length = std::abs((int64_t)query_length - (int64_t)target_length);
            double gap_cost = gap_length == 0 ? 0
                : 0.01 * seed_length * gap_length + 0.5 * log2(gap_length);
            uint64_t match_length = std::min(std::min(query_length,
                                                      target_length),
                                             seed_length);
            //std::cerr << "chain score is " << a.max_chain_score + match_length - gap_cost << std::endl;
            // round to 3 decimal digits, avoid problems with floating point instability causing chain truncation
            return std::round((a.max_chain_score + match_length - gap_cost) * 1000.0) / 1000.0 + query_overlap;
        }
    }
}

std::vector<superchain_t>
superchains(std::vector<chain_t>& chains,
            const uint64_t bandwidth) {
    // sort the chains by end coordinate in the query
    // score each with the previous N
    std::vector<chain_node_t> chain_nodes; chain_nodes.reserve(chains.size());
    // collect primary chains
    for (auto& chain : chains) {
        if (!chain.is_secondary) {
            chain_nodes.emplace_back(&chain);
        }
    }
    // sort by end in query
    std::sort(chain_nodes.begin(), chain_nodes.end(),
              [](const chain_node_t& a,
                 const chain_node_t& b) {
                  return a.chain->anchors.back()->query_end
                      < b.chain->anchors.back()->query_end;
              });
    // dynamic programming
    for (int64_t i = 0; i < chain_nodes.size(); ++i) {
        chain_node_t& chain_node_i = chain_nodes[i];
        //std::cerr << "anchor_i " << anchor_i.query_begin << ".." << anchor_i.query_end << std::endl;
        chain_node_i.max_superchain_score = chain_node_i.chain->score;
        int64_t min_j = bandwidth > i ? 0 : i - bandwidth;
        for (int64_t j = i-1; j >= min_j; --j) {
            chain_node_t& chain_node_j = chain_nodes[j];
            //anchor_t& anchor_j = anchors[j];
            //std::cerr << "anchor_j " << anchor_j.query_begin << ".." << anchor_j.query_end << std::endl;
            double proposed_score = score_chain_nodes(chain_node_j, chain_node_i);
            //std::cerr << "proposed " << proposed_score << " vs " << chain_node_i.max_superchain_score << std::endl;
            //std::cerr << "diff " << proposed_score - chain_node_i.max_superchain_score << std::endl;
            if (proposed_score > chain_node_i.max_superchain_score) {
                //std::cerr << "taken!" << std::endl;
                chain_node_i.max_superchain_score = proposed_score;
                chain_node_i.best_predecessor = &chain_node_j;
            }
        }
    }
    // collect superchains
    std::vector<superchain_t> superchains;
    int64_t i = chain_nodes.size()-1;
    while (i >= 0) {
        chain_node_t* a = &chain_nodes[i];
        //std::cerr << "best predecessor " << a->best_predecessor << " "<< a->max_chain_score <<std::endl;
        if (a->best_predecessor != nullptr
            && a->max_superchain_score) { // HMMMM
            //std::cerr << "adding superchain" << std::endl;
            superchains.emplace_back();
            auto& curr_superchain = superchains.back();
            curr_superchain.chains.push_back(a->chain);
            a->used = true;
            curr_superchain.score = a->max_superchain_score;
            do {
                chain_node_t* b = a->best_predecessor;
                curr_superchain.chains.push_back(b->chain);
                b->used = true;
                a->best_predecessor = nullptr; // mark done
                a = b; // swap
            } while (a->best_predecessor != nullptr);
            // investigate if a deque would be faster
            std::reverse(curr_superchain.chains.begin(),
                         curr_superchain.chains.end());
        }
        --i;
    }
    // add back in all chains that aren't part of a superchain
    for (auto& chain_node : chain_nodes) {
        if (!chain_node.used) {
            superchains.emplace_back();
            superchains.back().chains.push_back(chain_node.chain);
            superchains.back().score = chain_node.chain->score;
        }
    }
    //std::cerr << "chain count " <<chains.size() <<std::endl;
    //std::cerr << "superchain count " <<superchains.size() <<std::endl;
    // sort the superchains by score, descending
    std::sort(superchains.begin(), superchains.end(),
              [](const superchain_t& a,
                 const superchain_t& b) {
                  return a.score > b.score;
              });
    return superchains;
}

double score_chain_nodes(const chain_node_t& a,
                         const chain_node_t& b) {
    // do we overlap? --> 0
    if (a.chain->anchors.back()->query_end
        > b.chain->anchors.front()->query_begin) {
        return 0;
    } else {
        // otherwise, add scores
        return a.max_superchain_score + b.chain->score;
    }
}

}
