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
        //std::cerr << std::string(seq+i, kmer_length) << std::endl;
        index.for_values_of(
            seq + i, kmer_length,
            [&](const kmer_start_end_t& v) {
                anchors.emplace_back(v.begin, v.end,
                                     seq_pos::encode(i, false),
                                     seq_pos::encode(i+kmer_length, false));
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
        anchor_i.max_chain_score = seed_length;
        int64_t min_j = bandwidth > i ? 0 : i - bandwidth;
        for (int64_t j = i; j >= min_j; --j) {
            anchor_t& anchor_j = anchors[j];
            double proposed_score = score_anchors(anchor_j,
                                                  anchor_i,
                                                  seed_length,
                                                  max_gap);
            std::cerr <<"prop score " <<proposed_score <<std::endl;
            if (proposed_score > seed_length
                && proposed_score > anchor_i.max_chain_score) {
                anchor_i.max_chain_score = proposed_score;
                anchor_i.best_predecessor = &anchor_j;
            }
        }
    }
    // collect chains
    std::vector<chain_t> chains;
    int64_t i = anchors.size()-1;
    while (i >= 0) {
        anchor_t* a = &anchors[i];
        std::cerr << a->best_predecessor << " "<< a->max_chain_score <<std::endl;
        if (a->best_predecessor != nullptr
            && a->max_chain_score > seed_length) { //!curr_chain) {
            chains.emplace_back();
            auto& curr_chain = chains.back();
            curr_chain.score = a->max_chain_score;
            do {
                curr_chain.anchors.push_back(a);
                anchor_t* b = a->best_predecessor;
                a->best_predecessor = nullptr; // mark done
                a = b; // swap
            } while (a->best_predecessor != nullptr);
            // investigate if a deque would be faster
            std::reverse(curr_chain.anchors.begin(),
                         curr_chain.anchors.end());
        }
        --i;
    }
    std::cerr << "chain count " <<chains.size() <<std::endl;
    // sort the chains by score, descending
    std::sort(chains.begin(), chains.end(),
              [](const chain_t& a,
                 const chain_t& b) {
                  return a.score > b.score;
              });
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
                if (!other_chain->processed()) { // not sure we need this check, based on definitions and order of evaluation
                    const seq_pos_t& other_begin = tree.start(idx);
                    const seq_pos_t& other_end = tree.end(idx);
                    uint64_t other_length = other_end - other_begin;
                    seq_pos_t ovlp_begin = std::max(chain_begin, other_begin);
                    seq_pos_t ovlp_end = std::min(chain_begin, other_begin);
                    uint64_t ovlp_length = ovlp_end - ovlp_begin;
                    if ((double)other_length < (double)ovlp_length * secondary_chain_threshold) {
                        // this chain is secondary
                        other_chain->mapping_quality = 0;
                        other_chain->is_secondary = true;
                        if (best_secondary->score < other_chain->score
                            || best_secondary == nullptr) {
                            best_secondary = other_chain;
                        }
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
        uint64_t query_length = b.query_end - a.query_end;
        uint64_t target_length = b.target_end - a.target_end;
        if (std::max(query_length, target_length) > max_gap) {
            return -std::numeric_limits<double>::max();
            //return std::numeric_limits<double>::min();
        } else {
            uint64_t gap_length = std::abs((int64_t)query_length - (int64_t)target_length);
            double gap_cost = gap_length == 0 ? 0
                : 0.01 * seed_length * gap_length + 0.5 * log2(gap_length);
            uint64_t match_length = std::min(std::min(query_length,
                                                      target_length),
                                             seed_length);
            // wtf here?
            return -( a.max_chain_score + match_length - gap_cost );
        }
    }
}

}
