#include "chain.hpp"

namespace gyeet {

std::vector<chain_t>
chains(std::vector<anchor_t>& anchors,
       const uint64_t& seed_length,
       const uint64_t& max_gap,
       const uint64_t bandwidth) {
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
        uint64_t min_j = bandwidth > i ? 0 : i - bandwidth;
        for (int64_t j = i; j >= min_j; --j) {
            const anchor_t& anchor_j = anchors[j];
            anchor_i.max_chain_score
                = std::max(
                    anchor_i.max_chain_score,
                    score_anchors(anchor_j,
                                  anchor_i,
                                  seed_length,
                                  max_gap));
        }
    }
    // collect chains
    std::vector<chain_t> chains;
    return chains;
}

double score_anchors(const anchor_t& a,
                     const anchor_t& b,
                     const uint64_t& seed_length,
                     const uint64_t& max_gap) {
    if (a.query_end >= b.query_end) {
        return std::numeric_limits<double>::min();
    } else {
        uint64_t query_length = b.query_end - a.query_end;
        uint64_t target_length = b.target_end - a.target_end;
        if (std::max(query_length, target_length) > max_gap) {
            return std::numeric_limits<double>::min();
        } else {
            uint64_t gap_length = std::abs((int64_t)query_length - (int64_t)target_length);
            double gap_cost = gap_length == 0 ? 0
                : 0.01 * seed_length * gap_length + 0.5 * log2(gap_length);
            uint64_t match_length = std::min(std::min(query_length,
                                                      target_length),
                                             seed_length);
            return a.max_chain_score + match_length - gap_cost;
        }
    }
}

}
