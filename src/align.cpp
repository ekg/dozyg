#include "align.hpp"

namespace gyeet {

void for_handle_at_anchor_begin_in_chain(
    const chain_t& chain,
    const gyeet_index_t& index,
    const std::function<void(const handle_t&)>& func) {
    for (auto& anchor : chain.anchors) {
        func(index.get_handle_at(anchor->target_begin));
    }
}

}
