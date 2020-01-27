#include "align.hpp"

namespace gyeet {

void for_handle_at_anchor_begin_in_chain(
    const chain_t& chain,
    const gyeet_index_t& index,
    const std::function<void(const handle_t&)>& func) {
    handle_t last = max_handle();
    for (auto& anchor : chain.anchors) {
        handle_t curr = index.get_handle_at(anchor->target_begin);
        if (curr != last) {
            func(curr);
            last = curr;
        }
    }
}

/*
GAF format
https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf

Col     Type    Description
1       string  Query sequence name
2       int     Query sequence length
3       int     Query start (0-based; closed)
4       int     Query end (0-based; open)
5       char    Strand relative to the path: "+" or "-"
6       string  Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
7       int     Path length
8       int     Start position on the path (0-based)
9       int     End position on the path (0-based)
10      int     Number of residue matches
11      int     Alignment block length
12      int     Mapping quality (0-255; 255 for missing)
*/

void write_gaf_record(
    std::ostream& out,
    const chain_t& chain,
    const gyeet_index_t& index) {
}


}
