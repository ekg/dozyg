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

void write_chain_gaf(
    std::ostream& out,
    const chain_t& chain,
    const gyeet_index_t& index,
    const std::string& query_name,
    const uint64_t& query_length) {
    out << query_name << "\t"
        << query_length << "\t"
        << chain.anchors.front()->query_begin << "\t"
        << chain.anchors.back()->query_end << "\t"
        << "+" << "\t";  // we're always forward strand relative to our path
    uint64_t path_length = 0;
    for_handle_at_anchor_begin_in_chain(
        chain, index,
        [&out,&index,&path_length](const handle_t& h) {
            path_length += index.get_length(h);
            out << (handle_is_rev(h) ? "<" : ">") << to_id(h);
        });
    out << "\t"
        << path_length << "\t"
        << 0 << "\t" // fixme
        << path_length << "\t" // fixme
        << path_length << "\t" // fixme
        << path_length << "\t" // fixme
        << std::min((int)std::round(chain.mapping_quality), 254) << std::endl;

}

// alignment
// assume that chains are over linearized segments of the graph
// apply global alignment within these regions to obtain sequence to graph sequence mapping
// use this to determine the correct graph path (and check it)
// re-score using the graph topology
alignment_t align(
    const char* query,
    const chain_t& chain,
    const gyeet_index_t& index) {

    const char* query_begin = query + seq_pos::offset(chain.anchors.front()->query_begin);
    uint64_t query_length = chain.anchors.back()->query_end - chain.anchors.front()->query_begin;
    const char* target_begin = index.get_target(chain.anchors.front()->target_begin);
    uint64_t target_length = chain.anchors.back()->target_end - chain.anchors.front()->target_begin;
    EdlibAlignResult result = edlibAlign(query_begin, query_length, target_begin, target_length,
                                         edlibNewAlignConfig(query_length, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
        printf("%d\n", result.editDistance);
    }
    std::cerr << result.alignmentLength << std::endl;
    char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
    printf("%s\n", cigar);
    free(cigar);
    edlibFreeAlignResult(result);
}

void fill_alignment_path(
    alignment_t& aln,
    const gyeet_index_t& index,
    const unsigned char* const alignment,
    const int alignmentLength) {
    
}

void fill_alignment_cigar(
    alignment_t& aln,
    const unsigned char* const alignment,
    const int alignmentLength,
    const bool& extended_cigar) {
    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (!extended_cigar) {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }
    vector<char>* cigar = new vector<char>();
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                cigar->push_back('0' + numOfSameMoves % 10);
                numDigits++;
            }
            reverse(cigar->end() - numDigits, cigar->end());
            // Write code of move to cigar string.
            cigar->push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    assert(false);
                    delete cigar;
                    return 0;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
    cigar->push_back(0);  // Null character termination.
    aln.cigar = string(cigar); // save in our alignment
    delete cigar;
}



}
