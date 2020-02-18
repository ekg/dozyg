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
        << std::min((int)std::round(chain.mapping_quality), 254) << "\t"
        << "ta:Z:chain" << std::endl;

}

void write_superchain_gaf(
    std::ostream& out,
    const superchain_t& superchain,
    const gyeet_index_t& index,
    const std::string& query_name,
    const uint64_t& query_length) {
    out << query_name << "\t"
        << query_length << "\t"
        << superchain.chains.front()->anchors.front()->query_begin << "\t"
        << superchain.chains.back()->anchors.back()->query_end << "\t"
        << "+" << "\t";  // we're always forward strand relative to our path
    uint64_t path_length = 0;
    for (auto& chain : superchain.chains) {
        for_handle_at_anchor_begin_in_chain(
            *chain, index,
            [&out,&index,&path_length](const handle_t& h) {
                path_length += index.get_length(h);
                out << (handle_is_rev(h) ? "<" : ">") << to_id(h);
            });
    }
    out << "\t"
        << path_length << "\t"
        << 0 << "\t" // fixme
        << path_length << "\t" // fixme
        << path_length << "\t" // fixme
        << path_length << "\t" // fixme
        << 0 << "\t"
        << "ta:Z:superchain" << std::endl;
        //<< std::min((int)std::round(superchain.mapping_quality), 254) << std::endl;

}

void write_alignment_gaf(
    std::ostream& out,
    const alignment_t& aln,
    const gyeet_index_t& index) {
    out << aln.query_name << "\t"
        << aln.query_length << "\t"
        << aln.query_begin << "\t"
        << aln.query_end << "\t"
        << "+" << "\t";  // we're always forward strand relative to our path
    uint64_t path_length = 0;
    for (auto& h : aln.path) {
        path_length += index.get_length(h);
        out << (handle_is_rev(h) ? "<" : ">") << to_id(h);
    }
    out << "\t"
        << path_length << "\t"
        << 0 << "\t" // fixme
        << path_length << "\t" // fixme
        << path_length << "\t" // fixme
        << path_length << "\t" // fixme
        << std::min((int)std::round(aln.mapping_quality), 254) << "\t"
        << "as:i:" << aln.score << "\t"
        << "ta:Z:" << (aln.is_secondary ? "secondary" : "primary") << "\t"
        << "cs:f:" << aln.chain_score << "\t"
        << "ac:i:" << aln.anchor_count << "\t"
        << "cg:Z:" << aln.cigar << std::endl;
}

// alignment
// assume that chains are over linearized segments of the graph
// apply global alignment within these regions to obtain sequence to graph sequence mapping
// use this to determine the correct graph path (and check it)
// re-score using the graph topology
alignment_t align(
    const std::string& query_name,
    const uint64_t& query_total_length,
    const char* query,
    const chain_t& chain,
    const gyeet_index_t& index,
    const uint64_t& extra_bp,
    const uint64_t& max_edit_distance,
    const seq_pos_t& query_pos,
    const seq_pos_t& query_length,
    const seq_pos_t& target_pos,
    const seq_pos_t& target_length) {

    //seq_pos_t query_pos = chain.query_begin(); //anchors.front()->query_begin;
    //seq_pos_t target_pos = chain.target_begin; // 0; // XXX TODO take this from the superchain chain_node_t target start and end
    /*
    if (seq_pos::offset(target_pos) >= extra_bp) {
        target_pos -= extra_bp;
    } else {
        target_pos = 0;
    }
    */
    //std::cerr << seq_pos::offset(target_pos) << std::endl;
    const char* query_begin = query + seq_pos::offset(query_pos);
    //uint64_t query_length = chain.query_end() - query_pos;
    const char* target_begin = index.get_target(target_pos);
    /*
    seq_pos_t target_end = seq_pos::encode(std::min(index.seq_length, seq_pos::offset(chain.anchors.back()->target_begin) + extra_bp),
                                           seq_pos::is_rev(chain.anchors.back()->target_begin));
    */
    //seq_pos_t target_end = chain.target_end; // anchors.back()->target_end;
    //uint64_t target_length = target_end - target_pos;
    /*
    if (target_end > target_pos) {
        target_length = target_end - target_pos;
    } else {
        // XXXX haxx
        target_length = query_length;
    }
    */
    //std::cerr << "query start " << seq_pos::offset(query_pos) << " length " << query_length << " target start " << seq_pos::to_string(target_pos) << " length " << target_length << std::endl;
//" last anchor begin " <<  seq_pos::to_string(chain.anchors.back()->target_begin) << " length " << target_length << std::endl;
    EdlibAlignResult result = edlibAlign(query_begin, query_length, target_begin, target_length,
                                         edlibNewAlignConfig(max_edit_distance, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));

    //std::cerr << "numlocs " << result.numLocations << std::endl;

    alignment_t aln;
    aln.query_name = query_name;
    aln.query_length = query_total_length;
    if (result.status != EDLIB_STATUS_OK) {
        std::cerr << "[gyeet map] alignment failure" << std::endl;
        assert(false);
        return aln;
    }

    graph_relativize(aln, query_pos, target_pos, index, result.alignment, result.alignmentLength, true);
    edlibFreeAlignResult(result);

    aln.anchor_count = chain.anchors.size();
    aln.is_secondary = chain.is_secondary;
    aln.chain_score = chain.score;
    aln.mapping_quality = chain.mapping_quality;
    // these are fine
    aln.query_begin = seq_pos::offset(query_pos);
    aln.query_end = seq_pos::offset(chain.query_end()); //anchors.back()->query_end);
    // change to offsets on the path we've given
    //aln.target_begin = seq_pos::offset(target_pos);
    //aln.target_end = seq_pos::offset(chain.anchors.back()->target_end);

    return aln;
}

bool has_matches(const cigar_t& cigar) {
    for (auto& elem : cigar) {
        switch (elem.second) {
        case 'M':
        case '=':
        case 'X':
            return true;
            break;
        default:
            break;
        }
    }
    return false;
}

uint64_t insertion_length(const cigar_t& cigar) {
    uint64_t ins_len = 0;
    for (auto& elem : cigar) {
        if (elem.second == 'I') {
            ++ins_len;
        }
    }
    return ins_len;
}

std::ostream& operator<<(std::ostream& out, const cigar_t& cigar) {
    for (auto& elem : cigar) {
        out << elem.first << elem.second;
    }
    return out;
}

void extend_cigar_string(std::string& cigar_str, const cigar_t& cigar) {
    for (auto& elem : cigar) {
        cigar_str.append(std::to_string(elem.first) + std::to_string(elem.second));
    }
}

/*
std::string cigar_string(const cigar_t& cigar) {
    std::string s;
    for (auto& elem : cigar) {
        s.append(std::to_string(elem.first) + std::to_string(elem.second));
    }
    return s;
}
*/

void extend_cigar(cigar_t& cigar, const cigar_t& extension) {
    if (cigar.empty()) {
        cigar = extension;
    } else if (extension.size()) {
        if (cigar.back().second == extension.front().second) {
            cigar.back().first += extension.front().first;
        } else {
            cigar.push_back(extension.front());
        }
        for (uint64_t i = 1; i < extension.size(); ++i) {
            cigar.push_back(extension[i]);
        }
    }
}

void extend_cigar(cigar_t& cigar, const uint32_t& len, const char& type) {
    if (cigar.empty()) {
        cigar.push_back(std::make_pair(len, type));
    } else if (len) {
        if (cigar.back().second == type) {
            cigar.back().first += len;
        } else {
            cigar.push_back(std::make_pair(len, type));
        }
    }
}

void cigar_handle_skips(cigar_t& cigar) {
    if (cigar.empty()) return;
    if (cigar.size() > 1
        && cigar.front().second == 'I'
        && cigar.at(1).second == 'D') {
        auto tmp = cigar.front();
        cigar.front() = cigar.at(1);
        cigar.at(1) = tmp;
    }
    // TODO in the case of short cigars, this could just swap the elements back
    if (cigar.size() > 1
        && cigar.back().second == 'I'
        && cigar.at(cigar.size()-1).second == 'D') {
        auto tmp = cigar.back();
        cigar.back() = cigar.at(cigar.size()-1);
        cigar.at(cigar.size()-1) = tmp;
    }
    // now we will start and end with deletions
    // we can set them to N if so
    if (cigar.front().second == 'D') {
        cigar.front().second = 'N';
    }
    if (cigar.back().second == 'D') {
        cigar.back().second = 'N';
    }
}

int64_t score_cigar(
    const cigar_t& cigar,
    int64_t match,
    int64_t mismatch,
    int64_t gap_open,
    int64_t gap_extend) {
    int64_t score = 0;
    uint64_t i = 0;
    for (auto& elem : cigar) {
        switch (elem.second) {
        case 'M':
        case '=':
            score += match * elem.first;
            break;
        case 'X':
            score -= mismatch * elem.first;
            break;
        case 'D':
        case 'I':
            // we're doing a global alignment
            // but we want a semiglobal one
            // don't penalize indels at the beginning or end
            if (i > 0 && i < cigar.size()-1) {
                score -= gap_open + gap_extend * elem.first;
            }
            break;
        default:
            break;
        }
        ++i;
    }
    return score;
}


void graph_relativize(
    alignment_t& aln,
    seq_pos_t query_pos,
    seq_pos_t target_pos,
    const gyeet_index_t& index,
    const unsigned char* const alignment,
    const int alignmentLength,
    const bool extended_cigar) {
    std::string alignment_cigar;
    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (!extended_cigar) {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }
    // project the alignment into the graph subset that it touches
    // computing:
    // implied path in graph
    // cigar relative to this path
    // target (in graph path) begin and end
    // and score / edit distance
    seq_pos_t target_first_match = std::numeric_limits<uint64_t>::max();
    seq_pos_t target_last_match = 0;

    auto record =
        [&aln](const handle_t& curr,
               cigar_t& curr_cigar) {
            // are there any positional matches in the last handle cigar?
            // if so, add them to our path and cigar
            if (has_matches(curr_cigar)) {
                aln.path.push_back(curr);
                extend_cigar(aln.cigar, curr_cigar);
            } else {
                // save insertions in the query that might have occurred between
                // matches and deletions relative to the graph vector
                uint64_t ins_len = insertion_length(curr_cigar);
                if (ins_len) {
                    curr_cigar.clear();
                    curr_cigar.push_back(std::make_pair(ins_len, 'I'));
                    extend_cigar(aln.cigar, curr_cigar);
                }
            }
        };
    //seq_pos_t query_pos = chain.anchors.front()->query_begin;
    //seq_pos_t target_pos = chain.anchors.front()->target_begin;
    handle_t curr = max_handle();
    cigar_t curr_cigar;
    for (uint64_t i = 0; i < alignmentLength; i++) {
        // determine if we need to update our path
        handle_t handle = index.get_handle_at(target_pos);
        if (i == 0) {
            curr = handle;
        } else if (handle != curr) {
            record(curr, curr_cigar);
            curr_cigar.clear();
            curr = handle;
        }
        /*
        if (i == alignmentLength) {
            break;
        }
        */

        // extend cigar
        const uint8_t& move_code = alignment[i];
        //std::cerr << "Move code " << moveCodeToChar[move_code] << std::endl;
        if (!curr_cigar.empty()
            && curr_cigar.back().second == moveCodeToChar[move_code]) {
            ++curr_cigar.back().first;
        } else {
            curr_cigar.push_back(std::make_pair(1, moveCodeToChar[move_code]));
        }

        // update our query/target pointers
        switch (move_code) {
        case 0:
        case 3:
            if (target_first_match == std::numeric_limits<uint64_t>::max()) {
                target_first_match = target_pos;
            }
            target_last_match = target_pos;
            ++query_pos;
            ++target_pos;
            break;
        case 1:
            ++query_pos;
            break;
        case 2:
            ++target_pos;
            break;
        default:
            assert(false);
            break;
        }
    }
    // last step
    record(curr, curr_cigar);
    // swap D with N at the beginning and end of the cigar
    // if the first element is 'I' and the second 'D', swap them
    cigar_handle_skips(aln.cigar);
    aln.score = score_cigar(aln.cigar);
    // record our start and end on the first and last handles
    if (!aln.path.empty()) {
        aln.first_handle_from_begin = target_first_match - index.get_seq_pos(aln.path.front());
        aln.last_handle_to_end = index.get_seq_pos(aln.path.back()) + index.get_length(aln.path.back()) - target_last_match - 1;
    }
}

uint64_t edit_distance_estimate(const chain_t& chain, const double& max_mismatch_rate) {
    if (chain.anchors.empty()) {
        return 0;
    } else {
        int64_t query_length = chain.query_end() - chain.query_begin();
        int64_t target_length = chain.target_end - chain.target_begin;
        //int64_t query_length = chain.anchors.back()->query_end - chain.anchors.front()->query_end;
        //int64_t target_length = chain.anchors.back()->target_end - chain.anchors.front()->target_end;
        return std::abs(query_length - target_length) + std::ceil(query_length * max_mismatch_rate);
    }
}

alignment_t superalign(
    const std::string& query_name,
    const uint64_t& query_total_length,
    const char* query,
    const superchain_t& superchain,
    const gyeet_index_t& index,
    const uint64_t& extra_bp,
    const double& max_mismatch_rate) {

    alignment_t superaln;
    superaln.query_name = query_name;
    superaln.query_length = query_total_length;
    superaln.query_begin = (superchain.chains.size() ? seq_pos::offset(superchain.chains.front()->query_begin()) : 0);
    superaln.query_end = (superchain.chains.size() ? seq_pos::offset(superchain.chains.back()->query_end()) : 0);
    superaln.mapping_quality = std::numeric_limits<double>::max();
    superaln.is_secondary = superchain.is_secondary;
    //uint64_t flanking_bp = std::ceil(index.kmer_length * (1 + max_mismatch_rate));
    for (uint64_t i = 0; i < superchain.chains.size(); ++i) {
        auto& chain = superchain.chains[i];
        // what's the gap from the last chain on the query?
        int64_t query_gap = (i > 0 ?
                             chain->query_begin() - superchain.chains[i-1]->query_end()
                             : chain->query_begin());
        //std::cerr << "query gap " << query_gap << std::endl;
        /*
        // tack on an insertion for unaligned intervening sequence
        if (query_gap < 0) {
            std::cerr << "[gyeet map] superchains overlap in query" << std::endl;
            assert(false);
            exit(1);
        } else if (query_gap > 0) {
            //extend_cigar(superaln.cigar, query_gap, 'I');
        }
        */
        if (-query_gap > chain->target_end - chain->target_begin) continue;
        seq_pos_t query_begin = chain->query_begin() - query_gap;
        assert(query_begin < chain->query_end());
        //std::cerr << "chain->target_begin " << seq_pos::to_string(chain->target_begin)  << std::endl;
        seq_pos_t target_begin = seq_pos::encode(std::max((int64_t)0, (int64_t)seq_pos::offset(chain->target_begin) - (int64_t)query_gap),
                                                 seq_pos::is_rev(chain->target_begin));
        //std::cerr << "target begin " << seq_pos::to_string(target_begin) << std::endl;
        if (query_begin >= chain->query_end()) continue;
        alignment_t aln
            = align(
                query_name,
                query_total_length,
                query,
                *chain,
                index,
                extra_bp,
                edit_distance_estimate(*chain, max_mismatch_rate),
                query_begin,
                chain->query_end() - query_begin,
                target_begin,
                chain->target_end - target_begin);
        // extend the superalignment
        // add deletions for distance to the target start
        if (aln.first_handle_from_begin > 0) {
            extend_cigar(superaln.cigar, aln.first_handle_from_begin, 'N');
        }
        // add the superalignment cigar
        extend_cigar(superaln.cigar, aln.cigar);
        //if (index.get_seq_pos(aln.path.back())
        if (aln.last_handle_to_end > 0) {
            extend_cigar(superaln.cigar, aln.last_handle_to_end, 'N');
        }
        // and path
        superaln.path.reserve(superaln.path.size() + aln.path.size());
        superaln.path.insert(superaln.path.end(), aln.path.begin(), aln.path.end());
        // and add to the score
        superaln.mapping_quality = std::min(superaln.mapping_quality, aln.mapping_quality); // hmm
        superaln.score += aln.score;
    }

    return superaln;
}

std::string alignment_cigar(
    const unsigned char* const alignment,
    const int alignmentLength,
    const bool extended_cigar) {
    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (!extended_cigar) {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }
    std::vector<char>* cigar = new std::vector<char>();
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
    // this copy should be optimized away via NRVO
    // but it might be even better to add things directly onto this string
    std::string s = std::string((const char*)&cigar->at(0));
    delete cigar;
    return s;
}



}
