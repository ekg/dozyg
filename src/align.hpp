#pragma once

#include <iostream>
#include "edlib.h"
#include "index.hpp"
#include "chain.hpp"

namespace gyeet {

void for_handle_at_anchor_begin_in_chain(
    const chain_t& chain,
    const gyeet_index_t& index,
    const std::function<void(const handle_t&)>& func);


void write_chain_gaf(
    std::ostream& out,
    const chain_t& chain,
    const gyeet_index_t& index,
    const std::string& query_name,
    const uint64_t& query_length);

}
