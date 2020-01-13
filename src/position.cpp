#include "position.hpp"

namespace gyeet {

std::ostream& operator<<(std::ostream& out, const pos_t& pos) {
    return out << id(pos) << (is_rev(pos) ? "-" : "+") << offset(pos);
}

}
