#pragma once
#ifndef _SPRESS_UTILS_H_
#define _SPRESS_UTILS_H_

#include <cstdint>
#include <hilbert.hpp>
#include <vector>

#include "spress-omp.hpp"

namespace spress {
namespace utils {

template<typename numeric_t, typename pos_t>
inline pos_t coord_to_dim(pos_t n, numeric_t v, numeric_t max, numeric_t min)
{
    if (v >= min && v < max) {
        return static_cast<pos_t>(
          floor((v - min) * (static_cast<numeric_t>(n) / (max - min)))
        );
    } else if (v == max) {
        return n - 1;
    } else {
        return static_cast<pos_t>(-1);
    }
}

template<typename numeric_t, typename pos_t>
inline numeric_t dim_to_coord(pos_t n, pos_t v, numeric_t max, numeric_t min)
{
    return static_cast<numeric_t>(
      min + ((v + 0.5) * ((max - min) / static_cast<numeric_t>(n)))
    );
}

} // namespace utils
} // namespace spress

#endif